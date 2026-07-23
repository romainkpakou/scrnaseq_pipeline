#!/usr/bin/env nextflow

// =============================================================================
// SCRNASEQ_PIPELINE : WORKFLOW PRINCIPAL
// =============================================================================
//
// Pipeline Nextflow DSL2 pour l'analyse single-cell RNA-seq, de la matrice
// brute 10X Genomics au rapport HTML final.
//
// ARCHITECTURE
//   Chaque etape analytique est encapsulee dans un module independant appelant
//   un script R parametrable. Les modules communiquent par des canaux : la
//   sortie de l'un devient l'entree du suivant, ce qui etablit les dependances
//   d'execution sans qu'aucun ordonnancement ne soit code explicitement.
//
// CONFIGURATION A DEUX NIVEAUX
//   Le pipeline distingue deux espaces de parametres qu'il importe de ne pas
//   confondre :
//     - les parametres Nextflow (params.*), issus de nextflow.config et de la
//       ligne de commande, qui gouvernent l'execution : chemins d'entree et de
//       sortie, profil d'execution
//     - la configuration d'analyse (analysis_cfg), lue depuis le fichier YAML,
//       qui gouverne la science : seuils, resolutions, methodes
//   Les scripts R ne connaissent que la seconde ; le workflow lit les deux.
//
// USAGE
//   nextflow run main.nf \
//       --config params/example_pbmc.yaml \
//       --input_path data/raw/pbmc3k/hg19 \
//       --output_path results
//
// AUTEUR
//   Romain KPAKOU
// =============================================================================

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// IMPORT DES MODULES
// -----------------------------------------------------------------------------

include { QC_NORMALIZATION } from './modules/qc_normalization.nf'
include { PCA_SELECTION }    from './modules/pca_selection.nf'
include { CLUSTERING_UMAP }  from './modules/clustering_umap.nf'
include { FIND_MARKERS }     from './modules/find_markers.nf'
include { ANNOTATION }       from './modules/annotation.nf'
include { TRAJECTORIES }     from './modules/trajectories.nf'
include { GENERATE_REPORT }  from './modules/generate_report.nf'


workflow {

    // -------------------------------------------------------------------------
    // 1. RESOLUTION DES ENTREES
    // -------------------------------------------------------------------------

    input_data  = file(params.input_path)
    config_file = file(params.config)

    // Lecture de la configuration d'analyse cote Nextflow. La bibliotheque
    // snakeyaml est embarquee dans la distribution Nextflow, aucune dependance
    // supplementaire n'est donc requise.
    //
    // Cette lecture est indispensable : les conditions d'execution des modules
    // optionnels dependent du contenu du YAML, que params.* ne connait pas.
    //
    // Contrainte de syntaxe : cette instruction doit imperativement se trouver
    // a l'interieur du bloc workflow. Placee au niveau du script, elle serait
    // rejetee a la compilation par les versions recentes de Nextflow, qui
    // n'admettent que des declarations a ce niveau.
    analysis_cfg = new org.yaml.snakeyaml.Yaml().load(file(params.config).text)

    // -------------------------------------------------------------------------
    // 2. DECLARATION DES SCRIPTS R
    // -------------------------------------------------------------------------
    // Les scripts sont declares comme fichiers et transmis en entree de chaque
    // process, plutot que references par un chemin absolu dans la commande.
    //
    // Cette approche est ce qui rend le cache fiable. Nextflow calcule
    // l'empreinte d'une tache a partir de ses entrees declarees : un script
    // simplement reference dans la chaine de commande verrait ses modifications
    // ignorees, et une reprise avec -resume reutiliserait silencieusement des
    // resultats perimes. En le declarant comme entree, toute modification du
    // code invalide automatiquement le cache du module concerne et de tous ceux
    // qui en dependent.

    bin_qc      = file("${projectDir}/bin/01_qc_normalization.R")
    bin_pca     = file("${projectDir}/bin/02_pca_selection.R")
    bin_clust   = file("${projectDir}/bin/03_clustering_umap.R")
    bin_markers = file("${projectDir}/bin/04_find_markers.R")
    bin_annot   = file("${projectDir}/bin/05_annotation.R")
    bin_traj    = file("${projectDir}/bin/06_trajectories.R")
    bin_report  = file("${projectDir}/bin/07_generate_report.R")

    // -------------------------------------------------------------------------
    // 3. MODULES OBLIGATOIRES
    // -------------------------------------------------------------------------
    // Ces quatre etapes constituent le socle de toute analyse single-cell et
    // s'executent systematiquement. Le chainage est explicite : chaque module
    // recoit l'objet produit par le precedent, ce qui suffit a Nextflow pour
    // etablir l'ordre d'execution.

    qc_out      = QC_NORMALIZATION(input_data, config_file, bin_qc)
    pca_out     = PCA_SELECTION(qc_out.rds, config_file, bin_pca)
    clust_out   = CLUSTERING_UMAP(pca_out.rds, config_file, bin_clust)
    markers_out = FIND_MARKERS(clust_out.rds, config_file, bin_markers)

    // -------------------------------------------------------------------------
    // 4. AGREGATION DES SORTIES POUR LE RAPPORT
    // -------------------------------------------------------------------------
    // Trois canaux collecteurs accumulent progressivement les sorties de tous
    // les modules. L'operateur mix fusionne plusieurs canaux en un seul, sans
    // garantie d'ordre, ce qui convient ici puisque le rapport identifie les
    // fichiers par leur nom et non par leur position.
    //
    // Ces canaux sont enrichis au fil des blocs conditionnels ci-dessous, de
    // sorte que le rapport ne recoive que les sorties effectivement produites.

    all_figures = qc_out.figures.mix(pca_out.figures, clust_out.figures,
                                     markers_out.figures)
    all_tables  = markers_out.markers.mix(markers_out.top_markers)
    all_logs    = qc_out.log.mix(pca_out.log, clust_out.log, markers_out.log)

    // -------------------------------------------------------------------------
    // 5. ANNOTATION CELLULAIRE (CONDITIONNELLE)
    // -------------------------------------------------------------------------
    // L'operateur de navigation sure ?. protege contre l'absence de la section
    // dans le YAML : la condition est evaluee a faux sans lever d'erreur, ce
    // qui rend le pipeline tolerant a une configuration minimale.

    if (analysis_cfg.annotation?.perform) {

        // Le chemin de la base de marqueurs est relatif dans le YAML, ce qui le
        // rend lisible et portable. Il est resolu ici en chemin absolu depuis
        // la racine du projet avant d'etre transmis au process.
        markers_db = file("${projectDir}/" + analysis_cfg.annotation.markers_file)

        annot_out = ANNOTATION(clust_out.rds, config_file, markers_db, bin_annot)

        all_figures = all_figures.mix(annot_out.figures)
        all_tables  = all_tables.mix(annot_out.summary)
        all_logs    = all_logs.mix(annot_out.log)

        // ---------------------------------------------------------------------
        // 6. TRAJECTOIRES CELLULAIRES (CONDITIONNELLES)
        // ---------------------------------------------------------------------
        // Ce bloc est imbrique dans le precedent, et non place a son suite :
        // l'analyse de trajectoire selectionne les populations par leur nom
        // biologique, information que seule l'annotation produit. Cette
        // imbrication traduit donc une dependance reelle et non un choix de
        // presentation.

        if (analysis_cfg.trajectories?.perform) {

            traj_out = TRAJECTORIES(annot_out.rds, config_file, bin_traj)

            all_figures = all_figures.mix(traj_out.figures)
            all_tables  = all_tables.mix(traj_out.dynamic, traj_out.summary)
            all_logs    = all_logs.mix(traj_out.log)
        }
    }

    // -------------------------------------------------------------------------
    // 7. RAPPORT HTML FINAL
    // -------------------------------------------------------------------------
    // L'operateur collect transforme chaque canal, qui emet ses fichiers un a
    // un, en une emission unique contenant l'ensemble. Cette conversion a deux
    // effets : elle contraint le rapport a attendre l'achevement de tous les
    // modules amont, et elle garantit que tous les fichiers sont deposes
    // ensemble dans un meme repertoire de travail, condition necessaire pour
    // que le template les y retrouve.

    if (analysis_cfg.report) {

        template = file("${projectDir}/assets/report_template.Rmd")

        GENERATE_REPORT(all_figures.collect(),
                        all_tables.collect(),
                        all_logs.collect(),
                        config_file,
                        template,
                        bin_report)
    }
}
