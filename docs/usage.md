# Guide d'utilisation

## Prerequis

| Composant | Version minimale | Verification |
|---|---|---|
| Java | 17 | `java -version` |
| Nextflow | 22.04 | `nextflow -version` |
| R | 4.3 | `R --version` |
| pandoc | 1.12.3 | `pandoc --version` |

Packages R requis : Seurat (>= 5.0), slingshot, SingleCellExperiment,
DelayedMatrixStats, optparse, yaml, dplyr, ggplot2, patchwork, rmarkdown, knitr.

L'environnement Conda fourni installe l'ensemble de ces dependances :

    conda env create -f environment/environment.yml

## Execution

Lancement standard :

    nextflow run main.nf \
        --config params/example_pbmc.yaml \
        --input_path data/raw/pbmc3k/hg19 \
        --output_path results

Avec environnement Conda :

    nextflow run main.nf -profile conda \
        --config params/example_pbmc.yaml \
        --input_path data/raw/pbmc3k/hg19

Sur cluster SLURM :

    nextflow run main.nf -profile slurm,conda \
        --config params/example_pbmc.yaml \
        --input_path /chemin/vers/donnees

Reprise apres interruption. L'option -resume reutilise les resultats des etapes
deja calculees. Les scripts R etant declares comme entrees des process, toute
modification de code invalide automatiquement le cache du module concerne.

    nextflow run main.nf -resume --config params/example_pbmc.yaml --input_path ...

## Donnees d'entree

Le pipeline attend un repertoire au format 10X Genomics contenant :

- matrix.mtx : matrice de comptage au format creux
- barcodes.tsv : identifiants des cellules
- genes.tsv ou features.tsv : identifiants des genes, selon la version de
  Cell Ranger

## Configuration de l'analyse

Tous les parametres scientifiques sont regroupes dans un unique fichier YAML.

Pour analyser un nouveau jeu de donnees, copier le fichier d'exemple et
l'adapter :

    cp params/example_pbmc.yaml params/mon_analyse.yaml

Les points a ajuster en priorite sont l'espece, les seuils de controle qualite,
la resolution de clustering, et le chemin vers la base de marqueurs canoniques
correspondant au tissu etudie.

## Sorties

    results/
    |-- report.html              rapport final consolide
    |-- figures/                 figures au format PNG, 300 dpi
    |-- objects/                 objets intermediaires (.rds)
    |-- tables/                  tables de resultats (.csv)
    |-- logs/                    journaux d'execution par module
    +-- pipeline_info/           trace, rapport d'execution, chronologie

## Modules optionnels

Deux etapes peuvent etre desactivees depuis la configuration :

    annotation:
      perform: false      # desactive l'annotation et, par dependance, les trajectoires

    trajectories:
      perform: false      # desactive la seule analyse de trajectoire

L'analyse de trajectoire depend de l'annotation, car elle selectionne les
populations d'interet par leur nom biologique. Desactiver l'annotation desactive
donc necessairement les trajectoires.

## Resolution des problemes courants

pandoc version 1.12.3 or higher is required
  pandoc est absent du PATH. Une session RStudio le fournit implicitement, ce
  qui n'est pas le cas d'une execution en ligne de commande.
  Correctif : sudo apt install pandoc, ou utiliser le profil conda.

aucun package nomme 'slingshot'
  Dependance Bioconductor manquante.
  Correctif : BiocManager::install(c("slingshot", "DelayedMatrixStats"))

Un cluster annote Unknown
  Aucun type de la base n'atteint le seuil de confiance pour ce cluster.
  Consulter les scores detailles dans results/logs/annotation_log.txt. Deux
  causes possibles : la base de marqueurs ne couvre pas ce type cellulaire, ou
  le cluster correspond a un etat inattendu meritant une inspection manuelle.
