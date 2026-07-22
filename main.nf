#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_NORMALIZATION } from './modules/qc_normalization.nf'
include { PCA_SELECTION }    from './modules/pca_selection.nf'
include { CLUSTERING_UMAP }  from './modules/clustering_umap.nf'
include { FIND_MARKERS }     from './modules/find_markers.nf'
include { ANNOTATION }       from './modules/annotation.nf'
include { TRAJECTORIES }     from './modules/trajectories.nf'

workflow {
    input_data  = file(params.input_path)
    config_file = file(params.config)

    analysis_cfg = new org.yaml.snakeyaml.Yaml().load(file(params.config).text)

    qc_out      = QC_NORMALIZATION(input_data, config_file)
    pca_out     = PCA_SELECTION(qc_out.rds, config_file)
    clust_out   = CLUSTERING_UMAP(pca_out.rds, config_file)
    markers_out = FIND_MARKERS(clust_out.rds, config_file)

    if (analysis_cfg.annotation?.perform) {
        markers_db = file("${projectDir}/" + analysis_cfg.annotation.markers_file)
        annot_out  = ANNOTATION(clust_out.rds, config_file, markers_db)

        // Trajectoires : conditionnelles, alimentees par l'annotation
        if (analysis_cfg.trajectories?.perform) {
            TRAJECTORIES(annot_out.rds, config_file)
        }
    }
}
