#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_NORMALIZATION } from './modules/qc_normalization.nf'
include { PCA_SELECTION }    from './modules/pca_selection.nf'
include { CLUSTERING_UMAP }  from './modules/clustering_umap.nf'

workflow {
    input_data  = file(params.input_path)
    config_file = file(params.config)

    qc_out    = QC_NORMALIZATION(input_data, config_file)
    pca_out   = PCA_SELECTION(qc_out.rds, config_file)
    CLUSTERING_UMAP(pca_out.rds, config_file)
}
