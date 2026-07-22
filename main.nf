#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_NORMALIZATION } from './modules/qc_normalization.nf'
include { PCA_SELECTION }    from './modules/pca_selection.nf'

workflow {
    input_data  = file(params.input_path)
    config_file = file(params.config)

    // Etape 1 : QC
    qc_out = QC_NORMALIZATION(input_data, config_file)

    // Etape 2 : ACP, alimentee par la sortie .rds du QC
    PCA_SELECTION(qc_out.rds, config_file)
}
