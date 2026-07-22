#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC_NORMALIZATION } from './modules/qc_normalization.nf'

workflow {
    input_data = file(params.input_path)
    config_file = file(params.config)
    QC_NORMALIZATION(input_data, config_file)
}
