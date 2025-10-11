#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MERGE_DATA {
    tag ""

    input:
    path clinical_csv
    path methylation_csv
    val trimmed_flag

    output:
    path "merged_data*.csv"

    publishDir "${workflow.projectDir}/data", mode: 'copy'

    script:
    """
    python3 ${workflow.projectDir}/src/data_processing.py \
        --clinical ${clinical_csv} \
        --methylation ${methylation_csv} \
        --output . \
        ${trimmed_flag ? "--trimmed" : ""}
    """
}
workflow {
    // Channels pointing to the input CSVs
    clinical_ch = Channel.fromPath('data/clinical.csv')
    methylation_ch = Channel.fromPath('data/methylation.csv')

    // Run the preprocessing process
    MERGE_DATA(
        clinical_ch,
        methylation_ch,
        params.trimmed
    )
}