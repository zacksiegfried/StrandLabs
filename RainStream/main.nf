#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MERGE_DATA {
    tag ""

    publishDir "${workflow.projectDir}/data", mode: 'copy'

    input:
    path clinical_csv
    path methylation_csv
    val trimmed_flag

    output:
    path "merged_data*.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/data_processing.py \
        --clinical ${clinical_csv} \
        --methylation ${methylation_csv} \
        --output . \
        ${trimmed_flag ? "--trimmed" : ""}
    """
}

process SUMMARIZE_REPS {
    tag ""

    publishDir "${workflow.projectDir}/data", mode: 'copy'

    input:
    path merged_trim_csv

    output:
    path "merged_data_trim_condensed.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/rep_handling.py \
        --input ${merged_trim_csv} \
        --output .
    """
}

workflow {
    // Channels pointing to the input CSVs
    clinical_ch = Channel.fromPath('data/clinical.csv')
    methylation_ch = Channel.fromPath('data/methylation.csv')

    // Run the preprocessing process
    merged_data_trim = MERGE_DATA(
        clinical_ch,
        methylation_ch,
        params.trimmed
    )

    // Run the summarize replicates process
    SUMMARIZE_REPS(merged_data_trim)
}