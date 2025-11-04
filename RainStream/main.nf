#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process METHYL_MERGE_DATA {
    tag ""

    publishDir "${workflow.projectDir}/data", mode: 'copy'

    input:
    path clinical_csv
    path methylation_csv

    output:
    path "methyl_data_trim.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/methyl_data_processing.py \
        --clinical ${clinical_csv} \
        --methylation ${methylation_csv} \
        --output .
    """
}


process METHYL_SUMMARIZE_REPS {
    tag ""

    publishDir "${workflow.projectDir}/data", mode: 'copy'

    input:
    path methyl_trim_csv

    output:
    path "methyl_data_trim_condensed.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/methyl_rep_handling.py \
        --input ${methyl_trim_csv} \
        --output .
    """
}


process METHYL_WIDE_FORMAT {
    tag ""

    publishDir "${workflow.projectDir}/data", mode: 'copy'

    input:
    path methyl_condensed_csv

    output:
    path "methyl_data_wide.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/methyl_wide_formatting.py \
        --input ${methyl_condensed_csv} \
        --output .
    """
}


process METHYL_FEATURE_SELECTION {
    tag ""

    publishDir "${workflow.projectDir}/output", mode: 'copy'

    input:
    path methyl_wide_csv

    output:
    path "feature_importance_scores.csv"
    path "top_markers_plot.png", optional: true

    script:
    """
    python3 ${workflow.projectDir}/src/methyl_feature_selection.py \
        --input ${methyl_wide_csv} \
        --output . \
        --plot
    """
}

workflow {
    // Channels pointing to the input CSVs
    clinical_ch = Channel.fromPath('data/clinical.csv')
    methylation_ch = Channel.fromPath('data/methylation.csv')

    // Run the preprocessing process
    methyl_data_trim = METHYL_MERGE_DATA(
        clinical_ch,
        methylation_ch
    )

    // Run the summarize replicates process
    methyl_data_condensed = METHYL_SUMMARIZE_REPS(methyl_data_trim)

    // Run the wide format conversion process
    methyl_data_wide = METHYL_WIDE_FORMAT(methyl_data_condensed)
    
    // Add feature selection
    METHYL_FEATURE_SELECTION(methyl_data_wide)
}