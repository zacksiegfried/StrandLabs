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


process METHYL_PRECISION_PROFILE {
    tag ""

    publishDir "${workflow.projectDir}/output", mode: 'copy'

    input:
    path methyl_trim_csv

    output:
    path "methyl_precision_profile_data.csv"
    path "marker_precision_profiles/*.png", optional: true

    script:
    """
    python3 ${workflow.projectDir}/src/methyl_precision_profile.py \
        --input ${methyl_trim_csv} \
        --output . \
        --plot
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


process METHYL_CANCER_DETECTION {
    tag ""

    publishDir "${workflow.projectDir}/output", mode: 'copy'

    input:
    path methyl_wide_csv

    output:
    path "cancer_detection_predictions.csv"
    path "cancer_detection_model_summary.csv"
    path "cancer_detection_roc_curve.png", optional: true

    script:
    """
    python3 ${workflow.projectDir}/src/methyl_cancer_detection.py \
        --input ${methyl_wide_csv} \
        --output . \
        --plot
    """
}


process METHYL_CANCER_TYPE {
    tag ""

    publishDir "${workflow.projectDir}/output", mode: 'copy'

    input:
    path methyl_wide_csv

    output:
    path "cancer_type_predictions.csv"
    path "cancer_type_model_summary.csv"
    path "cancer_type_confusion_matrix.png", optional: true

    script:
    """
    python3 ${workflow.projectDir}/src/cancer_type_classification.py \
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

    // Assess precision across analytical measuring range (uses replicate-level data)
    METHYL_PRECISION_PROFILE(methyl_data_trim)

    // Run the summarize replicates process
    methyl_data_condensed = METHYL_SUMMARIZE_REPS(methyl_data_trim)

    // Run the wide format conversion process
    methyl_data_wide = METHYL_WIDE_FORMAT(methyl_data_condensed)
    
    // Rank markers by importance against cancer_yn
    METHYL_FEATURE_SELECTION(methyl_data_wide)

    // Stage 1: binary cancer / non-cancer classifier
    METHYL_CANCER_DETECTION(methyl_data_wide)

    // Stage 2: tissue of origin classifier (cancer samples only)
    METHYL_CANCER_TYPE(methyl_data_wide)
}