#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FILTER_VARIANTS {
    tag "${input_csv.simpleName}"

    publishDir "${workflow.projectDir}/data", mode: 'copy'

    input:
    path input_csv

    output:
    path "filtered_data.csv", emit: filtered_csv
    path "filtered_variants.txt"

    script:
    """
    python3 ${workflow.projectDir}/src/filter_variants.py \
        --input ${input_csv} \
        --output . \
        --maf-col "${params.maf_col}" \
        --mutid-col "${params.mutid_col}" \
        --donor-col "${params.donor_col}" \
        --dilution-col "${params.dilution_col}" \
        --donors ${params.donors.split(',').join(' ')}
    """
}


process PRECISION_PROFILE {
    tag "${input_csv.simpleName}"

    publishDir "${params.outdir}/precision", mode: 'copy'

    input:
    path input_csv

    output:
    path "plots/*.png"
    path "precision_stats.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/precision_profile.py \
        --input ${input_csv} \
        --output . \
        --maf-col "${params.maf_col}" \
        --mutid-col "${params.mutid_col}" \
        --donor-col "${params.donor_col}" \
        --dilution-col "${params.dilution_col}" \
        --donors ${params.donors.split(',').join(' ')}
    """
}


process HIT_RATE_LOD {
    tag "${input_csv.simpleName}"

    publishDir "${params.outdir}/lod", mode: 'copy'

    input:
    path input_csv

    output:
    path "plots/*.png"
    path "hit_rate_stats.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/hit_rate_lod.py \
        --input ${input_csv} \
        --output . \
        --alt-n-col "${params.alt_n_col}" \
        --mutid-col "${params.mutid_col}" \
        --donor-col "${params.donor_col}" \
        --dilution-col "${params.dilution_col}" \
        --donors ${params.donors.split(',').join(' ')} \
        --hit-threshold ${params.hit_threshold}
    """
}


workflow {
    // Validate required inputs
    if (!params.input) {
        error "ERROR: --input is required. Please provide an input CSV file."
    }
    if (!params.donors) {
        error "ERROR: --donors is required. Please provide a comma-separated list of donor IDs."
    }
    if (!params.alt_n_col) {
        error "ERROR: --alt_n_col is required. Please provide the allele count column name."
    }
    if (!params.maf_col) {
        error "ERROR: --maf_col is required. Please provide the MAF column name."
    }
    if (!params.mutid_col) {
        error "ERROR: --mutid_col is required. Please provide the mutation ID column name."
    }
    if (!params.donor_col) {
        error "ERROR: --donor_col is required. Please provide the donor ID column name."
    }
    if (!params.dilution_col) {
        error "ERROR: --dilution_col is required. Please provide the dilution level column name."
    }

    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)

    // Filter variants
    filtered_ch = FILTER_VARIANTS(input_ch).filtered_csv

    // Run analyses on filtered data
    PRECISION_PROFILE(filtered_ch)
    HIT_RATE_LOD(filtered_ch)
}
