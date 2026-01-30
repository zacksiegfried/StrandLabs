#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PRECISION_PROFILE {
    tag "${input_csv.simpleName}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path input_csv

    output:
    path "plots/*.png"
    path "${params.prefix}_stats.csv"

    script:
    """
    python3 ${workflow.projectDir}/src/precision_profile.py \
        --input ${input_csv} \
        --output . \
        --sdmc-col "${params.sdmc_col}" \
        --mutid-col "${params.mutid_col}" \
        --donor-col "${params.donor_col}" \
        --group-col "${params.group_col}" \
        --prefix "${params.prefix}" \
        --donors ${params.donors.split(',').join(' ')}
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

    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)

    // Run precision profile
    PRECISION_PROFILE(input_ch)
}
