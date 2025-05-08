#!/usr/bin/env nextflow

params.input = "data/*.csv"
params.outdir = "output"

process ScoreVariants{

    tag "$sample_id"

    input:
    path csv_file

    output:
    path "${csv_file.simpleName}_pscore.csv"

    script:
    sample_id = csv_file.simpleName
    """
    python3 bin/score_variants.py $csv_file > ${sample_id}_pscore.csv
    """

    stub:
    """
    
    """
}

workflow {
    Channel
        .fromPath(params.input)
        .set { variant_files }

    ScoreVariants(variant_files)
}