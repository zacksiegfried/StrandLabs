#!/usr/bin/env nextflow

process pscore_algorithm{

    input:
    path input_file

    output:
    path "${params.output_dir}/${input_file.getBaseName()}_pscore.csv"

    script:
    """
    python3 bin/scoring_algorithm.py --input $input_file --output ${params.output_dir}
    """

}

workflow {
    Channel
        .fromPath("${params.input}/*.csv")
        .set { variant_files }

    pscore_algorithm(variant_files)
}