#!/usr/bin/env nextflow

process cleanup {
    """
    echo 'Cleaning up temporary Nextflow files...'
    rm -rf .nextflow .nextflow.log .nextflow.log.* work/
    echo 'Cleanup complete.'
    """
}

workflow {
    cleanup()
}