#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_ASV } from './modules/local/validate_asv/main'

params.counts = null
params.taxonomy = null
params.outdir = 'results'

workflow {

    if (!params.counts || !params.taxonomy) {
        error "Both --counts and --taxonomy are required."
    }

    log.info "ASV Explorer v0.2.0"
    log.info "Counts: ${params.counts}"
    log.info "Taxonomy: ${params.taxonomy}"
    log.info "Output directory: ${params.outdir}"

    counts_ch = Channel.fromPath(
        params.counts,
        checkIfExists: true
    )

    taxonomy_ch = Channel.fromPath(
        params.taxonomy,
        checkIfExists: true
    )

    VALIDATE_ASV(counts_ch, taxonomy_ch)

}