#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfdata-omics/nanotranseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfdata-omics/nanotranseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NANOTRANSEQ  } from './workflows/nanotranseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_nanotranseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_nanotranseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_nanotranseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')
params.minimap2_index = getGenomeAttribute('minimap2')
params.gtf = getGenomeAttribute('gtf')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFDATAOMICS_NANOTRANSEQ {

    take:
    samplesheet // channel: samplesheet read in from --input
    fasta       // channel: fasta read in from --fasta
    gtf         // channel: gtf read in from --gtf
    direct_rna  // channel: direct_rna read in from --direct_rna
    minimap2_index  // channel: minimap2_index read in from --minimap2_index
    cpat_hexamer // channel: cpat_hexamer
    cpat_logit   // channel: cpat_logit

    main:

    //
    // WORKFLOW: Run pipeline
    //
    NANOTRANSEQ (
        samplesheet,
        fasta,
        gtf,
        direct_rna,
        minimap2_index,
        cpat_hexamer,
        cpat_logit
    )
    emit:
    multiqc_report = NANOTRANSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.fasta,
        params.direct_rna,
        params.minimap2_index,
        params.help,
        params.help_full,
        params.show_hidden
    )

    ch_gtf = params.gtf ? Channel.fromPath(params.gtf) : Channel.empty()
    ch_cpat_hexamer = params.cpat_hexamer ? Channel.fromPath(params.cpat_hexamer) : Channel.empty()
    ch_cpat_logit = params.cpat_logit ? Channel.fromPath(params.cpat_logit) : Channel.empty()

    //
    // WORKFLOW: Run main workflow
    //
    NFDATAOMICS_NANOTRANSEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.fasta,
        ch_gtf,
        PIPELINE_INITIALISATION.out.direct_rna,
        PIPELINE_INITIALISATION.out.minimap2_index,
        ch_cpat_hexamer,
        ch_cpat_logit
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFDATAOMICS_NANOTRANSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
