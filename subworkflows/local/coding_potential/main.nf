//
// Check coding potential of transcript isoforms
//

include { GTF_TO_FASTA } from '../../modules/local/gtf_to_fasta/main'
include { GFFCOMPARE } from '../../modules/local/gffcompare/main'
include { FILTER_GTF } from '../../modules/local/filter_gtf/main'
include { CPAT } from '../../modules/local/cpat/main'
include { PLEK } from '../../modules/local/plek/main'
include { FEELNC_CODPOT } from '../../modules/local/feelnc_codpot/main'
include { COMBINE_CODING_POTENTIAL } from '../../modules/local/combine_coding_potential/main'

workflow CODING_POTENTIAL {
    take:
    ch_gtf          // channel: [ val(meta), [ gtf ] ]
    ch_fasta        // channel: [ fasta ] (genome)
    ch_ref_gtf      // channel: [ gtf ] (reference annotation)
    ch_cpat_hexamer // channel: [ hexamer ]
    ch_cpat_logit   // channel: [ logit ]

    main:

    ch_versions = Channel.empty()

    //
    // Classify transcripts against reference
    //
    GFFCOMPARE (
        ch_gtf,
        ch_fasta,
        ch_ref_gtf
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

    //
    // Filter for novel isoforms only (j, i, o, u, x)
    //
    FILTER_GTF (
        GFFCOMPARE.out.annotated_gtf
    )
    ch_versions = ch_versions.mix(FILTER_GTF.out.versions)

    //
    // Extract transcript sequences (only for novel ones)
    //
    GTF_TO_FASTA (
        FILTER_GTF.out.gtf,
        ch_fasta
    )
    ch_versions = ch_versions.mix(GTF_TO_FASTA.out.versions)

    //
    // Run CPAT
    //
    CPAT (
        GTF_TO_FASTA.out.fasta,
        ch_cpat_hexamer,
        ch_cpat_logit
    )
    ch_versions = ch_versions.mix(CPAT.out.versions)

    //
    // Run PLEK
    //
    PLEK (
        GTF_TO_FASTA.out.fasta
    )
    ch_versions = ch_versions.mix(PLEK.out.versions)

    //
    // Run FEELnc
    //
    FEELNC_CODPOT (
        GTF_TO_FASTA.out.fasta,
        ch_fasta
    )
    ch_versions = ch_versions.mix(FEELNC_CODPOT.out.versions)

    //
    // Combine predictions
    //

    // Join all results by meta
    // Note: We use FILTER_GTF.out.gtf here instead of original ch_gtf

    ch_combined_inputs = CPAT.out.cpat_results
        .join(FEELNC_CODPOT.out.feelnc_results)
        .join(PLEK.out.plek_results)
        .join(FILTER_GTF.out.gtf)
        .join(GTF_TO_FASTA.out.fasta)

    // ch_combined_inputs structure: [ meta, cpat, feelnc, plek, gtf, transcript_fasta ]

    COMBINE_CODING_POTENTIAL (
        ch_combined_inputs.map { [it[0], it[1]] }, // cpat
        ch_combined_inputs.map { [it[0], it[2]] }, // feelnc
        ch_combined_inputs.map { [it[0], it[3]] }, // plek
        ch_combined_inputs.map { [it[0], it[4]] }, // gtf
        ch_combined_inputs.map { [it[0], it[5]] }  // fasta
    )
    ch_versions = ch_versions.mix(COMBINE_CODING_POTENTIAL.out.versions)

    emit:
    lncrna_gtf   = COMBINE_CODING_POTENTIAL.out.lncrna_gtf
    lncrna_fasta = COMBINE_CODING_POTENTIAL.out.lncrna_fasta
    summary      = COMBINE_CODING_POTENTIAL.out.summary
    report       = COMBINE_CODING_POTENTIAL.out.report
    versions     = ch_versions
}
