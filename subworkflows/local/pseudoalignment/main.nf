include { OARFISH            } from '../../../modules/local/oarfish/main'
include { SALMON_INDEX       } from '../../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT       } from '../../../modules/nf-core/salmon/quant/main'
include { TX2GENE_GTF        } from '../../../modules/local/tx2gene_gtf/main'
include { TXIMPORT           } from '../../../modules/local/tximport_transcript/main'

//
// ONT-native transcript quantification for DTU.
// Default is oarfish: an EM over long-read transcriptome alignments, reliable on noisy
// long reads. Salmon is available as an opt-in (`--pseudoaligner salmon`) but experimental
// on ONT — its k-mer selective alignment is unreliable on noisy long reads. tximport
// (type-aware) builds gene/transcript matrices from either quantifier.
//
workflow PSEUDOALIGNMENT {

    take:
    reference_gtf       // channel: reference GTF: path(gtf)
    transcriptome_bam   // channel: [meta, bam] reads-to-transcriptome (minimap2 -ax map-ont), used by oarfish
    reads               // channel: [meta, fastq] raw reads, used by salmon
    fasta               // channel: [meta, genome_fasta], used by salmon as decoy
    transcript_fasta    // channel: path(transcript_fasta), used by salmon

    main:
    versions = Channel.empty()

    //
    // Transcript quantification: oarfish (default) or salmon (opt-in, experimental).
    //
    if (params.pseudoaligner == 'salmon') {
        SALMON_INDEX(
            fasta.collect { it[1] },
            transcript_fasta
        )
        versions = versions.mix(SALMON_INDEX.out.versions)

        SALMON_QUANT(
            reads,
            SALMON_INDEX.out.index.first(),
            reference_gtf,
            transcript_fasta.first(),
            Channel.value(false),
            Channel.value(false)
        )
        versions = versions.mix(SALMON_QUANT.out.versions)

        ch_quant   = SALMON_QUANT.out.results   // [meta, results_dir]
        quant_type = Channel.value('salmon')
    }
    else {
        OARFISH(
            transcriptome_bam
        )
        versions = versions.mix(OARFISH.out.versions)

        ch_quant   = OARFISH.out.quant          // [meta, *.quant]
        quant_type = Channel.value('oarfish')
    }

    //
    // Build tx2gene mapping from the reference GTF.
    //
    TX2GENE_GTF(
        reference_gtf.map { gtf -> tuple([id: 'reference'], gtf) }
    )
    versions = versions.mix(TX2GENE_GTF.out.versions)

    //
    // Import all samples' quants into gene/transcript matrices (type-aware).
    //
    ch_quants = ch_quant
        .map { it[1] }
        .collect()
        .map { quants -> tuple([id: 'all_samples'], quants) }

    TXIMPORT(
        ch_quants,
        TX2GENE_GTF.out.tx2gene,
        quant_type
    )
    versions = versions.mix(TXIMPORT.out.versions)

    tpm_gene                  = TXIMPORT.out.tpm_gene
    counts_gene               = TXIMPORT.out.counts_gene
    counts_gene_length_scaled = TXIMPORT.out.counts_gene_length_scaled
    counts_gene_scaled        = TXIMPORT.out.counts_gene_scaled
    lengths_gene              = TXIMPORT.out.lengths_gene
    tpm_transcript            = TXIMPORT.out.tpm_transcript
    counts_transcript         = TXIMPORT.out.counts_transcript
    lengths_transcript        = TXIMPORT.out.lengths_transcript

    emit:
    versions

    tpm_gene
    counts_gene
    counts_gene_length_scaled
    counts_gene_scaled
    lengths_gene
    tpm_transcript
    counts_transcript
    lengths_transcript
    quant = ch_quant

}
