include { NANOPOLISH_INDEX        } from '../../../modules/local/nanopolish_index/main'
include { NANOPOLISH_EVENTALIGN   } from '../../../modules/local/nanopolish_eventalign/main'
include { SAMTOOLS_SORT           } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'

//
// Direct-RNA signal preparation: build the nanopolish index (reads <-> raw FAST5 signal)
// and the shared eventalign artifact (signal aligned to the transcriptome).
// Emitted `indexed` feeds NANOPOLISH_POLYA; emitted `eventalign` feeds M6ANET (RNA
// modifications) and RNA_METHYLATION (xpore).
//
workflow SIGNAL_ANALYSIS {
    take:
    signal              // channel: [ meta, reads, fast5_dir ]
    transcriptome_bam   // channel: [ meta, bam ] name-sorted reads-to-transcriptome (natural minimap2 order)
    transcript_fasta    // channel: path(transcript_fasta)
    run_eventalign      // value:   run the (expensive) eventalign step (modifications/methylation)

    main:
    versions = Channel.empty()

    //
    // Build the nanopolish index once (reads <-> raw signal), reused downstream.
    //
    NANOPOLISH_INDEX(signal)
    versions = versions.mix(NANOPOLISH_INDEX.out.versions)
    ch_indexed = NANOPOLISH_INDEX.out.indexed   // [ meta, reads, index, fast5 ]

    //
    // Signal-to-transcriptome alignment (eventalign): shared by the modification analyses.
    // Skipped when only poly(A) is requested (poly(A) uses the genome alignment instead).
    //
    ch_eventalign = Channel.empty()
    if (run_eventalign) {
        // nanopolish needs a coordinate-sorted + indexed BAM, but the transcriptome BAM
        // arrives name-sorted (for oarfish). Re-sort by coordinate, then index.
        SAMTOOLS_SORT(transcriptome_bam, [[], [], []], '')

        SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

        // [ meta, reads, index, fast5, bam, bai ] -- the fast5 dir travels with the
        // reads: `nanopolish index` writes the readdb with paths relative to it
        // (`<dir>/<read>.fast5`), so eventalign cannot open the signal without it staged.
        ch_eventalign_in = ch_indexed
            .join(SAMTOOLS_SORT.out.bam)
            .join(SAMTOOLS_INDEX.out.index)

        NANOPOLISH_EVENTALIGN(
            ch_eventalign_in,
            transcript_fasta
        )
        versions = versions.mix(NANOPOLISH_EVENTALIGN.out.versions)
        ch_eventalign = NANOPOLISH_EVENTALIGN.out.eventalign
    }

    emit:
    indexed    = ch_indexed      // channel: [ meta, reads, index, fast5 ] for NANOPOLISH_POLYA
    eventalign = ch_eventalign   // channel: [ meta, path(eventalign.txt) ]
    versions
}
