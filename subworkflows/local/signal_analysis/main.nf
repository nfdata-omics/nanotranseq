include { NANOPOLISH_INDEX        } from '../../../modules/local/nanopolish_index/main'

//
// Direct-RNA signal preparation: build the nanopolish index (reads <-> raw FAST5 signal)
// and the shared eventalign artifact (signal aligned to the transcriptome).
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
    ch_indexed = NANOPOLISH_INDEX.out.indexed   // [meta, reads, index, fast5]

    emit:
    indexed    = ch_indexed      // channel: [meta, reads, index, fast5] for NANOPOLISH_POLYA
    versions
}
