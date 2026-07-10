include { SAMTOOLS_FAIDX                                } from '../../../modules/nf-core/samtools/faidx/main'
include { GFFCOMPARE as GFFCOMPARE_NOVEL_TRANSCRIPTS    } from '../../../modules/nf-core/gffcompare/main'
include { GFFREAD    as GFFREAD_NOVEL_TRANSCRIPTS       } from '../../../modules/nf-core/gffread/main'
include { FILTER_GFFCOMPARE_CLASSES                     } from '../../../modules/local/filter_gffcompare_classes'
include { SQANTI3                                       } from '../../../modules/local/sqanti3/main'

workflow NOVEL_TRANSCRIPTS {
    take:
        merged_gtf      // StringTie merged GTF
        fasta           // reference FASTA: tuple(meta), path(fasta)
        reference_gtf   // reference GTF: path(gtf)

    main:
        versions = Channel.empty()

        ch_meta = Channel.value([id: 'novel_transcripts'])

        //
        // Index the reference FASTA for gffcompare.
        //
        SAMTOOLS_FAIDX(
            fasta.map { meta, fasta_file -> tuple(meta, fasta_file, []) },
            false
        )

        ch_gffcompare_fasta = fasta
            .join(SAMTOOLS_FAIDX.out.fai)
            .map { meta, fasta_file, fai ->
                tuple([id: meta.toString()], fasta_file, fai)
            }

        //
        // Compare merged StringTie transcripts against the reference annotation.
        //
        ch_merged_gtf_for_gffcompare = merged_gtf
            .combine(ch_meta)
            .map { gtf, meta -> tuple(meta, [gtf]) }

        ch_reference_gtf_for_gffcompare = reference_gtf
            .combine(ch_meta)
            .map { gtf, meta -> tuple(meta, gtf) }

        GFFCOMPARE_NOVEL_TRANSCRIPTS(
            ch_merged_gtf_for_gffcompare,
            ch_gffcompare_fasta,
            ch_reference_gtf_for_gffcompare
        )
        versions = versions.mix(GFFCOMPARE_NOVEL_TRANSCRIPTS.out.versions)

        //
        // Keep candidate novel transcripts by gffcompare class code.
        //
        ch_merged_gtf_for_filter = merged_gtf
            .combine(ch_meta)
            .map { gtf, meta -> tuple(meta, gtf) }

        FILTER_GFFCOMPARE_CLASSES(
            ch_merged_gtf_for_filter,
            GFFCOMPARE_NOVEL_TRANSCRIPTS.out.tmap,
            Channel.value(params.novel_class_codes)
        )
        versions = versions.mix(FILTER_GFFCOMPARE_CLASSES.out.versions)

        //
        // Classify + filter isoform artifacts (intra-priming, RT-switching)
        // with SQANTI3. Its corrected GTF replaces the raw filtered GTF downstream.
        //

        ch_fasta_for_sqanti = fasta.map { _meta, fasta_file -> fasta_file }

        ch_selected_gtf = FILTER_GFFCOMPARE_CLASSES.out.filtered_gtf
        sqanti_classification = channel.empty()

        if (params.run_sqanti) {
            SQANTI3(
                FILTER_GFFCOMPARE_CLASSES.out.filtered_gtf,
                reference_gtf,
                ch_fasta_for_sqanti
            )
            versions = versions.mix(SQANTI3.out.versions)
            ch_selected_gtf = SQANTI3.out.corrected_gtf
            sqanti_classification = SQANTI3.out.classification
        }

        //
        // Extract transcript sequences for the selected novel candidates.
        //
        ch_fasta_for_gffread = fasta.map { meta, fasta_file -> fasta_file }

        GFFREAD_NOVEL_TRANSCRIPTS(
            ch_selected_gtf,
            ch_fasta_for_gffread
        )
        versions = versions.mix(GFFREAD_NOVEL_TRANSCRIPTS.out.versions)

    emit:
        tmap                    = GFFCOMPARE_NOVEL_TRANSCRIPTS.out.tmap
        annotated_gtf           = GFFCOMPARE_NOVEL_TRANSCRIPTS.out.annotated_gtf
        novel_gtf               = ch_selected_gtf
        novel_tmap              = FILTER_GFFCOMPARE_CLASSES.out.filtered_tmap
        class_summary           = FILTER_GFFCOMPARE_CLASSES.out.class_summary
        novel_fasta             = GFFREAD_NOVEL_TRANSCRIPTS.out.gffread_fasta
        sqanti_classification   = sqanti_classification
        versions                = versions
}
