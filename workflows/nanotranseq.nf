/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RAW_READS_QC                    } from '../subworkflows/local/raw_read_qc/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { CDNA_QC                         } from '../subworkflows/local/cdna_qc/main'
include { ALIGNMENT                       } from '../subworkflows/local/alignment/main'
include { BEDTOOLS_BIGWIG                 } from '../subworkflows/local/bedtools_bigwig/main'
include { STRINGTIE_FEATURECOUNTS         } from '../subworkflows/local/stringtie_featurecounts/main'
include { NOVEL_TRANSCRIPTS               } from '../subworkflows/local/novel_transcripts/main'
include { IDENTIFY_NOVEL_PROTEIN_CODING   } from '../subworkflows/local/identify_novel_protein_coding/main'
include { EXTRACT_MRNA_SEQUENCES          } from '../modules/local/extract_mrna_sequences/main'
include { EXTRACT_CDS_SEQUENCES           } from '../modules/local/extract_cds_sequences/main'
include { EXTRACT_LNCRNA_SEQUENCES        } from '../modules/local/extract_lncrna_sequences/main'
include { DIFFERENTIAL_ANALYSIS           } from '../subworkflows/nfdata-omics/deseq2_analysis/main'
include { PSEUDOALIGNMENT                 } from '../subworkflows/local/pseudoalignment/main'
include { SIGNAL_ANALYSIS                 } from '../subworkflows/local/signal_analysis/main'
include { NANOPOLISH_POLYA                } from '../modules/local/nanopolish_polya/main'
include { TRANSCRIPT_USAGE                } from '../subworkflows/local/transcript_usage/main'
include { MINIMAP2_ALIGN as MINIMAP2_TRANSCRIPTOME } from '../modules/nf-core/minimap2/align/main'
include { UNTAR                           } from '../modules/nf-core/untar/main'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_nanotranseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOTRANSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta       // channel: fasta read in from --fasta
    ch_gtf          // channel: gtf file read in from --gtf
    ch_gene_id      // channel: attributed gene ID in the GTF file
    ch_gene_attributes     // channel: extra gene attributes in the GTF file
    ch_transcript_fasta     // channel: transcript fasta file read in from --transcript_fasta
    ch_direct_rna           // channel: direct_rna read in from --direct_rna
    ch_minimap2_index   // channel: index read in from --minimap2_index
    ch_formula       // channel: formula read in from --deseq2_formula
    ch_comparison    // channel: comparison read in from --deseq2_comparison
    ch_fdr_threshold // channel: fdr_threshold read in from --deseq2_fdr_threshold
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    //
    // MODULE: Run QC on raw reads
    //
    RAW_READS_QC (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(RAW_READS_QC.out.versions)

    //
    // Run Chopper if cDNA sequencing was performed
    if (!params.direct_rna) {

        CDNA_QC (
            ch_samplesheet,
        )

        ch_versions = ch_versions.mix(CDNA_QC.out.versions)

    }

    // If cDNA was performed, use CHOPPER's output as reads. If not, use raw data
    ch_reads = params.direct_rna ? ch_samplesheet : CDNA_QC.out.reads

    // Genome alignment (reads-to-genome BAM+BAI), set inside the featurecounts branch.
    ch_genome_bam = Channel.empty()

    //
    // Transcriptome alignment (map-ont), built once and reused by oarfish quant
    //
    ch_transcriptome_bam     = Channel.empty()   // [ meta, bam ]
    ch_transcriptome_bam_bai = Channel.empty()   // [ meta, bam, bai ]
    // Transcriptome BAM is needed by oarfish quant and by the eventalign step, which the
    // m6anet (RNA modifications) and xpore (RNA methylation) analyses consume.
    // poly(A) uses the genome alignment, so it does not require the transcriptome BAM.
    def needs_transcriptome_bam = params.quantification_tool == 'oarfish' ||
                                  params.quantification_tool == 'salmon' ||
                                  params.quantification_tool == 'both'
    if (needs_transcriptome_bam) {
        MINIMAP2_TRANSCRIPTOME(
            ch_reads,
            ch_transcript_fasta.map { [ [id: 'transcriptome'], it ] }.first(),
            true,   // bam_format
            '',  // bam_index_extension
            false,  // cigar_paf_format
            false   // cigar_bam
        )
        ch_versions = ch_versions.mix(MINIMAP2_TRANSCRIPTOME.out.versions)
        ch_transcriptome_bam     = MINIMAP2_TRANSCRIPTOME.out.bam
        //ch_transcriptome_bam_bai = MINIMAP2_TRANSCRIPTOME.out.bam.join(MINIMAP2_TRANSCRIPTOME.out.index)
    }

    //
    // Run alignment if either `featurecounts` or `both` is selected as quantification tool
    //
    if (params.quantification_tool == 'featurecounts' || params.quantification_tool == 'both') {

        // Run alignment with Minimap2
        ALIGNMENT(
            ch_reads,
            ch_fasta,
            ch_minimap2_index
        )
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

        // Genome BAM+BAI for signal-level poly(A) estimation
        ch_genome_bam = ALIGNMENT.out.minimap2_bam.join(ALIGNMENT.out.minimap2_bai)

        // Generate BigWig files for visualisation in genome browsers
        BEDTOOLS_BIGWIG(
            ch_fasta,
            ALIGNMENT.out.minimap2_bam
        )

        // Assemble and quantify
        STRINGTIE_FEATURECOUNTS(
            ch_fasta,
            ALIGNMENT.out.minimap2_bam,
            ch_gtf,
        )
        ch_versions = ch_versions.mix(STRINGTIE_FEATURECOUNTS.out.versions)

        // Classify merged StringTie transcripts against the reference annotation
        // and extract candidate novel transcript sequences.
        NOVEL_TRANSCRIPTS(
            STRINGTIE_FEATURECOUNTS.out.merged_gtf,
            ch_fasta,
            ch_gtf
        )
        ch_versions = ch_versions.mix(NOVEL_TRANSCRIPTS.out.versions)

        if (params.run_coding_potential) {
            if (params.skip_cpat || params.skip_feelnc || params.skip_plek) {
                exit 1, "Coding-potential skips are not supported yet because COMBINE_PREDICTIONS currently expects CPAT, FEELnc and PLEK outputs."
            }
            has_cpat_models = params.cpat_hexamer && params.cpat_logit_model
            has_partial_cpat_models = params.cpat_hexamer || params.cpat_logit_model
            if (has_partial_cpat_models && !has_cpat_models) {
                exit 1, "Coding-potential analysis requires both --cpat_hexamer and --cpat_logit_model when using pre-built CPAT models."
            }

            EXTRACT_MRNA_SEQUENCES(
                ch_fasta.map { meta, fasta_file -> fasta_file },
                ch_gtf
            )
            ch_versions = ch_versions.mix(EXTRACT_MRNA_SEQUENCES.out.versions)

            EXTRACT_CDS_SEQUENCES(
                ch_fasta.map { meta, fasta_file -> fasta_file },
                ch_gtf
            )
            ch_versions = ch_versions.mix(EXTRACT_CDS_SEQUENCES.out.versions)

            if (!params.cpat_training_noncoding_fasta) {
                EXTRACT_LNCRNA_SEQUENCES(
                    ch_fasta.map { meta, fasta_file -> fasta_file },
                    ch_gtf,
                    Channel.value(params.lncrna_biotypes)
                )
                ch_versions = ch_versions.mix(EXTRACT_LNCRNA_SEQUENCES.out.versions)
            }

            ch_cpat_training_coding_fasta = params.cpat_training_coding_fasta ?
                Channel.fromPath(params.cpat_training_coding_fasta) :
                EXTRACT_CDS_SEQUENCES.out.fasta
            ch_cpat_training_noncoding_fasta = params.cpat_training_noncoding_fasta ?
                Channel.fromPath(params.cpat_training_noncoding_fasta) :
                EXTRACT_LNCRNA_SEQUENCES.out.fasta
            ch_feelnc_mrna_fasta = params.feelnc_mrna_fasta ?
                Channel.fromPath(params.feelnc_mrna_fasta) :
                EXTRACT_MRNA_SEQUENCES.out.fasta

            IDENTIFY_NOVEL_PROTEIN_CODING(
                NOVEL_TRANSCRIPTS.out.tmap,
                NOVEL_TRANSCRIPTS.out.novel_gtf,
                NOVEL_TRANSCRIPTS.out.novel_fasta,
                ch_cpat_training_coding_fasta,
                ch_cpat_training_noncoding_fasta,
                ch_feelnc_mrna_fasta
            )
            ch_versions = ch_versions.mix(IDENTIFY_NOVEL_PROTEIN_CODING.out.versions)
        }

        // Create metadata for DESeq2
        ch_metadata = channel
            .fromPath(params.input)
            .map { samplesheet ->
                def metadata = file(workDir + '/metadata.tsv')
                metadata.text = samplesheet.text.readLines()
                    .collect { row -> row.replace(',', '\t') }
                    .join('\n') + '\n'
                tuple([id: 'deseq2'], metadata)
            }

        // Create counts file for DESeq2
        STRINGTIE_FEATURECOUNTS.out.featurecounts_genes_out
            .map { meta, featurecounts_file ->
                def counts_file = file(workDir + '/counts_for_deseq2.tsv')
                def lines = featurecounts_file.text.readLines()

                // Remove comment line if present
                if (lines[0].startsWith('#')) lines = lines[1..-1]

                // Remove unnecessary columns and strip '.bam' suffix from header
                counts_file.text = lines.collect { line ->
                    def cols = line.split('\t') as List
                    def row = ([cols[0]] + cols[6..-1]).join('\t')
                    line == lines[0] ? row.replaceAll('\\.bam', '') : row
                }.join('\n') + '\n'

                tuple([id: meta], counts_file)
            }
            .set { ch_counts_for_deseq2 }

        // Differential analysis
        DIFFERENTIAL_ANALYSIS(
            ch_counts_for_deseq2,
            ch_metadata,
            ch_formula,
            ch_comparison,
            ch_fdr_threshold
        )
    }

    //
    // Run pseudoalignment if either `salmon` or `both` is selected as quantification tool
    //
    if (params.quantification_tool == 'salmon' || params.quantification_tool == 'both') {

        PSEUDOALIGNMENT(
            ch_gtf,
            ch_transcriptome_bam,
            ch_reads,
            ch_fasta,
            ch_transcript_fasta,
        )
        ch_versions = ch_versions.mix(PSEUDOALIGNMENT.out.versions)

        //
        // SUBWORKFLOW: Transcript Usage
        //
        TRANSCRIPT_USAGE(
            PSEUDOALIGNMENT.out.counts_transcript,
            PSEUDOALIGNMENT.out.tpm_transcript,
            ch_fasta.map{ it[1] },
            ch_gtf,
            ch_reads
        )
        ch_versions = ch_versions.mix(TRANSCRIPT_USAGE.out.versions)

    }

    //
    // SUBWORKFLOW: Direct-RNA signal preparation (nanopolish index).
    // Runs automatically when poly(A) is enabled
    //
    if (params.run_polya) {
        // Read the fast5 signal directory per sample from the samplesheet
        def sheet_dir = file(params.input).parent
        ch_fast5_in = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                if (!row.fast5) {
                    error("--run_polya requires a 'fast5' column in the samplesheet (sample: ${row.sample})")
                }
                def f5 = row.fast5 ==~ /^(\/|[a-zA-Z][a-zA-Z0-9+.-]*:\/\/).*/ ?
                    file(row.fast5, checkIfExists: true) :
                    file("${sheet_dir}/${row.fast5}", checkIfExists: true)
                tuple(row.sample, f5)
            }
            // The fast5 column takes either a directory or a tar archive of one. The
            // archive form exists because Nextflow cannot stage a remote directory, so
            // signal data hosted outside the repo has to travel as a single file.
            .branch { _id, f5 ->
                archive: f5.name ==~ /.*\.(tar\.gz|tgz|tar)$/
                dir    : true
            }

        // Unpack only the archives; directories pass through untouched.
        UNTAR( ch_fast5_in.archive.map { id, f5 -> tuple([ id: id ], f5) } )

        ch_fast5 = ch_fast5_in.dir
            .mix( UNTAR.out.untar.map { meta, untarred -> tuple(meta.id, untarred) } )

        // Attach fast5 to reads: [ meta, reads, fast5 ]
        ch_signal = ch_reads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(ch_fast5)
            .map { _id, meta, reads, fast5 -> tuple(meta, reads, fast5) }

        SIGNAL_ANALYSIS(
            ch_signal,
            ch_transcriptome_bam,
            ch_transcript_fasta,
            params.run_rna_modifications || params.run_rna_methylation
        )
        ch_versions = ch_versions.mix(SIGNAL_ANALYSIS.out.versions)

        //
        // Poly(A) tail length (per read, genome alignment). Single module, called directly.
        //
        if (params.run_polya) {
            ch_polya_in = SIGNAL_ANALYSIS.out.indexed.join(ch_genome_bam)   // [ meta, reads, index, fast5, bam, bai ]
            NANOPOLISH_POLYA(
                ch_polya_in,
                ch_fasta.map { it[1] }
            )
            ch_versions = ch_versions.mix(NANOPOLISH_POLYA.out.versions)
        }

    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'nanotranseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    multiqc_report = RAW_READS_QC.out.multiqc_report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
