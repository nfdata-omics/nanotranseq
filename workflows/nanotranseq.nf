/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RAW_READS_QC           } from '../subworkflows/local/raw_read_qc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { DIRECT_RNA_QC          } from '../subworkflows/local/direct_rna_qc/main'
include { ALIGNMENT              } from '../subworkflows/local/alignment/main'
include { STRINGTIE_STRINGTIE    } from '../modules/nf-core/stringtie/stringtie/main'
include { CODING_POTENTIAL       } from '../subworkflows/local/coding_potential/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanotranseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOTRANSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta       // channel: fasta read in from --fasta
    ch_gtf         // channel: gtf read in from --gtf
    ch_direct_rna // channel: direct_rna read in from --direct_rna
    ch_minimap2_index   // channel: index read in from --minimap2_index
    ch_cpat_hexamer // channel: cpat_hexamer
    ch_cpat_logit   // channel: cpat_logit

    main:

    ch_versions = channel.empty()

    //
    // MODULE: Run QC on raw reads
    //
    RAW_READS_QC (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(RAW_READS_QC.out.versions)

    //
    // Run Chopper if direct RNA sequencing was performed
    if (ch_direct_rna) {

        DIRECT_RNA_QC (
            ch_samplesheet,
        )

        ch_versions = ch_versions.mix(DIRECT_RNA_QC.out.versions)

    }

    // If direct RNA was performed, use CHOPPER's output as reads. If not, use raw data
    ch_reads = ch_direct_rna ? DIRECT_RNA_QC.out.reads : ch_samplesheet

    //
    // Run alignment
    //
    ALIGNMENT(ch_reads,
              ch_fasta,
              ch_minimap2_index
              )

    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    //
    // Assembly: StringTie
    //
    STRINGTIE_STRINGTIE(
        ALIGNMENT.out.minimap2_bam,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

    //
    // Coding Potential
    //
    CODING_POTENTIAL(
        STRINGTIE_STRINGTIE.out.transcript_gtf,
        ch_fasta,
        ch_gtf,
        ch_cpat_hexamer,
        ch_cpat_logit
    )
    ch_versions = ch_versions.mix(CODING_POTENTIAL.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
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
            storeDir: "${params.outdir}/pipeline_info",
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
