//
// QUALITY CONTROL FOR DIRECT RNA WORKFLOW
//

include { CHOPPER     } from '../../../modules/nf-core/chopper/main'
include { FASTQC      } from '../../../modules/nf-core/fastqc/main'
include { TOULLIGQC   } from '../../../modules/nf-core/toulligqc/main'
include { NANOPLOT    } from '../../../modules/nf-core/nanoplot/main'
include { MULTIQC     } from '../../../modules/nf-core/multiqc/main'

workflow CDNA_QC {
    take:
    raw_reads   // raw reads input channel

    main:
    versions = Channel.empty()

    // Run CHOPPER
    CHOPPER (
        raw_reads,
        channel.value(file("no_fasta", checkIfExists: false))
    )

    // Create channel to store chopper's output
    reads = CHOPPER.out.fastq

    // Run FASTQC
    FASTQC(reads)
    fastqc_zip = FASTQC.out.zip
    fastqc_html = FASTQC.out.html

    MULTIQC(
          FASTQC.out.zip
              .collect { it[1] }
              .map { files -> tuple([:], files, [], [], [], []) }
    )

    // Run TOULLIGQC
    TOULLIGQC(reads)
    toulligqc_report_data = TOULLIGQC.out.report_data
    toulligqc_report_html = TOULLIGQC.out.report_html
    toulligqc_plots_html  = TOULLIGQC.out.plots_html
    toulligqc_plotly_js   = TOULLIGQC.out.plotly_js

    // Run Nanoplot
    NANOPLOT(reads)
    nanoplot_html = NANOPLOT.out.html
    nanoplot_png  = NANOPLOT.out.png
    nanoplot_txt  = NANOPLOT.out.txt


    // Collect versions for all tools used in this workflow
    // MULTIQC emits an eval-tuple (not versions.yml) and can't use the versions topic
    // (deadlock: its input depends on the topic). Convert it to a YAML string so
    // softwareVersionsToYAML can parse it like the file-based versions.
    ch_multiqc_version = MULTIQC.out.versions.map { process, tool, ver -> "\"${process}\":\n  ${tool}: ${ver}\n" }
    versions = versions.mix(ch_multiqc_version, TOULLIGQC.out.versions, NANOPLOT.out.versions, CHOPPER.out.versions)

    emit:

    reads

    fastqc_zip
    fastqc_html

    multiqc_report = MULTIQC.out.report

    toulligqc_report_data
    toulligqc_report_html
    toulligqc_plots_html
    toulligqc_plotly_js

    nanoplot_html
    nanoplot_png
    nanoplot_txt

    versions
}
