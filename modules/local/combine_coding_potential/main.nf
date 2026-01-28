process COMBINE_CODING_POTENTIAL {
    tag "combine"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta1), path(cpat_results)
    tuple val(meta2), path(feelnc_results)
    tuple val(meta3), path(plek_results)
    tuple val(meta4), path(gtf)
    tuple val(meta5), path(fasta)

    output:
    tuple val(meta4), path("*.final_lncrna.gtf")   , emit: lncrna_gtf
    tuple val(meta5), path("*.final_lncrna.fa")    , emit: lncrna_fasta
    path "*.prediction_summary.tsv"                , emit: summary
    path "*.lncrna_report.txt"                     , emit: report
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta4.id}"
    def consensus_mode = task.ext.args ?: 'majority'  // majority, strict, lenient
    """
    combine_coding_potential.py \\
        --cpat $cpat_results \\
        --feelnc $feelnc_results \\
        --plek $plek_results \\
        --gtf $gtf \\
        --fasta $fasta \\
        --mode $consensus_mode \\
        --prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
