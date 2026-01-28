process GTF_TO_FASTA {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"

    input:
    tuple val(meta), path(gtf)
    path fasta

    output:
    tuple val(meta), path("*.transcripts.fa"), emit: fasta
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract transcript sequences from genome using GTF coordinates
    gffread \\
        -w ${prefix}.transcripts.fa \\
        -g $fasta \\
        $args \\
        $gtf
    if ! grep -q '^>' ${prefix}.transcripts.fa; then
        echo "ERROR: No '>' headers in ${prefix}.transcripts.fa; upstream filtering may have produced an empty GTF." 1>&2
        exit 2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcripts.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
