process PLEK_RUN {
    tag "plek_run"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plek:1.2--py39he88f293_9' :
        'biocontainers/plek:1.2--py39he88f293_9' }"

    input:
    tuple val(meta), path(fasta)        // FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta

    output:
    tuple val(meta), path("*.plek.txt") ,  emit: plek_raw
    tuple val(meta), path(fasta)        ,  emit: fasta
    path "versions.yml"                 ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    PLEK.py \\
        -fasta $fasta \\
        -out ${prefix}.plek.txt \\
        -thread $task.cpus \\
        $args || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: "1.2"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plek.txt
    rm -f ${fasta}
    touch ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: "1.2"
    END_VERSIONS
    """
}
