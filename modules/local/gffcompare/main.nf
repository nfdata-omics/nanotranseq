process GFFCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h4a3c31e_3' :
        'biocontainers/gffcompare:0.12.6--h4a3c31e_3' }"

    input:
    tuple val(meta), path(gtf)
    path fasta
    path reference_gtf

    output:
    tuple val(meta), path("*.annotated.gtf"), emit: annotated_gtf
    tuple val(meta), path("*.stats")        , emit: stats
    tuple val(meta), path("*.tracking")     , emit: tracking
    tuple val(meta), path("*.loci")         , emit: loci
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gffcompare \\
        -r $reference_gtf \\
        -s $fasta \\
        -o $prefix \\
        $gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(gffcompare --version 2>&1 | grep -o 'gffcompare v[0-9.]*' | sed 's/gffcompare v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotated.gtf
    touch ${prefix}.stats
    touch ${prefix}.tracking
    touch ${prefix}.loci

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: "0.12.6"
    END_VERSIONS
    """
}
