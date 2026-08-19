process OARFISH {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/oarfish:0.10.0--hd727d2a_0' :
        'biocontainers/oarfish:0.10.0--hd727d2a_0' }"

    input:
    tuple val(meta), path(bam)   // reads aligned to the transcriptome (minimap2 -ax map-ont)

    output:
    tuple val(meta), path("*.quant"), emit: quant
    tuple val(meta), path("*.meta_info.json"), emit: meta_info
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    oarfish \\
        --threads $task.cpus \\
        --alignments $bam \\
        --output $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version 2>&1 | sed -n 's/^oarfish //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.quant
    echo '{}' > ${prefix}.meta_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: "0.6.2"
    END_VERSIONS
    """
}
