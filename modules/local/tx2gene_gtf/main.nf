process TX2GENE_GTF {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(gtf)   // reference GTF

    output:
    tuple val(meta), path("*.tx2gene.tsv"), emit: tx2gene
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'tx2gene.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'transcript_id\\tgene_id\\tgene_name\\n' > ${prefix}.tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3"
    END_VERSIONS
    """
}
