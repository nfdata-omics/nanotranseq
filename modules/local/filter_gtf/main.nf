process FILTER_GTF {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.novel.gtf"), emit: gtf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Filter for novel class codes:
    # j: Potentially novel isoform (at least one splice junction is shared with reference)
    # i: A transfrag falling entirely within a reference intron
    # o: Generic exonic overlap with a reference transcript
    # u: Unknown, intergenic transcript
    # x: Exonic overlap with reference on the opposite strand
    
    grep -E 'class_code "[jioux]"' $gtf > ${prefix}.novel.gtf || true
    
    # If file is empty (no novel transcripts), create empty file to avoid errors downstream
    if [ ! -s ${prefix}.novel.gtf ]; then
        echo "WARNING: No novel transcripts found." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | grep -oP '[0-9.]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.novel.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: "3.4"
    END_VERSIONS
    """
}
