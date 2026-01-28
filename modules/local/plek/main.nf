process PLEK {
    tag "plek"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plek:1.2--py39he88f293_9' :
        'biocontainers/plek:1.2--py311h8ddd9a4_10' }"

    input:
    tuple val(meta), path(fasta)        // From ch_filtered_exons_fa

    output:
    tuple val(meta), path("*.plek.tsv"), emit: plek_results
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    # Run PLEK
    PLEK.py \\
        -fasta $fasta \\
        -out ${prefix}.plek.txt \\
        -thread $task.cpus \\
        $args

    # Convert to TSV format
    convert_plek_output.py \\
        --input ${prefix}.plek.txt \\
        --output ${prefix}.plek.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: \$(PLEK.py -v 2>&1 | grep -oP 'PLEK \\K[0-9.]+' || echo "1.2")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plek.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: "1.2"
    END_VERSIONS
    """
}
