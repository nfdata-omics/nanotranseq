process CPAT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cpat:3.0.5--py39hff726c5_3' :
        'biocontainers/cpat:3.0.5--py39hff726c5_4' }"

    input:
    tuple val(meta), path(fasta)            // From ch_filtered_exons_fa
    path hexamer                            // From ch_cpat_hexamer
    path logit_model                        // From ch_cpat_logit

    output:
    tuple val(meta), path("*.cpat.tsv"), emit: cpat_results
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run CPAT prediction
    cpat.py \\
        -g $fasta \\
        -x $hexamer \\
        -d $logit_model \\
        -o ${prefix}.cpat \\
        $args

    # Rename output
    mv ${prefix}.cpat ${prefix}.cpat.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpat: \$(cpat.py --version 2>&1 | grep -oP 'CPAT-\\K[0-9.]+' || echo "3.0.5")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cpat.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpat: "3.0.5"
    END_VERSIONS
    """
}
