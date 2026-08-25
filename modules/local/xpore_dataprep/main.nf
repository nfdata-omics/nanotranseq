process XPORE_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

    // Tag verified against quay.io/biocontainers (2026-08-24).
    conda "bioconda::xpore=2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh5e36f6f_0' :
        'biocontainers/xpore:2.1--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(eventalign)   // NANOPOLISH_EVENTALIGN.out.eventalign (transcriptomic coords)

    output:
    tuple val(meta), path("dataprep_${meta.id}"), emit: dataprep   // per-sample dir consumed by XPORE_DIFFMOD
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    xpore dataprep \\
        --eventalign $eventalign \\
        --out_dir dataprep_${prefix} \\
        --n_processes $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed -n 's/^xpore //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p dataprep_${prefix}
    touch dataprep_${prefix}/data.json dataprep_${prefix}/data.index dataprep_${prefix}/data.readcount dataprep_${prefix}/eventalign.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: "2.1"
    END_VERSIONS
    """
}
