process NANOPOLISH_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::nanopolish=0.14.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.14.0--h773013f_3' :
        'biocontainers/nanopolish:0.14.0--h773013f_3' }"

    input:
    tuple val(meta), path(reads), path(fast5_dir)   // basecalled reads + raw signal dir

    output:
    tuple val(meta), path(reads), path("${reads}.index*"), path(fast5_dir)  , emit: indexed   // reads + sidecar index files
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    nanopolish index \\
        -d $fast5_dir \\
        $args \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$(nanopolish --version 2>&1 | sed -n 's/^nanopolish version //p')
    END_VERSIONS
    """

    stub:
    """
    touch ${reads}.index ${reads}.index.fai ${reads}.index.gzi ${reads}.index.readdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: "0.14.0"
    END_VERSIONS
    """
}
