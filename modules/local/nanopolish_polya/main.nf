process NANOPOLISH_POLYA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::nanopolish=0.14.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.14.0--h773013f_3' :
        'biocontainers/nanopolish:0.14.0--h773013f_3' }"

    input:
    tuple val(meta), path(reads), path(index), path(fast5_dir), path(bam), path(bai)  // reads + nanopolish index sidecars, reads-to-genome BAM + index
    path fasta                                                       // reference genome FASTA

    output:
    tuple val(meta), path("*.polya.tsv"), emit: polya
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    nanopolish polya \\
        --threads $task.cpus \\
        --reads $reads \\
        --bam $bam \\
        --genome $fasta \\
        $args \\
        > ${prefix}.polya.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$(nanopolish --version 2>&1 | sed -n 's/^nanopolish version //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "readname\\tcontig\\tposition\\tleader_start\\tadapter_start\\tpolya_start\\ttranscript_start\\tread_rate\\tpolya_length\\tqc_tag" > ${prefix}.polya.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: "0.14.0"
    END_VERSIONS
    """
}
