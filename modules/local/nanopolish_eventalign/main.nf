process NANOPOLISH_EVENTALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::nanopolish=0.14.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.14.0--h773013f_3' :
        'biocontainers/nanopolish:0.14.0--h773013f_3' }"

    input:
    tuple val(meta), path(reads), path(index), path(fast5_dir), path(bam), path(bai)  // reads + nanopolish index sidecars + raw signal dir, reads-to-transcriptome BAM + index
    path transcript_fasta                                            // transcriptome FASTA (reads must be aligned to this)

    output:
    tuple val(meta), path("*.eventalign.txt"), emit: eventalign
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    nanopolish eventalign \\
        --threads $task.cpus \\
        --reads $reads \\
        --bam $bam \\
        --genome $transcript_fasta \\
        --scale-events \\
        --signal-index \\
        --summary ${prefix}.eventalign.summary.txt \\
        $args \\
        > ${prefix}.eventalign.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$(nanopolish --version 2>&1 | sed -n 's/^nanopolish version //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "contig\\tposition\\treference_kmer\\tread_index\\tstrand\\tevent_index\\tevent_level_mean" > ${prefix}.eventalign.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: "0.14.0"
    END_VERSIONS
    """
}
