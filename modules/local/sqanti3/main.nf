process SQANTI3 {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/anaconesalab/sqanti3:v6.0.1"

    input:
    tuple val(meta), path(isoforms_gtf)   // novel / merged transcript GTF to classify
    path reference_gtf                    // reference annotation GTF
    path fasta                            // reference genome FASTA

    output:
    tuple val(meta), path("*_classification.txt"), emit: classification
    tuple val(meta), path("*_corrected.gtf")     , emit: corrected_gtf
    tuple val(meta), path("*.pdf")               , emit: report, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sqanti3_qc.py \\
        --isoforms $isoforms_gtf \\
        --refGTF $reference_gtf \\
        --refFasta $fasta \\
        -t $task.cpus \\
        -o $prefix \\
        -d . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py --version 2>&1 | sed -n 's/.*SQANTI3 //p' | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_classification.txt
    touch ${prefix}_corrected.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: "6.0"
    END_VERSIONS
    """
}
