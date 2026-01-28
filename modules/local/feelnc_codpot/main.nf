process FEELNC_CODPOT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/feelnc:0.2--pl526_0' :
        'biocontainers/feelnc:0.2--pl526_0' }"

    input:
    tuple val(meta), path(candidate_fasta)        // From ch_filtered_exons_fa
    path genome_fasta                             // From ch_fasta

    output:
    tuple val(meta), path("*.feelnc.tsv")          , emit: feelnc_results
    path "*.feelnc_codpot.txt"                     , emit: codpot_full
    path "*.lncRNA.fa"                             , emit: lncrna_fasta, optional: true
    path "*.mRNA.fa"                               , emit: mrna_fasta, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Pre-built mRNA reference (optional, improves accuracy)
    def mrna_ref = task.ext.mrna_ref ?: ""
    def mrna_option = mrna_ref ? "-m $mrna_ref" : ""

    """
    # Run FEELnc_codpot
    FEELnc_codpot.pl \\
        -i $candidate_fasta \\
        -a $genome_fasta \\
        $mrna_option \\
        --mode=shuffle \\
        --numtx=500 \\
        -o ${prefix}.feelnc_codpot.txt \\
        $args

    # Parse output to standard TSV format
    parse_feelnc_output.py \\
        --input ${prefix}.feelnc_codpot.txt \\
        --output ${prefix}.feelnc.tsv \\
        --fasta $candidate_fasta \\
        --prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        feelnc: \$(FEELnc_codpot.pl --version 2>&1 | grep -oP 'FEELnc version \\K[0-9.]+' || echo "0.2.1")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.feelnc.tsv
    touch ${prefix}.feelnc_codpot.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        feelnc: "0.2.1"
    END_VERSIONS
    """
}
