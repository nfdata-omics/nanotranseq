process TXIMPORT {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta%3A1.24.0--r44hdfd78af_0' :
        'biocontainers/bioconductor-tximeta:1.24.0--r44hdfd78af_0' }"

    input:
    tuple val(meta), path("quants/*")
    tuple val(meta2), path(tx2gene)
    val quant_type   // 'oarfish' (also handles 'salmon'/'kallisto' for reuse)

    output:
    tuple val(meta), path("*gene_tpm.tsv")                 , emit: tpm_gene
    tuple val(meta), path("*gene_counts.tsv")              , emit: counts_gene
    tuple val(meta), path("*gene_counts_length_scaled.tsv"), emit: counts_gene_length_scaled
    tuple val(meta), path("*gene_counts_scaled.tsv")       , emit: counts_gene_scaled
    tuple val(meta), path("*gene_lengths.tsv")             , emit: lengths_gene
    tuple val(meta), path("*transcript_tpm.tsv")           , emit: tpm_transcript
    tuple val(meta), path("*transcript_counts.tsv")        , emit: counts_transcript
    tuple val(meta), path("*transcript_lengths.tsv")       , emit: lengths_transcript
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'tximport.r'

    stub:
    """
    touch ${meta.id}.gene_tpm.tsv
    touch ${meta.id}.gene_counts.tsv
    touch ${meta.id}.gene_counts_length_scaled.tsv
    touch ${meta.id}.gene_counts_scaled.tsv
    touch ${meta.id}.gene_lengths.tsv
    touch ${meta.id}.transcript_tpm.tsv
    touch ${meta.id}.transcript_counts.tsv
    touch ${meta.id}.transcript_lengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
