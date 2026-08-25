process M6ANET {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::m6anet=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/m6anet:2.1.0--pyhdfd78af_0' :
        'biocontainers/m6anet:2.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(eventalign)   // NANOPOLISH_EVENTALIGN.out.eventalign

    output:
    tuple val(meta), path("*.m6a_sites.csv"), emit: sites
    path "*_m6anet_mqc.tsv"                 , emit: mqc
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    m6anet dataprep \\
        --eventalign $eventalign \\
        --out_dir dataprep \\
        --n_processes $task.cpus \\
        $args

    m6anet inference \\
        --input_dir dataprep \\
        --out_dir inference \\
        --n_processes $task.cpus \\
        $args2

    cp inference/data.site_proba.csv ${prefix}.m6a_sites.csv

    # One row per sample for the MultiQC report. Files sharing the `id` below are merged
    # into a single table, so every sample lands in the same section.
    # Tabs are written as \t: Groovy turns them into real tabs before the shell sees this.
    python3 - <<'PYTHON' > ${prefix}_m6anet_mqc.tsv
    import csv

    rows = list(csv.DictReader(open("${prefix}.m6a_sites.csv")))
    high = [r for r in rows if float(r["probability_modified"]) > 0.9]
    ratio = sum(float(r["mod_ratio"]) for r in high) / len(high) if high else 0.0

    print("# id: 'm6anet'")
    print("# section_name: 'm6anet m6A sites'")
    print("# description: 'DRACH sites called by m6anet on the nanopolish eventalign table. High-confidence sites are those with probability_modified above 0.9.'")
    print("# plot_type: 'table'")
    print("# pconfig:")
    print("#     id: 'm6anet_table'")
    print("#     namespace: 'm6anet'")
    print("Sample\tSites tested\tHigh-confidence sites\tMean modification ratio")
    print("${prefix}\t%d\t%d\t%.3f" % (len(rows), len(high), ratio))
    PYTHON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        m6anet: \$(m6anet --version 2>&1 | sed -n 's/^m6anet //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio" > ${prefix}.m6a_sites.csv

    cat <<-MQC > ${prefix}_m6anet_mqc.tsv
    # id: 'm6anet'
    # plot_type: 'table'
    Sample\tSites tested\tHigh-confidence sites\tMean modification ratio
    ${prefix}\t0\t0\t0.000
    MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        m6anet: "2.1.0"
    END_VERSIONS
    """
}
