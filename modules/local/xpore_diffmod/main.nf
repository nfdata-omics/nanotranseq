process XPORE_DIFFMOD {
    tag "$meta.id"
    label 'process_medium'

    // ponytail: verify container tag against biocontainers registry before first run
    conda "bioconda::xpore=2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh5e36f6f_0' :
        'biocontainers/xpore:2.1--pyh5e36f6f_0' }"

    input:
    // samples: list of [ condition, id ] parallel to dataprep_dirs (dir basename = dataprep_<id>)
    tuple val(meta), val(samples), path(dataprep_dirs)

    output:
    tuple val(meta), path("diffmod/diffmod.table")                          , emit: diffmod
    tuple val(meta), path("diffmod/majority_direction_kmer_diffmod.table")  , emit: majority, optional: true
    path "*_xpore_mqc.tsv"                                                   , emit: mqc
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Build the xpore config.yml line by line: group replicate dataprep dirs by condition.
    // Each replicate points at the staged dir basename `dataprep_<id>`.
    def byCond = [:]
    samples.each { cond, id -> byCond.get(cond, []).add(id) }
    def cfg = ['data:']
    byCond.each { cond, ids ->
        cfg.add("  ${cond}:")
        ids.eachWithIndex { id, i -> cfg.add("    rep${i + 1}: dataprep_${id}") }
    }
    cfg.add('out: diffmod')
    def writeCfg = cfg.collect { line -> "echo '${line}' >> ${prefix}_config.yml" }.join('\n    ')
    """
    rm -f ${prefix}_config.yml
    ${writeCfg}

    xpore diffmod \\
        --config ${prefix}_config.yml \\
        --n_processes $task.cpus \\
        $args

    # Summary row for the MultiQC report. The comparison columns are named after the
    # conditions in the config, so the p-value column is found by prefix, not by name.
    python3 - <<'PYTHON' > ${prefix}_xpore_mqc.tsv
    import csv

    rows = list(csv.DictReader(open("diffmod/diffmod.table")))
    pvals = [c for c in (rows[0] if rows else {}) if c.startswith("pval_")]
    col = pvals[0] if pvals else None
    sig = [r for r in rows if col and r[col] not in ("", "nan") and float(r[col]) < 0.05]
    comparison = col[len("pval_"):] if col else "n/a"

    print("# id: 'xpore'")
    print("# section_name: 'xPore differential modification'")
    print("# description: 'Per-site comparison of modification rates between conditions. Significant sites are those with an uncorrected p-value below 0.05.'")
    print("# plot_type: 'table'")
    print("# pconfig:")
    print("#     id: 'xpore_table'")
    print("#     namespace: 'xpore'")
    print("Comparison\tSites tested\tSignificant sites (p < 0.05)")
    print("%s\t%d\t%d" % (comparison, len(rows), len(sig)))
    PYTHON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed -n 's/^xpore //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p diffmod
    echo -e "id,position,kmer,diff_mod_rate_ko_vs_wt,pval_ko_vs_wt,z_score_ko_vs_wt" > diffmod/diffmod.table
    touch diffmod/majority_direction_kmer_diffmod.table
    cat <<-MQC > ${prefix}_xpore_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: "2.1"
    END_VERSIONS
    """
}
