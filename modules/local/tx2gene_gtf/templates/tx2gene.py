#!/usr/bin/env python
import re, gzip, platform


def op(f):
    return gzip.open(f, 'rt') if f.endswith('.gz') else open(f)


tx_re   = re.compile(r'transcript_id "([^"]+)"')
gene_re = re.compile(r'gene_id "([^"]+)"')
name_re = re.compile(r'gene_name "([^"]+)"')

seen = set()
with op("${gtf}") as fh, open("${prefix}.tx2gene.tsv", "w") as out:
    out.write("transcript_id\\tgene_id\\tgene_name\\n")
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 9 or cols[2] != "transcript":
            continue
        attrs = cols[8]
        tm, gm = tx_re.search(attrs), gene_re.search(attrs)
        if not tm or not gm:
            continue
        tx, gene = tm.group(1), gm.group(1)
        nm = name_re.search(attrs)
        name = nm.group(1) if nm else gene
        if tx in seen:
            continue
        seen.add(tx)
        out.write(f"{tx}\\t{gene}\\t{name}\\n")

with open("versions.yml", "w") as v:
    v.write('"${task.process}":\\n')
    v.write(f"    python: {platform.python_version()}\\n")
