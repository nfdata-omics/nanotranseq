#!/usr/bin/env python3

"""
Parse FEELnc_codpot output to standardized TSV format
"""

import argparse
import sys
from Bio import SeqIO


def parse_feelnc_output(input_file, output_file, fasta_file, prefix):
    """
    Parse FEELnc output and create standardized TSV

    FEELnc output format:
    transcript_id   lncRNA_probability   mRNA_probability   classification
    """

    results = {}

    print(f"Parsing FEELnc output: {input_file}", file=sys.stderr)

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('transcript'):
                continue

            fields = line.strip().split('\t')
            if len(fields) >= 4:
                transcript_id = fields[0]
                lncrna_prob = float(fields[1])
                mrna_prob = float(fields[2])
                classification = fields[3]

                results[transcript_id] = {
                    'lncrna_prob': lncrna_prob,
                    'mrna_prob': mrna_prob,
                    'classification': classification
                }

    print(f"Parsed {len(results)} transcripts", file=sys.stderr)

    # Write standardized TSV
    with open(output_file, 'w') as out:
        out.write("transcript_id\tclassification\tscore\tlncrna_prob\tmrna_prob\n")

        for transcript_id, data in results.items():
            out.write(f"{transcript_id}\t{data['classification']}\t"
                      f"{data['lncrna_prob']}\t{data['lncrna_prob']}\t"
                      f"{data['mrna_prob']}\n")

    # Separate lncRNA and mRNA sequences
    lncrna_ids = {tid for tid, data in results.items()
                  if data['classification'] == 'lncRNA'}
    mrna_ids = {tid for tid, data in results.items()
                if data['classification'] == 'mRNA'}

    print(f"lncRNAs: {len(lncrna_ids)}, mRNAs: {len(mrna_ids)}", file=sys.stderr)

    # Write lncRNA FASTA
    if lncrna_ids:
        with open(f"{prefix}.lncRNA.fa", 'w') as lnc_out:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in lncrna_ids:
                    SeqIO.write(record, lnc_out, "fasta")

    # Write mRNA FASTA
    if mrna_ids:
        with open(f"{prefix}.mRNA.fa", 'w') as mrna_out:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in mrna_ids:
                    SeqIO.write(record, mrna_out, "fasta")

    print(f"✓ Created {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Parse FEELnc output to standardized TSV format'
    )

    parser.add_argument('--input', required=True,
                        help='FEELnc_codpot output file')
    parser.add_argument('--output', required=True,
                        help='Output TSV file')
    parser.add_argument('--fasta', required=True,
                        help='Input FASTA file (for splitting)')
    parser.add_argument('--prefix', required=True,
                        help='Output prefix')

    args = parser.parse_args()

    try:
        parse_feelnc_output(args.input, args.output, args.fasta, args.prefix)
    except Exception as e:
        print(f"✗ Error parsing FEELnc output: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
