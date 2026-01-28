#!/usr/bin/env python3

"""
Build CPAT hexamer table and logit model from training data
This is optional - you can use pre-built models for human/mouse
"""

import argparse
import sys
import subprocess


def build_hexamer_table(coding_fasta, noncoding_fasta, output_prefix):
    """Build hexamer frequency table"""

    print("Building hexamer table...", file=sys.stderr)

    cmd = [
        'make_hexamer_tab.py',
        '-c', coding_fasta,
        '-n', noncoding_fasta,
        '>', f'{output_prefix}_Hexamer.tsv'
    ]

    subprocess.run(' '.join(cmd), shell=True, check=True)

    print(f"✓ Created {output_prefix}_Hexamer.tsv", file=sys.stderr)


def build_logit_model(coding_fasta, noncoding_fasta, hexamer_table, output_prefix):
    """Build logistic regression model"""

    print("Building logit model...", file=sys.stderr)

    cmd = [
        'make_logitModel.py',
        '-c', coding_fasta,
        '-n', noncoding_fasta,
        '-x', hexamer_table,
        '-o', f'{output_prefix}_logitModel'
    ]

    subprocess.run(cmd, check=True)

    print(f"✓ Created {output_prefix}_logitModel.RData", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Build CPAT model from training data',
        epilog="""
Training data requirements:
  - coding_fasta: Known protein-coding transcripts (e.g., from GENCODE)
  - noncoding_fasta: Known lncRNAs (e.g., from lncRNAdb, NONCODE)

Example:
  build_cpat_model.py \\
    --coding human_mRNA.fa \\
    --noncoding human_lncRNA.fa \\
    --prefix Human
        """
    )

    parser.add_argument('--coding', required=True,
                        help='FASTA file with known coding transcripts')
    parser.add_argument('--noncoding', required=True,
                        help='FASTA file with known non-coding transcripts')
    parser.add_argument('--prefix', required=True,
                        help='Output prefix (e.g., Human, Mouse)')

    args = parser.parse_args()

    try:
        # Build hexamer table
        hexamer_file = f"{args.prefix}_Hexamer.tsv"
        build_hexamer_table(args.coding, args.noncoding, args.prefix)

        # Build logit model
        build_logit_model(args.coding, args.noncoding, hexamer_file, args.prefix)

        print("\n✓ CPAT model building complete!", file=sys.stderr)
        print(f"  - Hexamer table: {args.prefix}_Hexamer.tsv", file=sys.stderr)
        print(f"  - Logit model: {args.prefix}_logitModel.RData", file=sys.stderr)

    except Exception as e:
        print(f"✗ Error building CPAT model: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
