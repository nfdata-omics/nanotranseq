#!/usr/bin/env python3

"""
Convert PLEK output to standardized TSV format
"""

import argparse
import sys


def convert_plek_output(input_file, output_file):
    """
    Convert PLEK output to TSV with columns:
    transcript_id, classification, score

    PLEK output format:
    >transcript_id
    Non-coding (or Coding)
    Score: 0.123456
    """

    results = []

    with open(input_file, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('>'):
            # Transcript ID
            transcript_id = line[1:].split()[0]

            # Classification (next line)
            if i + 1 < len(lines):
                classification = lines[i + 1].strip()

                # Score (line after classification, if present)
                score = "NA"
                if i + 2 < len(lines) and lines[i + 2].strip().startswith('Score:'):
                    score_line = lines[i + 2].strip()
                    score = score_line.split(':')[1].strip()

                # Standardize classification
                if 'Non-coding' in classification or 'lncRNA' in classification:
                    pred_class = 'Non-coding'
                elif 'Coding' in classification or 'mRNA' in classification:
                    pred_class = 'Coding'
                else:
                    pred_class = classification

                results.append({
                    'transcript_id': transcript_id,
                    'classification': pred_class,
                    'score': score
                })

                i += 3  # Skip to next entry
            else:
                i += 1
        else:
            i += 1

    # Write TSV
    with open(output_file, 'w') as out:
        out.write("transcript_id\tclassification\tscore\n")
        for result in results:
            out.write(f"{result['transcript_id']}\t{result['classification']}\t{result['score']}\n")

    print(f"Converted {len(results)} PLEK predictions to {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Convert PLEK output to TSV format'
    )

    parser.add_argument('--input', required=True, help='PLEK output file')
    parser.add_argument('--output', required=True, help='Output TSV file')

    args = parser.parse_args()

    convert_plek_output(args.input, args.output)


if __name__ == '__main__':
    main()
