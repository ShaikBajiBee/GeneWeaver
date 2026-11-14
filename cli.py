#!/usr/bin/env python3
"""
GeneWeaver Command-Line Interface (CLI).
Analyze DNA sequences or FASTA files directly from the terminal.
"""

import argparse
import json
import sys
from pathlib import Path
from analyzer import (
    analyze_sequence,
    parse_fasta,
    generate_text_report,
    SequenceError,
)


def main():
    parser = argparse.ArgumentParser(
        description="GeneWeaver CLI - Analyze DNA sequences, compute stats, and detect ORFs.",
        epilog="Example: python cli.py -s ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-s", "--sequence",
        help="Raw DNA sequence string to analyze.",
    )
    group.add_argument(
        "-f", "--file",
        help="Path to FASTA or text file containing DNA sequence.",
    )

    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results in JSON format.",
    )
    parser.add_argument(
        "-o", "--output",
        help="Save output report to specified file path.",
    )
    parser.add_argument(
        "--orfs",
        action="store_true",
        help="Display detected Open Reading Frames exclusively.",
    )

    args = parser.parse_args()

    raw_sequence = ""
    try:
        if args.sequence:
            raw_sequence = args.sequence
        elif args.file:
            path = Path(args.file)
            if not path.is_file():
                print(f"Error: File not found: {args.file}", file=sys.stderr)
                sys.exit(1)
            content = path.read_text(encoding="utf-8-sig")
            records = parse_fasta(content)
            if not records:
                print("Error: No sequence data found in file.", file=sys.stderr)
                sys.exit(1)
            raw_sequence = records[0]["sequence"]

        data = analyze_sequence(raw_sequence)

        if args.orfs:
            print(f"Detected {len(data['orfs'])} Open Reading Frame(s):")
            for idx, orf in enumerate(data["orfs"], 1):
                print(f"  [{idx}] Frame +{orf['frame']} | {orf['start']}-{orf['end']} "
                      f"({orf['length_bp']} bp, {orf['codons']} aa) -> {orf['protein']}")
            return

        if args.json:
            output_text = json.dumps(data, indent=2)
        else:
            output_text = generate_text_report(data)

        if args.output:
            out_path = Path(args.output)
            out_path.write_text(output_text, encoding="utf-8")
            print(f"Report saved successfully to {out_path}")
        else:
            print(output_text)

    except SequenceError as exc:
        print(f"Validation Error: {exc}", file=sys.stderr)
        sys.exit(2)
    except Exception as exc:
        print(f"Execution Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
