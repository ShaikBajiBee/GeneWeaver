#!/usr/bin/env python3
"""
GeneWeaver Command-Line Interface (CLI).
Analyze DNA sequences, detect ORFs & restriction sites, and compare sequence variants.
"""

import argparse
import json
import sys
from pathlib import Path
from analyzer import (
    analyze_sequence,
    parse_fasta,
    generate_text_report,
    compare_sequences,
    SequenceError,
)


def main():
    parser = argparse.ArgumentParser(
        description="GeneWeaver CLI - Analyze DNA sequences, compute stats, and compare sequence variants.",
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
        "--compare",
        help="Query sequence to align and compare against the reference sequence for SNP/variant detection.",
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
    parser.add_argument(
        "--enzymes",
        action="store_true",
        help="Display restriction enzyme digestion profile exclusively.",
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

        # If --compare option is provided
        if args.compare:
            cmp_data = compare_sequences(raw_sequence, args.compare)
            if args.json:
                print(json.dumps(cmp_data, indent=2))
            else:
                print("=" * 60)
                print("GENEWEAVER SEQUENCE VARIANT COMPARISON REPORT")
                print("=" * 60)
                print(f"Reference Length : {cmp_data['ref_length']} bp")
                print(f"Query Length     : {cmp_data['query_length']} bp")
                print(f"Sequence Identity: {cmp_data['identity_pct']}%")
                print(f"Hamming Distance : {cmp_data['hamming_distance']}")
                print(f"Total Mutations  : {cmp_data['mutations_count']}")
                print(f"  Transitions (Ts): {cmp_data['transitions']}")
                print(f"  Transversions (Tv): {cmp_data['transversions']}")
                print(f"  Ts/Tv Ratio     : {cmp_data['ts_tv_ratio']}")
                print(f"  Synonymous (Silent): {cmp_data['synonymous_count']}")
                print(f"  Missense (Amino acid change): {cmp_data['missense_count']}")
                print(f"  Nonsense (Stop codon): {cmp_data['nonsense_count']}")
                print("-" * 60)
                print("MUTATION DETAILS:")
                for m in cmp_data["mutations"]:
                    print(f"  Pos {m['position']}: {m['notation']} [{m['type']}] -> {m['impact']}")
                print("=" * 60)
            return

        data = analyze_sequence(raw_sequence)

        if args.orfs:
            print(f"Detected {len(data['orfs'])} Open Reading Frame(s):")
            for idx, orf in enumerate(data["orfs"], 1):
                print(f"  [{idx}] Frame +{orf['frame']} | {orf['start']}-{orf['end']} "
                      f"({orf['length_bp']} bp, {orf['codons']} aa) -> {orf['protein']}")
            return

        if args.enzymes:
            print("Restriction Enzyme Digestion Profile:")
            cutters = [e for e in data["restriction_sites"] if e["cuts_sequence"]]
            if not cutters:
                print("  No recognition sites found for standard enzymes.")
            else:
                for enz in cutters:
                    cuts = ", ".join(str(p) for p in enz["cut_positions"])
                    frags = ", ".join(f"{f}bp" for f in enz["fragments"])
                    print(f"  * {enz['name']} ({enz['site']}): {enz['cuts']} cut(s) at [{cuts}] -> Fragments: {frags}")
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
