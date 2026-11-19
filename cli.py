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
from primers import (
    design_pcr_primers,
    reverse_translate,
)
from batch_analyzer import (
    analyze_batch_fasta,
)
from gc_skew_analyzer import (
    calculate_sliding_gc_skew,
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
    parser.add_argument(
        "--primers",
        action="store_true",
        help="Design forward and reverse flanking PCR primers with Tm calculation.",
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Perform high-throughput comparative analysis and distance matrix on a multi-FASTA file.",
    )
    parser.add_argument(
        "--skew",
        action="store_true",
        help="Calculate sliding-window GC skew and predict replication origin (oriC).",
    )

    args = parser.parse_args()

    raw_sequence = ""
    try:
        if args.batch:
            if not args.file:
                print("Error: --batch requires a file input (-f / --file).", file=sys.stderr)
                sys.exit(1)
            fpath = Path(args.file)
            if not fpath.is_file():
                print(f"Error: File not found: {args.file}", file=sys.stderr)
                sys.exit(1)
            b_data = analyze_batch_fasta(fpath.read_text(encoding="utf-8-sig"))
            if args.json:
                print(json.dumps(b_data, indent=2))
            else:
                print("=" * 70)
                print("GENEWEAVER MULTI-FASTA BATCH COMPARATIVE ANALYSIS")
                print("=" * 70)
                print(f"Total Sequences Analyzed : {b_data['total_sequences']}")
                print(f"Average Sequence Length  : {b_data['metrics']['avg_length_bp']} bp (Min: {b_data['metrics']['min_length_bp']}, Max: {b_data['metrics']['max_length_bp']})")
                print(f"Average GC Content       : {b_data['metrics']['avg_gc_content']}%")
                print("-" * 70)
                print(f"{'#':<3} {'Header':<30} {'Length':<8} {'GC%':<7} {'Tm(°C)':<7} {'ORFs':<5} {'Enz':<4}")
                print("-" * 70)
                for r in b_data["records"]:
                    hdr = (r["header"][:28] + "..") if len(r["header"]) > 30 else r["header"]
                    print(f"{r['index']:<3} {hdr:<30} {r['length']:<8} {r['gc_content']:<7} {r['melting_temperature']:<7} {r['orf_count']:<5} {r['restriction_cutters']:<4}")
                print("-" * 70)
                if b_data.get("closest_pair"):
                    cp = b_data["closest_pair"]
                    print(f"Closest Pair       : {cp['seq1']} <-> {cp['seq2']} ({cp['identity_pct']}% match, {cp['differences']} diffs)")
                if b_data.get("most_divergent_pair"):
                    dp = b_data["most_divergent_pair"]
                    print(f"Most Divergent Pair: {dp['seq1']} <-> {dp['seq2']} ({dp['identity_pct']}% match, {dp['differences']} diffs)")
                print("=" * 70)
            return

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

        if args.skew:
            skew_data = calculate_sliding_gc_skew(raw_sequence)
            if args.json:
                print(json.dumps(skew_data, indent=2))
            else:
                print("=" * 60)
                print("SLIDING-WINDOW GC SKEW & REPLICATION ORIGIN (oriC) ANALYSIS")
                print("=" * 60)
                print(f"Sequence Length        : {skew_data['sequence_length']} bp")
                print(f"Window Size / Step     : {skew_data['window_size']} bp / {skew_data['step_size']} bp")
                print(f"Total Sliding Windows  : {skew_data['total_windows']}")
                print(f"Predicted Origin (oriC): Position {skew_data['predicted_origin_bp']} bp (Min Cumulative Skew: {skew_data['min_cumulative_skew']})")
                print(f"Predicted Terminus(ter): Position {skew_data['predicted_terminus_bp']} bp (Max Cumulative Skew: {skew_data['max_cumulative_skew']})")
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

        if args.primers:
            p_data = design_pcr_primers(raw_sequence)
            fwd = p_data["forward_primer"]
            rev = p_data["reverse_primer"]
            print("=" * 60)
            print("PCR FLANKING PRIMER DESIGN & TM ANALYSIS")
            print("=" * 60)
            print(f"Template Length   : {p_data['template_length']} bp")
            print(f"Expected Amplicon : {p_data['amplicon_size']} bp")
            print(f"Tm Difference     : {p_data['tm_difference']} °C (Compatible: {p_data['compatible_annealing']})")
            print("-" * 60)
            print(f"Forward Primer (5' -> 3'): {fwd['sequence']}")
            print(f"  Length: {fwd['length']} bp | GC: {fwd['gc_content']}% | Tm: {fwd['melting_temperature']} °C | GC Clamp: {fwd['has_gc_clamp']} ({fwd['quality']})")
            print(f"Reverse Primer (5' -> 3'): {rev['sequence']}")
            print(f"  Length: {rev['length']} bp | GC: {rev['gc_content']}% | Tm: {rev['melting_temperature']} °C | GC Clamp: {rev['has_gc_clamp']} ({rev['quality']})")
            print("=" * 60)
            return

        if args.skew:
            skew_data = calculate_sliding_gc_skew(raw_sequence)
            if args.json:
                print(json.dumps(skew_data, indent=2))
            else:
                print("=" * 60)
                print("SLIDING-WINDOW GC SKEW & REPLICATION ORIGIN (oriC) ANALYSIS")
                print("=" * 60)
                print(f"Sequence Length        : {skew_data['sequence_length']} bp")
                print(f"Window Size / Step     : {skew_data['window_size']} bp / {skew_data['step_size']} bp")
                print(f"Total Sliding Windows  : {skew_data['total_windows']}")
                print(f"Predicted Origin (oriC): Position {skew_data['predicted_origin_bp']} bp (Min Cumulative Skew: {skew_data['min_cumulative_skew']})")
                print(f"Predicted Terminus(ter): Position {skew_data['predicted_terminus_bp']} bp (Max Cumulative Skew: {skew_data['max_cumulative_skew']})")
                print("=" * 60)
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
