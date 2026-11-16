"""
Batch FASTA sequence analysis and pairwise genetic distance matrix computation.
Enables high-throughput multi-sequence comparisons and summary statistics.
"""

from typing import Dict, List, Any, Optional
from analyzer import (
    clean_sequence,
    parse_fasta,
    calculate_molecular_weight,
    calculate_melting_temperature,
    find_orfs,
    find_restriction_sites,
    SequenceError,
)


def compute_sequence_distance(seq1: str, seq2: str) -> Dict[str, Any]:
    """
    Compute Hamming distance and percentage identity between two aligned sequences.
    """
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return {"distance": 0, "identity_pct": 100.0}

    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    total_differences = mismatches + (max_len - min_len)
    identity_pct = round(((max_len - total_differences) / max_len) * 100.0, 2)

    return {
        "mismatches": total_differences,
        "identity_pct": max(0.0, identity_pct),
    }


def analyze_batch_fasta(fasta_content: str) -> Dict[str, Any]:
    """
    Parse a multi-record FASTA string and compute comparative bioinformatics statistics
    along with an N x N pairwise genetic distance matrix.
    """
    records = parse_fasta(fasta_content)
    if not records:
        raise SequenceError("No valid sequence records found in FASTA data.")

    summaries: List[Dict[str, Any]] = []
    lengths: List[int] = []
    gc_contents: List[float] = []

    for idx, rec in enumerate(records, 1):
        header = rec["header"]
        raw_seq = rec["sequence"]
        seq = clean_sequence(raw_seq)
        seq_len = len(seq)

        g_count = seq.count("G")
        c_count = seq.count("C")
        gc_pct = round(((g_count + c_count) / seq_len) * 100.0, 2) if seq_len > 0 else 0.0

        mw = calculate_molecular_weight(seq)
        tm = calculate_melting_temperature(seq)
        orfs = find_orfs(seq, min_codon_length=2)
        restriction = find_restriction_sites(seq)
        cutters_count = sum(1 for e in restriction if e["cuts_sequence"])

        lengths.append(seq_len)
        gc_contents.append(gc_pct)

        summaries.append({
            "index": idx,
            "header": header,
            "sequence": seq,
            "length": seq_len,
            "gc_content": gc_pct,
            "at_content": round(100.0 - gc_pct, 2),
            "melting_temperature": tm,
            "molecular_weight_ds": mw["dsDNA"],
            "orf_count": len(orfs),
            "restriction_cutters": cutters_count,
        })

    total_records = len(summaries)
    avg_length = round(sum(lengths) / total_records, 1) if total_records else 0
    avg_gc = round(sum(gc_contents) / total_records, 2) if total_records else 0.0

    # Build N x N pairwise distance matrix
    headers = [s["header"] for s in summaries]
    matrix: List[List[float]] = []
    closest_pair: Optional[Dict[str, Any]] = None
    most_divergent: Optional[Dict[str, Any]] = None

    for i in range(total_records):
        row: List[float] = []
        for j in range(total_records):
            if i == j:
                row.append(100.0)
            else:
                dist_data = compute_sequence_distance(summaries[i]["sequence"], summaries[j]["sequence"])
                ident = dist_data["identity_pct"]
                row.append(ident)

                # Track closest and most divergent pairs for i < j
                if i < j:
                    pair_info = {
                        "seq1": headers[i],
                        "seq2": headers[j],
                        "identity_pct": ident,
                        "differences": dist_data["mismatches"],
                    }
                    if closest_pair is None or ident > closest_pair["identity_pct"]:
                        closest_pair = pair_info
                    if most_divergent is None or ident < most_divergent["identity_pct"]:
                        most_divergent = pair_info
        matrix.append(row)

    return {
        "total_sequences": total_records,
        "metrics": {
            "avg_length_bp": avg_length,
            "min_length_bp": min(lengths),
            "max_length_bp": max(lengths),
            "avg_gc_content": avg_gc,
        },
        "records": summaries,
        "headers": headers,
        "distance_matrix": matrix,
        "closest_pair": closest_pair,
        "most_divergent_pair": most_divergent,
    }
