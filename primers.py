"""
Reverse-translation and PCR primer design utilities for GeneWeaver.
Provides degenerate codon consensus mapping and forward/reverse PCR primer analysis.
"""

from typing import Dict, List, Any, Tuple
import math

# Reverse genetic code mapping amino acid -> (consensus degenerate codon, all possible codons)
REVERSE_CODON_TABLE: Dict[str, Tuple[str, List[str]]] = {
    "A": ("GCN", ["GCT", "GCC", "GCA", "GCG"]),
    "C": ("TGY", ["TGT", "TGC"]),
    "D": ("GAY", ["GAT", "GAC"]),
    "E": ("GAR", ["GAA", "GAG"]),
    "F": ("TTY", ["TTT", "TTC"]),
    "G": ("GGN", ["GGT", "GGC", "GGA", "GGG"]),
    "H": ("CAY", ["CAT", "CAC"]),
    "I": ("ATH", ["ATT", "ATC", "ATA"]),
    "K": ("AAR", ["AAA", "AAG"]),
    "L": ("YTN", ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
    "M": ("ATG", ["ATG"]),
    "N": ("AAY", ["AAT", "AAC"]),
    "P": ("CCN", ["CCT", "CCC", "CCA", "CCG"]),
    "Q": ("CAR", ["CAA", "CAG"]),
    "R": ("MGN", ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
    "S": ("WSN", ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
    "T": ("ACN", ["ACT", "ACC", "ACA", "ACG"]),
    "V": ("GTN", ["GTT", "GTC", "GTA", "GTG"]),
    "W": ("TGG", ["TGG"]),
    "Y": ("TAY", ["TAT", "TAC"]),
    "*": ("TRR", ["TAA", "TAG", "TGA"]),
}

COMPLEMENT_MAP = str.maketrans("ATGC", "TACG")


def reverse_translate(peptide: str) -> Dict[str, Any]:
    """
    Reverse-translate an amino acid sequence into a degenerate DNA consensus sequence.
    Calculates total combinatorial possibilities (degeneracy fold).
    """
    clean_pep = "".join(peptide.split()).upper()
    if not clean_pep:
        raise ValueError("Peptide sequence cannot be empty.")

    consensus_codons = []
    total_combinations = 1
    breakdown = []

    for aa in clean_pep:
        if aa not in REVERSE_CODON_TABLE:
            raise ValueError(f"Unknown amino acid code: '{aa}'. Please use standard 1-letter IUPAC codes.")
        codon_deg, all_codons = REVERSE_CODON_TABLE[aa]
        consensus_codons.append(codon_deg)
        total_combinations *= len(all_codons)
        breakdown.append({
            "amino_acid": aa,
            "consensus": codon_deg,
            "degeneracy": len(all_codons),
            "possible_codons": all_codons,
        })

    consensus_dna = "".join(consensus_codons)

    return {
        "peptide": clean_pep,
        "peptide_length": len(clean_pep),
        "consensus_dna": consensus_dna,
        "dna_length": len(consensus_dna),
        "total_combinations": total_combinations,
        "breakdown": breakdown,
    }


def analyze_primer(primer_seq: str) -> Dict[str, Any]:
    """
    Analyze biophysical properties of a single primer oligonucleotide:
    Length, GC content, Melting Temperature (Tm), and 3' GC Clamp.
    """
    seq = "".join(primer_seq.split()).upper()
    length = len(seq)
    if length == 0:
        return {}

    g_count = seq.count("G")
    c_count = seq.count("C")
    gc_pct = round(((g_count + c_count) / length) * 100.0, 1)

    # Nearest-neighbor salt-adjusted formula for primers
    if length < 14:
        tm = 2.0 * (seq.count("A") + seq.count("T")) + 4.0 * (g_count + c_count)
    else:
        tm = 64.9 + 41.0 * ((g_count + c_count - 16.4) / length)

    # GC clamp: does the 3' end terminate with G or C? (Improves polymerase binding efficiency)
    has_gc_clamp = seq[-1] in ("G", "C") if length > 0 else False

    # Ideal criteria: Length 18-24 bp, GC 40-60%, Tm 52-65 °C, GC clamp present
    quality = "Optimal"
    notes = []
    if length < 18 or length > 28:
        notes.append(f"Length ({length} bp) outside typical 18-28 bp range")
    if gc_pct < 40 or gc_pct > 65:
        notes.append(f"GC content ({gc_pct}%) outside ideal 40-65% window")
    if tm < 50 or tm > 68:
        notes.append(f"Tm ({round(tm, 1)} °C) outside ideal 50-68 °C window")
    if not has_gc_clamp:
        notes.append("No 3' GC clamp (3' base is A or T)")

    if len(notes) > 1:
        quality = "Suboptimal"
    elif len(notes) == 1:
        quality = "Acceptable"

    return {
        "sequence": seq,
        "length": length,
        "gc_content": gc_pct,
        "melting_temperature": round(tm, 1),
        "has_gc_clamp": has_gc_clamp,
        "quality": quality,
        "notes": notes,
    }


def design_pcr_primers(template_dna: str, primer_len: int = 20) -> Dict[str, Any]:
    """
    Design forward and reverse flanking PCR primers for a target DNA template sequence.
    - Forward Primer: 5' end of sense strand (5' -> 3')
    - Reverse Primer: 5' -> 3' reverse complement of the 3' end of sense strand
    """
    seq = "".join(template_dna.split()).upper()
    if len(seq) < primer_len * 2:
        raise ValueError(f"Sequence length ({len(seq)} bp) is too short to design paired primers of length {primer_len} bp.")

    fwd_seq = seq[:primer_len]
    rev_raw = seq[-primer_len:]
    rev_seq = rev_raw.translate(COMPLEMENT_MAP)[::-1]

    fwd_info = analyze_primer(fwd_seq)
    rev_info = analyze_primer(rev_seq)

    # Tm difference between primers (ideal is < 3-5 °C)
    tm_diff = round(abs(fwd_info["melting_temperature"] - rev_info["melting_temperature"]), 1)
    amplicon_length = len(seq)

    return {
        "template_length": amplicon_length,
        "amplicon_size": amplicon_length,
        "forward_primer": fwd_info,
        "reverse_primer": rev_info,
        "tm_difference": tm_diff,
        "compatible_annealing": tm_diff <= 5.0,
    }
