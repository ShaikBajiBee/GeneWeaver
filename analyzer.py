"""DNA sequence validation, bioinformatics analysis, and transformation utilities."""

from typing import Dict, List, Optional, Tuple, Any


class SequenceError(ValueError):
    """Raised when a DNA sequence is empty or contains invalid characters."""


CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AMINO_ACID_NAMES: Dict[str, str] = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic Acid",
    "C": "Cysteine", "E": "Glutamic Acid", "Q": "Glutamine", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
    "*": "Stop",
}

VALID_BASES = set("ATGC")
COMPLEMENT = str.maketrans("ATGC", "TACG")

# Molecular weights in g/mol (anhydrous nucleotide monophosphates)
NUCLEOTIDE_MW: Dict[str, float] = {
    "A": 313.21,
    "T": 304.20,
    "C": 289.18,
    "G": 329.21,
}

AMINO_ACID_MW: Dict[str, float] = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10,
    "C": 121.16, "E": 147.13, "Q": 146.15, "G": 75.07,
    "H": 155.16, "I": 131.18, "L": 131.18, "K": 146.19,
    "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
    "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}

# Standard Restriction Endonucleases for cloning & mapping
RESTRICTION_ENZYMES: Dict[str, Dict[str, Any]] = {
    "EcoRI": {"site": "GAATTC", "cut_offset": 1, "source": "Escherichia coli"},
    "BamHI": {"site": "GGATCC", "cut_offset": 1, "source": "Bacillus amyloliquefaciens"},
    "HindIII": {"site": "AAGCTT", "cut_offset": 1, "source": "Haemophilus influenzae"},
    "NotI": {"site": "GCGGCCGC", "cut_offset": 2, "source": "Nocardia otitidis"},
    "PstI": {"site": "CTGCAG", "cut_offset": 5, "source": "Providencia stuartii"},
    "SmaI": {"site": "CCCGGG", "cut_offset": 3, "source": "Serratia marcescens"},
    "XhoI": {"site": "CTCGAG", "cut_offset": 1, "source": "Xanthomonas holcicola"},
    "TaqI": {"site": "TCGA", "cut_offset": 1, "source": "Thermus aquaticus"},
    "HaeIII": {"site": "GGCC", "cut_offset": 2, "source": "Haemophilus aegyptius"},
    "AluI": {"site": "AGCT", "cut_offset": 2, "source": "Arthrobacter luteus"},
}

# pKa values for peptide charge & isoelectric point (Lehninger scale)
PKA_VALUES = {
    "N_term": 9.69,
    "C_term": 2.34,
    "D": 3.86,  # Asp
    "E": 4.25,  # Glu
    "C": 8.33,  # Cys
    "Y": 10.07, # Tyr
    "H": 6.00,  # His
    "K": 10.53, # Lys
    "R": 12.48, # Arg
}


def clean_sequence(raw_sequence: str) -> str:
    """Remove whitespace and line breaks; normalize to uppercase DNA sequence."""
    if not raw_sequence:
        raise SequenceError("Please enter a DNA sequence or upload a FASTA file.")

    sequence = "".join(raw_sequence.split()).upper()
    if not sequence:
        raise SequenceError("The sequence contains only whitespace.")

    invalid = sorted(set(sequence) - VALID_BASES)
    if invalid:
        shown = ", ".join(invalid)
        raise SequenceError(
            f"Invalid DNA character(s): {shown}. Only standard bases (A, T, G, C) are accepted."
        )
    return sequence


def parse_fasta(contents: str) -> List[Dict[str, str]]:
    """
    Parse FASTA content supporting single or multiple records.
    Returns a list of dicts: [{'header': ..., 'sequence': ...}].
    """
    lines = [line.strip() for line in contents.splitlines() if line.strip()]
    if not lines:
        raise SequenceError("The FASTA file is empty.")

    # Plain text without FASTA header
    if not any(line.startswith(">") for line in lines):
        return [{"header": "Sequence 1", "sequence": "".join(lines)}]

    records: List[Dict[str, str]] = []
    current_header = ""
    current_seq: List[str] = []

    for line in lines:
        if line.startswith(">"):
            if current_header:
                seq_str = "".join(current_seq)
                if seq_str:
                    records.append({"header": current_header, "sequence": seq_str})
                current_seq = []
            current_header = line[1:].strip() or "Unnamed record"
        else:
            current_seq.append(line)

    if current_header and current_seq:
        records.append({"header": current_header, "sequence": "".join(current_seq)})

    if not records:
        raise SequenceError("No valid sequence data found after FASTA header(s).")

    return records


def reverse_complement(sequence: str) -> str:
    """Return reverse complement (5' -> 3' of the opposite strand)."""
    return sequence.translate(COMPLEMENT)[::-1]


def transcribe(sequence: str) -> str:
    """Transcribe DNA into mRNA (replace Thymine with Uracil)."""
    return sequence.replace("T", "U")


def translate(sequence: str) -> str:
    """Translate complete codons starting from index 0 until first stop codon."""
    protein: List[str] = []
    for index in range(0, len(sequence) - 2, 3):
        codon = sequence[index:index + 3]
        amino_acid = CODON_TABLE.get(codon, "?")
        if amino_acid == "*":
            protein.append("Stop")
            break
        protein.append(amino_acid)
    return " ".join(protein) if protein else "No complete codon found"


def calculate_molecular_weight(sequence: str) -> Dict[str, float]:
    """
    Calculate estimated molecular weights:
    - Single-stranded DNA (ssDNA in g/mol)
    - Double-stranded DNA (dsDNA in g/mol)
    """
    length = len(sequence)
    if length == 0:
        return {"ssDNA": 0.0, "dsDNA": 0.0}

    base_sum = sum(NUCLEOTIDE_MW[b] for b in sequence)
    ss_dna = base_sum - 61.96 if length > 1 else base_sum

    comp_sum = sum(NUCLEOTIDE_MW[b.translate(COMPLEMENT)] for b in sequence)
    comp_dna = comp_sum - 61.96 if length > 1 else comp_sum
    ds_dna = ss_dna + comp_dna

    return {
        "ssDNA": round(ss_dna, 2),
        "dsDNA": round(ds_dna, 2),
    }


def calculate_melting_temperature(sequence: str) -> float:
    """
    Estimate melting temperature (Tm in degrees Celsius):
    - Wallace rule for short sequences (< 14 bases): Tm = 2(A+T) + 4(G+C)
    - Salt-adjusted empirical rule for >= 14 bases: Tm = 64.9 + 41 * (yG+zC - 16.4) / N
    """
    length = len(sequence)
    if length == 0:
        return 0.0

    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")

    if length < 14:
        tm = 2 * (a_count + t_count) + 4 * (g_count + c_count)
    else:
        tm = 64.9 + 41.0 * ((g_count + c_count - 16.4) / length)

    return round(tm, 2)


def calculate_gc_skew(sequence: str) -> float:
    """Calculate GC Skew: (G - C) / (G + C). Range: -1.0 to 1.0."""
    g = sequence.count("G")
    c = sequence.count("C")
    if g + c == 0:
        return 0.0
    return round((g - c) / (g + c), 4)


def find_orfs(sequence: str, min_codon_length: int = 3) -> List[Dict[str, Any]]:
    """
    Locate Open Reading Frames (ORFs) across 3 forward reading frames (Frame 1, 2, 3).
    An ORF starts with start codon 'ATG' and terminates at 'TAA', 'TAG', or 'TGA'.
    """
    orfs: List[Dict[str, Any]] = []
    seq_len = len(sequence)

    for frame in range(3):
        i = frame
        while i < seq_len - 2:
            codon = sequence[i:i + 3]
            if codon == "ATG":
                start_pos = i
                orf_codons = ["M"]
                found_stop = False
                j = start_pos + 3

                while j < seq_len - 2:
                    curr_codon = sequence[j:j + 3]
                    aa = CODON_TABLE.get(curr_codon, "?")
                    if aa == "*":
                        found_stop = True
                        stop_pos = j + 3
                        orf_seq = sequence[start_pos:stop_pos]
                        if len(orf_codons) >= min_codon_length:
                            orfs.append({
                                "frame": frame + 1,
                                "start": start_pos + 1,
                                "end": stop_pos,
                                "length_bp": len(orf_seq),
                                "codons": len(orf_seq) // 3,
                                "dna": orf_seq,
                                "protein": "".join(orf_codons),
                            })
                        i = j
                        break
                    orf_codons.append(aa)
                    j += 3

                if not found_stop:
                    pass
            i += 3

    return orfs


def find_restriction_sites(sequence: str) -> List[Dict[str, Any]]:
    """
    Find recognition sites and compute digestion fragments for standard restriction endonucleases.
    """
    results: List[Dict[str, Any]] = []
    seq_len = len(sequence)

    for enzyme, info in RESTRICTION_ENZYMES.items():
        site = info["site"]
        offset = info["cut_offset"]
        cut_positions: List[int] = []

        # Find all occurrences of the recognition site
        start = 0
        while True:
            idx = sequence.find(site, start)
            if idx == -1:
                break
            cut_pos = idx + offset
            cut_positions.append(cut_pos)
            start = idx + 1

        num_cuts = len(cut_positions)

        # Calculate resulting fragment sizes for linear DNA
        fragments: List[int] = []
        if num_cuts == 0:
            fragments = [seq_len]
        else:
            prev = 0
            for pos in cut_positions:
                fragments.append(pos - prev)
                prev = pos
            fragments.append(seq_len - prev)

        fragments.sort(reverse=True)

        results.append({
            "name": enzyme,
            "site": site,
            "source": info["source"],
            "cuts": num_cuts,
            "cut_positions": cut_positions,
            "fragments": fragments,
            "cuts_sequence": num_cuts > 0,
        })

    # Sort enzymes by number of cuts (descending)
    results.sort(key=lambda x: x["cuts"], reverse=True)
    return results


def calculate_peptide_charge(sequence_aa: str, ph: float) -> float:
    """Compute the net electrical charge of an amino acid chain at a given pH."""
    # N-terminal positive charge
    charge = 1.0 / (1.0 + 10 ** (ph - PKA_VALUES["N_term"]))
    # C-terminal negative charge
    charge -= 1.0 / (1.0 + 10 ** (PKA_VALUES["C_term"] - ph))

    # Basic side-chains (positive when protonated)
    charge += sequence_aa.count("K") / (1.0 + 10 ** (ph - PKA_VALUES["K"]))
    charge += sequence_aa.count("R") / (1.0 + 10 ** (ph - PKA_VALUES["R"]))
    charge += sequence_aa.count("H") / (1.0 + 10 ** (ph - PKA_VALUES["H"]))

    # Acidic side-chains (negative when deprotonated)
    charge -= sequence_aa.count("D") / (1.0 + 10 ** (PKA_VALUES["D"] - ph))
    charge -= sequence_aa.count("E") / (1.0 + 10 ** (PKA_VALUES["E"] - ph))
    charge -= sequence_aa.count("C") / (1.0 + 10 ** (PKA_VALUES["C"] - ph))
    charge -= sequence_aa.count("Y") / (1.0 + 10 ** (PKA_VALUES["Y"] - ph))

    return charge


def calculate_isoelectric_point(protein_str: str) -> Dict[str, float]:
    """
    Determine the isoelectric point (pI) where net charge = 0 using bisection method.
    Also calculates net charge at physiological pH (7.4).
    """
    clean_aa = "".join(aa for aa in protein_str.split() if aa in AMINO_ACID_MW)
    if not clean_aa:
        return {"pi": 7.0, "charge_ph74": 0.0}

    # Bisection search between pH 0.0 and 14.0
    low = 0.0
    high = 14.0
    pi = 7.0

    for _ in range(30):
        mid = (low + high) / 2.0
        charge = calculate_peptide_charge(clean_aa, mid)
        if abs(charge) < 0.001:
            pi = mid
            break
        elif charge > 0:
            low = mid
        else:
            high = mid
        pi = mid

    charge_74 = calculate_peptide_charge(clean_aa, 7.4)

    return {
        "pi": round(pi, 2),
        "charge_ph74": round(charge_74, 2),
    }


def compare_sequences(ref_raw: str, query_raw: str) -> Dict[str, Any]:
    """
    Compare a query/mutant DNA sequence against a reference sequence.
    Detects single nucleotide polymorphisms (SNPs), transitions, transversions,
    and determines functional amino acid impact (synonymous, missense, nonsense).
    """
    ref = clean_sequence(ref_raw)
    query = clean_sequence(query_raw)

    min_len = min(len(ref), len(query))
    max_len = max(len(ref), len(query))
    length_diff = abs(len(ref) - len(query))

    purines = {"A", "G"}
    pyrimidines = {"C", "T"}

    mutations: List[Dict[str, Any]] = []
    transitions = 0
    transversions = 0
    synonymous = 0
    missense = 0
    nonsense = 0

    for i in range(min_len):
        r_base = ref[i]
        q_base = query[i]

        if r_base != q_base:
            is_transition = (r_base in purines and q_base in purines) or (
                r_base in pyrimidines and q_base in pyrimidines
            )
            mut_type = "Transition" if is_transition else "Transversion"
            if is_transition:
                transitions += 1
            else:
                transversions += 1

            codon_start = (i // 3) * 3
            ref_codon = ref[codon_start : codon_start + 3] if codon_start + 3 <= len(ref) else ""
            alt_codon = query[codon_start : codon_start + 3] if codon_start + 3 <= len(query) else ""

            ref_aa = CODON_TABLE.get(ref_codon, "-") if len(ref_codon) == 3 else "-"
            alt_aa = CODON_TABLE.get(alt_codon, "-") if len(alt_codon) == 3 else "-"

            if ref_aa != "-" and alt_aa != "-":
                if ref_aa == alt_aa:
                    impact = "Synonymous (Silent)"
                    synonymous += 1
                elif alt_aa == "*":
                    impact = "Nonsense (Stop gained)"
                    nonsense += 1
                else:
                    impact = f"Missense ({ref_aa} -> {alt_aa})"
                    missense += 1
            else:
                impact = "Non-coding / Partial codon"

            mutations.append({
                "position": i + 1,
                "ref_base": r_base,
                "alt_base": q_base,
                "notation": f"c.{i + 1}{r_base}>{q_base}",
                "type": mut_type,
                "ref_codon": ref_codon,
                "alt_codon": alt_codon,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "impact": impact,
            })

    matches = sum(1 for i in range(min_len) if ref[i] == query[i])
    identity_pct = round((matches / max_len) * 100.0, 2) if max_len > 0 else 100.0
    hamming_distance = len(mutations) + length_diff
    ts_tv_ratio = round(transitions / transversions, 2) if transversions > 0 else float(transitions)

    return {
        "ref_length": len(ref),
        "query_length": len(query),
        "length_matches": len(ref) == len(query),
        "identity_pct": identity_pct,
        "hamming_distance": hamming_distance,
        "mutations_count": len(mutations),
        "transitions": transitions,
        "transversions": transversions,
        "ts_tv_ratio": ts_tv_ratio,
        "synonymous_count": synonymous,
        "missense_count": missense,
        "nonsense_count": nonsense,
        "mutations": mutations,
    }



def generate_text_report(data: Dict[str, Any]) -> str:
    """Generate a clean, formatted plain-text bioinformatics report."""
    seq = data["sequence"]
    counts = data["counts"]
    mw = data["molecular_weight"]
    pi_data = data.get("isoelectric_point", {})

    lines = [
        "=" * 60,
        "GENEWEAVER DNA SEQUENCE ANALYSIS REPORT",
        "=" * 60,
        f"Sequence Length : {data['length']} bp",
        f"GC Content      : {data['gc_content']}%",
        f"AT Content      : {data['at_content']}%",
        f"GC Skew         : {data['gc_skew']}",
        f"Melting Temp Tm : {data['melting_temperature']} °C",
        f"Mol Weight (ss) : {mw['ssDNA']} g/mol",
        f"Mol Weight (ds) : {mw['dsDNA']} g/mol",
        f"Isoelectric (pI): {pi_data.get('pi', 'N/A')}",
        f"Charge (pH 7.4) : {pi_data.get('charge_ph74', 'N/A')}",
        "-" * 60,
        "NUCLEOTIDE COMPOSITION",
        f"  Adenine (A)   : {counts['A']} ({round(counts['A']/data['length']*100, 2)}%)",
        f"  Thymine (T)   : {counts['T']} ({round(counts['T']/data['length']*100, 2)}%)",
        f"  Guanine (G)   : {counts['G']} ({round(counts['G']/data['length']*100, 2)}%)",
        f"  Cytosine (C)  : {counts['C']} ({round(counts['C']/data['length']*100, 2)}%)",
        "-" * 60,
        "TRANSFORMATION PRODUCTS",
        "Original DNA (5' -> 3'):",
        seq,
        "",
        "Reverse Complement (5' -> 3'):",
        data['reverse_complement'],
        "",
        "mRNA Transcript (5' -> 3'):",
        data['rna'],
        "",
        "Primary Protein Translation:",
        data['protein'],
        "-" * 60,
        f"OPEN READING FRAMES IDENTIFIED: {len(data['orfs'])}",
    ]

    for idx, orf in enumerate(data["orfs"], 1):
        lines.append(
            f"  [ORF {idx}] Frame +{orf['frame']} | Pos: {orf['start']}-{orf['end']} "
            f"({orf['length_bp']} bp, {orf['codons']} aa) -> {orf['protein']}"
        )

    lines.append("-" * 60)
    lines.append("RESTRICTION ENZYME DIGESTION SUMMARY:")
    cutters = [e for e in data.get("restriction_sites", []) if e["cuts_sequence"]]
    if cutters:
        for enz in cutters:
            cuts_str = ", ".join(str(p) for p in enz["cut_positions"])
            frags_str = ", ".join(f"{f}bp" for f in enz["fragments"])
            lines.append(f"  {enz['name']} ({enz['site']}): {enz['cuts']} cut(s) at [{cuts_str}] -> Fragments: {frags_str}")
    else:
        lines.append("  No restriction sites for standard enzymes detected in sequence.")

    lines.append("=" * 60)
    lines.append("Report generated by GeneWeaver · Educational Bioinformatics System")
    lines.append("=" * 60)

    return "\n".join(lines)


def analyze_sequence(raw_sequence: str) -> Dict[str, Any]:
    """Execute complete analysis pipeline on a raw sequence string."""
    sequence = clean_sequence(raw_sequence)
    length = len(sequence)
    counts = {base: sequence.count(base) for base in "ATGC"}
    gc_content = ((counts["G"] + counts["C"]) / length) * 100
    at_content = 100.0 - gc_content
    mw = calculate_molecular_weight(sequence)
    tm = calculate_melting_temperature(sequence)
    skew = calculate_gc_skew(sequence)
    orfs = find_orfs(sequence)
    prot = translate(sequence)
    pi_data = calculate_isoelectric_point(prot)
    restriction = find_restriction_sites(sequence)

    percentages = {
        base: round((count / length) * 100, 1)
        for base, count in counts.items()
    }

    cutting_enzymes = [e for e in restriction if e["cuts_sequence"]]

    return {
        "sequence": sequence,
        "length": length,
        "counts": counts,
        "percentages": percentages,
        "gc_content": round(gc_content, 2),
        "at_content": round(at_content, 2),
        "gc_skew": skew,
        "melting_temperature": tm,
        "molecular_weight": mw,
        "reverse_complement": reverse_complement(sequence),
        "rna": transcribe(sequence),
        "protein": prot,
        "codons": length // 3,
        "orfs": orfs,
        "has_orfs": len(orfs) > 0,
        "isoelectric_point": pi_data,
        "restriction_sites": restriction,
        "cutting_enzymes": cutting_enzymes,
        "has_cutters": len(cutting_enzymes) > 0,
    }
