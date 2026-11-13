"""DNA sequence validation and analysis utilities."""

class SequenceError(ValueError):
    """Raised when a DNA sequence is empty or contains invalid characters."""


CODON_TABLE = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}

VALID_BASES = set("ATGC")
COMPLEMENT = str.maketrans("ATGC", "TACG")


def clean_sequence(raw_sequence):
    """Remove whitespace and FASTA-style line breaks; normalize to uppercase."""
    sequence = "".join(raw_sequence.split()).upper()
    if not sequence:
        raise SequenceError("Please enter a DNA sequence or upload a FASTA file.")

    invalid = sorted(set(sequence) - VALID_BASES)
    if invalid:
        shown = ", ".join(invalid)
        raise SequenceError(
            f"Invalid DNA character(s): {shown}. Only A, T, G and C are accepted."
        )
    return sequence


def parse_fasta(contents):
    """Read the first FASTA record from text content."""
    lines = [line.strip() for line in contents.splitlines() if line.strip()]
    if not lines:
        raise SequenceError("The FASTA file is empty.")
    if not any(line.startswith(">") for line in lines):
        # Also accept a plain text sequence file.
        return "".join(lines)

    sequence_lines = []
    found_header = False
    for line in lines:
        if line.startswith(">"):
            if found_header:
                break
            found_header = True
        elif found_header:
            sequence_lines.append(line)
    if not sequence_lines:
        raise SequenceError("No sequence data was found after the FASTA header.")
    return "".join(sequence_lines)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def transcribe(sequence):
    return sequence.replace("T", "U")


def translate(sequence):
    """Translate complete codons from the first base; ignore a partial final codon."""
    protein = []
    for index in range(0, len(sequence) - 2, 3):
        amino_acid = CODON_TABLE[sequence[index:index + 3]]
        if amino_acid == "*":
            protein.append("Stop")
            break
        protein.append(amino_acid)
    return " ".join(protein) if protein else "No complete codon found"


def analyze_sequence(raw_sequence):
    sequence = clean_sequence(raw_sequence)
    length = len(sequence)
    counts = {base: sequence.count(base) for base in "ATGC"}
    gc_content = ((counts["G"] + counts["C"]) / length) * 100

    return {
        "sequence": sequence,
        "length": length,
        "counts": counts,
        "gc_content": round(gc_content, 2),
        "at_content": round(100 - gc_content, 2),
        "reverse_complement": reverse_complement(sequence),
        "rna": transcribe(sequence),
        "protein": translate(sequence),
        "codons": (length // 3),
    }
