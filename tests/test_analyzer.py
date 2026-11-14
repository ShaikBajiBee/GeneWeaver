"""Unit tests for GeneWeaver bioinformatics analysis engine."""

import unittest
from analyzer import (
    clean_sequence,
    parse_fasta,
    reverse_complement,
    transcribe,
    translate,
    calculate_molecular_weight,
    calculate_melting_temperature,
    calculate_gc_skew,
    find_orfs,
    analyze_sequence,
    SequenceError,
)


class TestAnalyzer(unittest.TestCase):

    def setUp(self):
        # Known test sequence: Start codon ATG -> GCC (Ala) -> ATT (Ile) -> TGA (Stop)
        self.sample_dna = "ATGGCCATTTGA"

    def test_clean_sequence_valid(self):
        cleaned = clean_sequence("  atg  gcc\natt tga ")
        self.assertEqual(cleaned, "ATGGCCATTTGA")

    def test_clean_sequence_empty(self):
        with self.assertRaises(SequenceError):
            clean_sequence("")
        with self.assertRaises(SequenceError):
            clean_sequence("   \n\t  ")

    def test_clean_sequence_invalid_characters(self):
        with self.assertRaises(SequenceError) as ctx:
            clean_sequence("ATGCZN123")
        self.assertIn("Invalid DNA character(s)", str(ctx.exception))

    def test_reverse_complement(self):
        # 5'-ATGC-3' -> complement 5'-TACG-3' -> reverse 5'-GCAT-3'
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("AAAA"), "TTTT")

    def test_transcribe(self):
        self.assertEqual(transcribe("ATGGCCATTTGA"), "AUGGCCAUUUGA")

    def test_translate(self):
        # ATG -> M, GCC -> A, ATT -> I, TGA -> Stop
        result = translate("ATGGCCATTTGA")
        self.assertEqual(result, "M A I Stop")

    def test_translate_no_stop(self):
        result = translate("ATGGCCATT")
        self.assertEqual(result, "M A I")

    def test_molecular_weight(self):
        mw = calculate_molecular_weight("ATGC")
        self.assertIsInstance(mw, dict)
        self.assertIn("ssDNA", mw)
        self.assertIn("dsDNA", mw)
        self.assertGreater(mw["dsDNA"], mw["ssDNA"])

    def test_melting_temperature_short(self):
        # Wallace rule: 2*(A+T) + 4*(G+C)
        # ATGC: 2*2 + 4*2 = 12
        tm = calculate_melting_temperature("ATGC")
        self.assertEqual(tm, 12.0)

    def test_melting_temperature_long(self):
        long_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        tm = calculate_melting_temperature(long_seq)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 50.0)

    def test_gc_skew(self):
        # Equal G and C: skew should be 0.0
        self.assertEqual(calculate_gc_skew("GGCC"), 0.0)
        # All G: skew should be 1.0
        self.assertEqual(calculate_gc_skew("GGGG"), 1.0)
        # All C: skew should be -1.0
        self.assertEqual(calculate_gc_skew("CCCC"), -1.0)

    def test_find_orfs(self):
        # ATG (M) GCC (A) ATT (I) GTA (V) TAA (Stop) -> 4 amino acids
        dna = "ATGGCCATTGTATAA"
        orfs = find_orfs(dna, min_codon_length=2)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(orfs[0]["start"], 1)
        self.assertEqual(orfs[0]["end"], 15)
        self.assertEqual(orfs[0]["protein"], "MAIV")

    def test_parse_fasta_single(self):
        fasta = ">TestSeq\nATGG\nCCAT\nTTGA\n"
        records = parse_fasta(fasta)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["header"], "TestSeq")
        self.assertEqual(records[0]["sequence"], "ATGGCCATTTGA")

    def test_parse_fasta_multi(self):
        fasta = ">Seq1\nATGC\n>Seq2\nCGTA\n"
        records = parse_fasta(fasta)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["header"], "Seq1")
        self.assertEqual(records[0]["sequence"], "ATGC")
        self.assertEqual(records[1]["header"], "Seq2")
        self.assertEqual(records[1]["sequence"], "CGTA")

    def test_analyze_sequence_full_pipeline(self):
        report = analyze_sequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        self.assertEqual(report["length"], 39)
        self.assertIn("counts", report)
        self.assertIn("gc_content", report)
        self.assertIn("percentages", report)
        self.assertIn("melting_temperature", report)
        self.assertIn("molecular_weight", report)
        self.assertTrue(report["has_orfs"])


if __name__ == "__main__":
    unittest.main()
