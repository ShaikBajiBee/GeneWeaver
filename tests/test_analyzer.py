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
    find_restriction_sites,
    calculate_isoelectric_point,
    compare_sequences,
    analyze_sequence,
    SequenceError,
)


class TestAnalyzer(unittest.TestCase):

    def setUp(self):
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
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("AAAA"), "TTTT")

    def test_transcribe(self):
        self.assertEqual(transcribe("ATGGCCATTTGA"), "AUGGCCAUUUGA")

    def test_translate(self):
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
        tm = calculate_melting_temperature("ATGC")
        self.assertEqual(tm, 12.0)

    def test_melting_temperature_long(self):
        long_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        tm = calculate_melting_temperature(long_seq)
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 50.0)

    def test_gc_skew(self):
        self.assertEqual(calculate_gc_skew("GGCC"), 0.0)
        self.assertEqual(calculate_gc_skew("GGGG"), 1.0)
        self.assertEqual(calculate_gc_skew("CCCC"), -1.0)

    def test_find_orfs(self):
        dna = "ATGGCCATTGTATAA"
        orfs = find_orfs(dna, min_codon_length=2)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(orfs[0]["start"], 1)
        self.assertEqual(orfs[0]["end"], 15)
        self.assertEqual(orfs[0]["protein"], "MAIV")

    def test_restriction_sites_ecori(self):
        dna = "AAAAAGAATTCGGGGGGGGG"
        sites = find_restriction_sites(dna)
        ecori = next((e for e in sites if e["name"] == "EcoRI"), None)
        self.assertIsNotNone(ecori)
        self.assertEqual(ecori["cuts"], 1)
        self.assertEqual(ecori["cut_positions"], [6])
        self.assertEqual(ecori["fragments"], [14, 6])

    def test_restriction_sites_no_cutters(self):
        dna = "AAAAAAAAAAAAAAAA"
        sites = find_restriction_sites(dna)
        for e in sites:
            self.assertEqual(e["cuts"], 0)
            self.assertEqual(e["fragments"], [16])

    def test_isoelectric_point(self):
        prot = "K K K K K"
        pi_data = calculate_isoelectric_point(prot)
        self.assertGreater(pi_data["pi"], 9.0)

        prot_acidic = "D D D D D"
        pi_acidic = calculate_isoelectric_point(prot_acidic)
        self.assertLess(pi_acidic["pi"], 4.5)

    def test_compare_sequences_identical(self):
        res = compare_sequences("ATGGCC", "ATGGCC")
        self.assertEqual(res["hamming_distance"], 0)
        self.assertEqual(res["identity_pct"], 100.0)
        self.assertEqual(len(res["mutations"]), 0)

    def test_compare_sequences_transition_transversion(self):
        # Position 1: A -> G (Transition: purine to purine)
        # Position 2: T -> G (Transversion: pyrimidine to purine)
        ref = "AT"
        alt = "GG"
        res = compare_sequences(ref, alt)
        self.assertEqual(res["transitions"], 1)
        self.assertEqual(res["transversions"], 1)
        self.assertEqual(res["mutations_count"], 2)

    def test_compare_sequences_synonymous_impact(self):
        # TTT (Phe) -> TTC (Phe) is synonymous
        res = compare_sequences("TTT", "TTC")
        self.assertEqual(res["synonymous_count"], 1)
        self.assertIn("Synonymous", res["mutations"][0]["impact"])

    def test_parse_fasta_single(self):
        fasta = ">TestSeq\nATGG\nCCAT\nTTGA\n"
        records = parse_fasta(fasta)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["header"], "TestSeq")
        self.assertEqual(records[0]["sequence"], "ATGGCCATTTGA")

    def test_analyze_sequence_full_pipeline(self):
        report = analyze_sequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        self.assertEqual(report["length"], 39)
        self.assertIn("counts", report)
        self.assertIn("gc_content", report)
        self.assertIn("percentages", report)
        self.assertIn("isoelectric_point", report)
        self.assertIn("restriction_sites", report)
        self.assertTrue(report["has_orfs"])


if __name__ == "__main__":
    unittest.main()
