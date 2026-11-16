"""Unit tests for GeneWeaver reverse-translation and primer design module."""

import unittest
from primers import (
    reverse_translate,
    analyze_primer,
    design_pcr_primers,
)


class TestPrimers(unittest.TestCase):

    def test_reverse_translate_basic(self):
        # Methionine (M) has only 1 codon (ATG)
        res = reverse_translate("M")
        self.assertEqual(res["consensus_dna"], "ATG")
        self.assertEqual(res["total_combinations"], 1)

    def test_reverse_translate_multi_aa(self):
        # M (1 codon: ATG) + W (1 codon: TGG) -> 1 combination
        res = reverse_translate("MW")
        self.assertEqual(res["consensus_dna"], "ATGTGG")
        self.assertEqual(res["total_combinations"], 1)

        # M (1) + F (2: TTT, TTC) -> 2 combinations, consensus ATGTTC/TTT -> ATGTTY
        res2 = reverse_translate("MF")
        self.assertEqual(res2["consensus_dna"], "ATGTTY")
        self.assertEqual(res2["total_combinations"], 2)

    def test_reverse_translate_invalid_character(self):
        with self.assertRaises(ValueError):
            reverse_translate("MXZ")

    def test_analyze_primer_optimal(self):
        # 20 bp primer, 50% GC, 3' C clamp
        primer = "ATGCGATCGATCGATCGATC"
        info = analyze_primer(primer)
        self.assertEqual(info["length"], 20)
        self.assertTrue(info["has_gc_clamp"])
        self.assertGreater(info["melting_temperature"], 50.0)

    def test_design_pcr_primers(self):
        # 60 bp DNA template
        template = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGC"
        res = design_pcr_primers(template, primer_len=18)
        self.assertEqual(res["amplicon_size"], 60)
        self.assertIn("forward_primer", res)
        self.assertIn("reverse_primer", res)
        self.assertIsInstance(res["tm_difference"], float)


if __name__ == "__main__":
    unittest.main()
