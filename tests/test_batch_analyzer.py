"""Unit tests for GeneWeaver batch FASTA processing and distance matrix analysis."""

import unittest
from batch_analyzer import (
    compute_sequence_distance,
    analyze_batch_fasta,
)
from analyzer import SequenceError


class TestBatchAnalyzer(unittest.TestCase):

    def setUp(self):
        self.multi_fasta = """
>Seq1_Homo_sapiens
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>Seq2_Pan_troglodytes
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAA
>Seq3_Mus_musculus
ATGGCCATGGTAATGGGCCGCTGAAAGGGTACCCGATAG
"""

    def test_compute_sequence_distance_identical(self):
        res = compute_sequence_distance("ATGC", "ATGC")
        self.assertEqual(res["mismatches"], 0)
        self.assertEqual(res["identity_pct"], 100.0)

    def test_compute_sequence_distance_mismatch(self):
        res = compute_sequence_distance("ATGC", "ATGA")
        self.assertEqual(res["mismatches"], 1)
        self.assertEqual(res["identity_pct"], 75.0)

    def test_analyze_batch_fasta_metrics(self):
        batch = analyze_batch_fasta(self.multi_fasta)
        self.assertEqual(batch["total_sequences"], 3)
        self.assertEqual(len(batch["headers"]), 3)
        self.assertIn("Seq1_Homo_sapiens", batch["headers"][0])

        metrics = batch["metrics"]
        self.assertEqual(metrics["avg_length_bp"], 39.0)
        self.assertEqual(metrics["min_length_bp"], 39)
        self.assertEqual(metrics["max_length_bp"], 39)

    def test_analyze_batch_distance_matrix(self):
        batch = analyze_batch_fasta(self.multi_fasta)
        matrix = batch["distance_matrix"]
        self.assertEqual(len(matrix), 3)
        self.assertEqual(len(matrix[0]), 3)

        # Diagonal must be 100% identity
        for i in range(3):
            self.assertEqual(matrix[i][i], 100.0)

        # Closest pair check
        self.assertIsNotNone(batch["closest_pair"])
        self.assertGreater(batch["closest_pair"]["identity_pct"], 90.0)

    def test_analyze_batch_empty_fails(self):
        with self.assertRaises(SequenceError):
            analyze_batch_fasta("")


if __name__ == "__main__":
    unittest.main()
