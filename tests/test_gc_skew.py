"""Unit tests for sliding-window GC skew analysis and origin prediction."""

import unittest
from gc_skew_analyzer import (
    calculate_sliding_gc_skew,
    generate_skew_svg,
)


class TestGCSkew(unittest.TestCase):

    def setUp(self):
        # 40 bp test sequence
        self.seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    def test_calculate_sliding_gc_skew_basics(self):
        res = calculate_sliding_gc_skew(self.seq, window_size=10, step_size=5)
        self.assertEqual(res["sequence_length"], len(self.seq))
        self.assertGreater(res["total_windows"], 0)
        self.assertEqual(len(res["positions"]), len(res["window_skews"]))
        self.assertEqual(len(res["window_skews"]), len(res["cumulative_skews"]))
        self.assertIn("predicted_origin_bp", res)
        self.assertIn("predicted_terminus_bp", res)

    def test_calculate_sliding_gc_skew_short(self):
        res = calculate_sliding_gc_skew("ATGC", window_size=10)
        self.assertEqual(res["sequence_length"], 4)
        self.assertEqual(len(res["window_skews"]), 1)

    def test_generate_skew_svg(self):
        res = calculate_sliding_gc_skew(self.seq, window_size=10, step_size=5)
        svg = generate_skew_svg(res)
        self.assertTrue(svg.startswith("<svg"))
        self.assertTrue(svg.endswith("</svg>"))
        self.assertIn("polyline", svg)
        self.assertIn("Cumulative GC Skew", svg)


if __name__ == "__main__":
    unittest.main()
