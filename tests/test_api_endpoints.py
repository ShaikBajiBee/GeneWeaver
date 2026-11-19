"""
Integration test suite for GeneWeaver REST API endpoints.
Verifies status codes, JSON response schemas, and error handling across all routes.
"""

import unittest
from app import app


class TestAPIEndpoints(unittest.TestCase):

    def setUp(self):
        self.client = app.test_client()
        self.client.testing = True

    def test_health_check(self):
        res = self.client.get("/api/health")
        self.assertEqual(res.status_code, 200)
        data = res.get_json()
        self.assertEqual(data["status"], "healthy")
        self.assertIn("version", data)

    def test_analyze_endpoint_valid(self):
        payload = {"sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"}
        res = self.client.post("/api/analyze", json=payload)
        self.assertEqual(res.status_code, 200)
        data = res.get_json()
        self.assertTrue(data["success"])
        report = data["data"]
        self.assertEqual(report["length"], 39)
        self.assertIn("gc_content", report)
        self.assertIn("molecular_weight", report)
        self.assertIn("orfs", report)

    def test_analyze_endpoint_invalid_sequence(self):
        payload = {"sequence": "ATGCZ123"}
        res = self.client.post("/api/analyze", json=payload)
        self.assertEqual(res.status_code, 422)
        data = res.get_json()
        self.assertFalse(data["success"])
        self.assertIn("error", data)

    def test_compare_endpoint(self):
        payload = {
            "reference": "ATGGCCATTGTA",
            "query": "ATGGCCATGGTA"
        }
        res = self.client.post("/api/compare", json=payload)
        self.assertEqual(res.status_code, 200)
        data = res.get_json()
        self.assertTrue(data["success"])
        self.assertEqual(data["data"]["mutations_count"], 1)

    def test_primers_endpoint(self):
        payload = {
            "sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGC",
            "primer_len": 18
        }
        res = self.client.post("/api/primers", json=payload)
        self.assertEqual(res.status_code, 200)
        data = res.get_json()
        self.assertTrue(data["success"])
        self.assertIn("forward_primer", data["data"])
        self.assertIn("reverse_primer", data["data"])

    def test_reverse_translate_endpoint(self):
        payload = {"peptide": "MAIV"}
        res = self.client.post("/api/reverse-translate", json=payload)
        self.assertEqual(res.status_code, 200)
        data = res.get_json()
        self.assertTrue(data["success"])
        self.assertIn("consensus_dna", data["data"])

    def test_batch_analyze_endpoint(self):
        fasta = ">Seq1\nATGC\n>Seq2\nCGTA\n"
        res = self.client.post("/api/batch-analyze", json={"fasta": fasta})
        self.assertEqual(res.status_code, 200)
        data = res.get_json()
        self.assertTrue(data["success"])
        self.assertEqual(data["data"]["total_sequences"], 2)


if __name__ == "__main__":
    unittest.main()
