import tempfile
import unittest
from pathlib import Path

from motifminer.pipeline import parse_blast_table


class BlastParsingTests(unittest.TestCase):
    def test_filters_and_keeps_best_hsp(self):
        rows = (
            "A\t50\t80\t1\t80\t1e-20\t100\t100\t110\n"
            "A\t55\t90\t1\t90\t1e-30\t150\t100\t110\n"
            "B\t15\t90\t1\t90\t1e-10\t120\t100\t100\n"
            "C\t40\t40\t1\t40\t1e-10\t90\t100\t100\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "hits.tsv"
            path.write_text(rows, encoding="utf-8")
            hits = parse_blast_table(
                path,
                minimum_query_coverage=0.70,
                minimum_identity=0.20,
                maximum_identity=1.00,
            )
        self.assertEqual([hit.accession for hit in hits], ["A"])
        self.assertEqual(hits[0].bitscore, 150.0)
        self.assertEqual(hits[0].query_coverage, 0.90)


if __name__ == "__main__":
    unittest.main()
