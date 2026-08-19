import tempfile
import unittest
from pathlib import Path

from motifminer.analysis import analyze_alignment, reference_position_map
from motifminer.io import BindingSite, SequenceRecord, read_binding_sites


class PositionMappingTests(unittest.TestCase):
    def test_reference_positions_map_around_gaps(self):
        self.assertEqual(reference_position_map("MA--LTY"), {1: 1, 2: 2, 3: 5, 4: 6, 5: 7})


class ConservationTests(unittest.TestCase):
    def test_nongap_frequencies_and_exchangeability(self):
        reference = SequenceRecord("ref", "ref", "MAYR")
        records = [
            SequenceRecord("ref", "ref", "MA-YR"),
            SequenceRecord("s2", "s2", "MA-FR"),
            SequenceRecord("s3", "s3", "MA-YR"),
            SequenceRecord("s4", "s4", "MA-WR"),
            SequenceRecord("s5", "s5", "MA--R"),
        ]
        site = BindingSite(3, "Y", frozenset("YF"))
        result = analyze_alignment(reference, [site], records)[0]
        self.assertEqual(result.alignment_column, 4)
        self.assertEqual(result.total_n, 5)
        self.assertEqual(result.nongap_n, 4)
        self.assertEqual(result.gap_n, 1)
        self.assertAlmostEqual(result.exact_conservation, 0.5)
        self.assertEqual(result.automatic_residues, ["F", "Y"])
        self.assertEqual(result.extended_residues, ["F", "W", "Y"])
        self.assertAlmostEqual(result.user_grouped_conservation, 0.75)

    def test_binding_site_mismatch_is_rejected(self):
        reference = SequenceRecord("ref", "ref", "MAYR")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sites.csv"
            path.write_text("position,residue\n3,F\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "supplied F, reference has Y"):
                read_binding_sites(path, reference)


if __name__ == "__main__":
    unittest.main()

