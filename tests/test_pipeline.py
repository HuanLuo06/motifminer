import tempfile
import unittest
from pathlib import Path

from unittest.mock import patch

from motifminer.io import BindingSite, SequenceRecord, locate_unique_subsequence, map_sites_to_domain
from motifminer.pipeline import parse_blast_table, run_local_blast, write_homolog_regions


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

    def test_local_format_keeps_coordinates_and_subject_sequence(self):
        rows = "A\t50\t80\t1\t80\t11\t90\t1e-20\t100\t100\t120\tACD-EF\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "hits.tsv"
            path.write_text(rows, encoding="utf-8")
            hits = parse_blast_table(
                path,
                minimum_query_coverage=0.70,
                minimum_identity=0.20,
                maximum_identity=1.00,
                maximum_evalue=1e-5,
                local_format=True,
            )
            output = Path(directory) / "regions.fasta"
            self.assertEqual(write_homolog_regions(hits, output), 1)
            self.assertIn("ACDEF", output.read_text(encoding="utf-8"))
        self.assertEqual((hits[0].subject_start, hits[0].subject_end), (11, 90))

    @patch("motifminer.pipeline.subprocess.run")
    def test_local_blast_is_not_remote_and_uses_threads(self, run):
        run_local_blast(
            "query.fasta", "hits.tsv", blastp="blastp", database="/db/proteins",
            maximum_hits=100, evalue=1e-5, threads=8,
        )
        command = run.call_args.args[0]
        self.assertNotIn("-remote", command)
        self.assertEqual(command[command.index("-num_threads") + 1], "8")
        self.assertIn("sseq", command[command.index("-outfmt") + 1])


class DomainMappingTests(unittest.TestCase):
    def test_unique_domain_and_site_mapping(self):
        full = SequenceRecord("full", "full", "AACDEFGHII")
        domain = SequenceRecord("domain", "domain", "CDEFG")
        start = locate_unique_subsequence(full, domain)
        sites = [BindingSite(4, "D", frozenset("D"))]
        self.assertEqual(start, 3)
        self.assertEqual(map_sites_to_domain(sites, start, len(domain.sequence))[0].position, 2)

    def test_repeated_domain_is_rejected(self):
        full = SequenceRecord("full", "full", "ACDACD")
        domain = SequenceRecord("domain", "domain", "ACD")
        with self.assertRaisesRegex(ValueError, "occurs 2 times"):
            locate_unique_subsequence(full, domain)

    def test_site_outside_domain_is_rejected(self):
        sites = [BindingSite(2, "A", frozenset("A"))]
        with self.assertRaisesRegex(ValueError, "outside the domain"):
            map_sites_to_domain(sites, 3, 5)


if __name__ == "__main__":
    unittest.main()
