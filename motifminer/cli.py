from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from motifminer.analysis import analyze_alignment
from motifminer.io import read_binding_sites, read_fasta, read_reference
from motifminer.pipeline import run_local_pipeline, run_pipeline


FIELDS = ["position", "reference_residue", "alignment_column", "total_n", "nongap_n", "gap_n",
          "unknown_n", "gap_fraction", "top_residues", "exact_conservation", "automatic_residues",
          "automatic_grouped_conservation", "extended_residues", "extended_grouped_conservation",
          "user_residues", "user_grouped_conservation", "alignment_quality", "decision"]


def _write_outputs(results, output_directory: Path, parameters: dict) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    with (output_directory / "conservation.json").open("w", encoding="utf-8") as handle:
        json.dump({"parameters": parameters, "sites": [r.to_dict() for r in results]}, handle, indent=2); handle.write("\n")
    with (output_directory / "conservation.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t"); writer.writeheader()
        for result in results:
            row = result.to_dict()
            row["top_residues"] = ";".join(f"{x['residue']}:{x['frequency']:.4f}" for x in result.top_residues)
            for field in ("automatic_residues", "extended_residues", "user_residues"):
                row[field] = ",".join(row[field])
            writer.writerow({field: row[field] for field in FIELDS})


def _analysis_arguments(parser):
    parser.add_argument("--top-n", type=int, default=4)
    parser.add_argument("--min-exchange-frequency", type=float, default=0.05)
    parser.add_argument("--conservation-threshold", type=float, default=0.80)
    parser.add_argument("--max-gap-fraction", type=float, default=0.20)


def _blast_filter_arguments(parser):
    parser.add_argument("--max-hits", type=int, default=5000)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--min-query-coverage", type=float, default=0.70, help="Fraction from 0 to 1")
    parser.add_argument("--min-identity", type=float, default=0.20, help="Fraction from 0 to 1")
    parser.add_argument("--max-identity", type=float, default=1.00, help="Fraction from 0 to 1")
    parser.add_argument("--blastp", default="blastp", help="blastp executable or full path")
    parser.add_argument("--mafft", default="mafft", help="MAFFT executable or full path")
    parser.add_argument("--threads", type=int, default=1, help="CPU threads for BLAST and MAFFT")
    parser.add_argument("--force", action="store_true", help="Re-run completed expensive stages")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="motifminer")
    sub = parser.add_subparsers(dest="command", required=True)
    learn = sub.add_parser("learn", help="Analyze a user-provided protein MSA")
    learn.add_argument("--reference", required=True, help="Single-record reference FASTA")
    learn.add_argument("--sites", required=True, help="Binding-site CSV in reference coordinates")
    learn.add_argument("--msa", required=True, help="Aligned FASTA containing the reference")
    learn.add_argument("--output", required=True); _analysis_arguments(learn)

    run = sub.add_parser("run", help="Legacy remote NCBI BLAST/EFetch workflow")
    run.add_argument("--reference", required=True); run.add_argument("--sites", required=True)
    run.add_argument("--output", required=True); run.add_argument("--email", required=True)
    run.add_argument("--database", default="nr"); _blast_filter_arguments(run)

    local = sub.add_parser("local", help="Local BLAST workflow for any compatible protein database")
    local.add_argument("--reference", required=True, help="Full canonical reference protein FASTA")
    local.add_argument("--sites", required=True, help="Binding-site CSV using full-protein positions")
    local.add_argument("--output", required=True)
    local.add_argument("--query-mode", required=True, choices=("domain", "full"))
    local.add_argument("--domain", help="Domain/sensor FASTA; required in domain mode")
    local.add_argument("--database", required=True, help="Local makeblastdb protein database prefix/path")
    local.add_argument("--database-name", help="Optional display label (for metadata only)")
    _blast_filter_arguments(local)
    return parser


def _analysis_parameters(args):
    return {"top_n": args.top_n, "minimum_exchange_frequency": args.min_exchange_frequency,
            "conservation_threshold": args.conservation_threshold, "maximum_gap_fraction": args.max_gap_fraction}


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "learn":
        reference = read_reference(args.reference); sites = read_binding_sites(args.sites, reference)
        records = read_fasta(args.msa, aligned=True); parameters = _analysis_parameters(args)
        results = analyze_alignment(reference, sites, records, **parameters)
        metadata = {**parameters, "alignment_sequences": len(records)}
    elif args.command == "run":
        results, metadata = run_pipeline(reference_path=args.reference, sites_path=args.sites,
            output_directory=args.output, email=args.email, blastp_program=args.blastp,
            mafft_program=args.mafft, database=args.database, maximum_hits=args.max_hits,
            evalue=args.evalue, minimum_query_coverage=args.min_query_coverage,
            minimum_identity=args.min_identity, maximum_identity=args.max_identity,
            threads=args.threads, force=args.force)
    else:
        results, metadata = run_local_pipeline(reference_path=args.reference, sites_path=args.sites,
            output_directory=args.output, query_mode=args.query_mode, domain_path=args.domain,
            database=args.database, database_name=args.database_name, blastp_program=args.blastp,
            mafft_program=args.mafft, maximum_hits=args.max_hits, evalue=args.evalue,
            minimum_query_coverage=args.min_query_coverage, minimum_identity=args.min_identity,
            maximum_identity=args.max_identity, threads=args.threads, force=args.force)
    output = Path(args.output); _write_outputs(results, output, metadata)
    print(f"Analyzed {len(results)} binding sites across {metadata['alignment_sequences']} aligned sequences")
    print(f"Results: {output.resolve()}")


if __name__ == "__main__":
    main()
