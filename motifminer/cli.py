from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from motifminer.analysis import analyze_alignment
from motifminer.io import read_binding_sites, read_fasta, read_reference
from motifminer.pipeline import run_pipeline


def _write_outputs(results, output_directory: Path, parameters: dict) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    report = {
        "parameters": parameters,
        "sites": [result.to_dict() for result in results],
    }
    with (output_directory / "conservation.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    fields = [
        "position",
        "reference_residue",
        "alignment_column",
        "total_n",
        "nongap_n",
        "gap_n",
        "unknown_n",
        "gap_fraction",
        "top_residues",
        "exact_conservation",
        "automatic_residues",
        "automatic_grouped_conservation",
        "extended_residues",
        "extended_grouped_conservation",
        "user_residues",
        "user_grouped_conservation",
        "alignment_quality",
        "decision",
    ]
    with (output_directory / "conservation.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in results:
            row = result.to_dict()
            row["top_residues"] = ";".join(
                f"{item['residue']}:{item['frequency']:.4f}" for item in result.top_residues
            )
            for field in ("automatic_residues", "extended_residues", "user_residues"):
                row[field] = ",".join(row[field])
            writer.writerow({field: row[field] for field in fields})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="motifminer")
    subparsers = parser.add_subparsers(dest="command", required=True)
    learn = subparsers.add_parser("learn", help="Analyze binding-site conservation in an MSA")
    learn.add_argument("--reference", required=True, help="Single-record reference FASTA")
    learn.add_argument("--sites", required=True, help="Binding-site CSV")
    learn.add_argument("--msa", required=True, help="Aligned protein FASTA containing the reference")
    learn.add_argument("--output", required=True, help="Output directory")
    learn.add_argument("--top-n", type=int, default=4)
    learn.add_argument("--min-exchange-frequency", type=float, default=0.05)
    learn.add_argument("--conservation-threshold", type=float, default=0.80)
    learn.add_argument("--max-gap-fraction", type=float, default=0.20)

    run = subparsers.add_parser(
        "run", help="Run remote ClusteredNR search, retrieve homologs, align, and analyze"
    )
    run.add_argument("--reference", required=True, help="Single-record reference FASTA")
    run.add_argument("--sites", required=True, help="Binding-site CSV")
    run.add_argument("--output", required=True, help="Output directory")
    run.add_argument("--email", required=True, help="Email sent to NCBI with E-utilities requests")
    run.add_argument("--database", default="nr")
    run.add_argument("--max-hits", type=int, default=5000)
    run.add_argument("--evalue", type=float, default=1e-5)
    run.add_argument("--min-query-coverage", type=float, default=0.70)
    run.add_argument("--min-identity", type=float, default=0.20)
    run.add_argument("--max-identity", type=float, default=1.00)
    run.add_argument("--blastp", default="blastp", help="blastp executable or full path")
    run.add_argument("--mafft", default="mafft", help="MAFFT executable or full path")
    run.add_argument("--threads", type=int, default=1, help="CPU threads used by MAFFT")
    run.add_argument("--force", action="store_true", help="Re-run completed stages")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "learn":
        reference = read_reference(args.reference)
        sites = read_binding_sites(args.sites, reference)
        records = read_fasta(args.msa, aligned=True)
        parameters = {
            "top_n": args.top_n,
            "minimum_exchange_frequency": args.min_exchange_frequency,
            "conservation_threshold": args.conservation_threshold,
            "maximum_gap_fraction": args.max_gap_fraction,
        }
        results = analyze_alignment(reference, sites, records, **parameters)
        output = Path(args.output)
        _write_outputs(results, output, parameters)
        print(f"Analyzed {len(results)} binding sites across {len(records)} aligned sequences")
        print(f"Results: {output.resolve()}")
    elif args.command == "run":
        results, metadata = run_pipeline(
            reference_path=args.reference,
            sites_path=args.sites,
            output_directory=args.output,
            email=args.email,
            blastp_program=args.blastp,
            mafft_program=args.mafft,
            database=args.database,
            maximum_hits=args.max_hits,
            evalue=args.evalue,
            minimum_query_coverage=args.min_query_coverage,
            minimum_identity=args.min_identity,
            maximum_identity=args.max_identity,
            threads=args.threads,
            force=args.force,
        )
        output = Path(args.output)
        parameters = {
            "top_n": 4,
            "minimum_exchange_frequency": 0.05,
            "conservation_threshold": 0.80,
            "maximum_gap_fraction": 0.20,
            **metadata,
        }
        _write_outputs(results, output, parameters)
        print(f"Analyzed {len(results)} binding sites across {metadata['alignment_sequences']} sequences")
        print(f"Results: {output.resolve()}")


if __name__ == "__main__":
    main()
