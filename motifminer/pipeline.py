from __future__ import annotations

import json
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict
from pathlib import Path

from motifminer.analysis import analyze_alignment
from motifminer.io import read_binding_sites, read_fasta, read_reference


@dataclass(frozen=True)
class BlastHit:
    accession: str
    identity: float
    query_coverage: float
    evalue: float
    bitscore: float
    subject_length: int


def require_program(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(
            f"Required program {name!r} was not found. Install it or provide its full path."
        )
    return path


def parse_blast_table(
    path: str | Path,
    *,
    minimum_query_coverage: float,
    minimum_identity: float,
    maximum_identity: float,
) -> list[BlastHit]:
    best: dict[str, BlastHit] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Unexpected BLAST table format on line {line_number}")
            accession, pident, length, qstart, qend, evalue, bitscore, qlen, slen = fields
            aligned_query_length = abs(int(qend) - int(qstart)) + 1
            hit = BlastHit(
                accession=accession,
                identity=float(pident) / 100.0,
                query_coverage=aligned_query_length / int(qlen),
                evalue=float(evalue),
                bitscore=float(bitscore),
                subject_length=int(slen),
            )
            if not (
                hit.query_coverage >= minimum_query_coverage
                and minimum_identity <= hit.identity <= maximum_identity
            ):
                continue
            previous = best.get(accession)
            if previous is None or hit.bitscore > previous.bitscore:
                best[accession] = hit
    return sorted(best.values(), key=lambda hit: (-hit.bitscore, hit.accession))


def run_remote_blast(
    reference_path: str | Path,
    output_path: str | Path,
    *,
    blastp: str,
    database: str,
    maximum_hits: int,
    evalue: float,
) -> None:
    outfmt = "6 sacc pident length qstart qend evalue bitscore qlen slen"
    command = [
        blastp,
        "-remote",
        "-query",
        str(reference_path),
        "-db",
        database,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(maximum_hits),
        "-outfmt",
        outfmt,
        "-out",
        str(output_path),
    ]
    subprocess.run(command, check=True)


def fetch_protein_fasta(
    accessions: list[str],
    output_path: str | Path,
    *,
    email: str,
    batch_size: int = 200,
) -> None:
    output = Path(output_path)
    with output.open("w", encoding="utf-8") as handle:
        for start in range(0, len(accessions), batch_size):
            batch = accessions[start : start + batch_size]
            payload = urllib.parse.urlencode(
                {
                    "db": "protein",
                    "id": ",".join(batch),
                    "rettype": "fasta",
                    "retmode": "text",
                    "email": email,
                    "tool": "motifminer",
                }
            ).encode()
            request = urllib.request.Request(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                data=payload,
                headers={"User-Agent": f"MotifMiner/0.2 ({email})"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                handle.write(response.read().decode("utf-8"))
            if start + batch_size < len(accessions):
                time.sleep(0.4)


def make_alignment_input(reference_path: Path, homolog_path: Path, output_path: Path) -> int:
    reference = read_reference(reference_path)
    homologs = read_fasta(homolog_path)
    seen_sequences = {reference.sequence}
    retained = []
    for record in homologs:
        if record.sequence not in seen_sequences:
            retained.append(record)
            seen_sequences.add(record.sequence)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{reference.identifier}\n{reference.sequence}\n")
        for record in retained:
            handle.write(f">{record.identifier}\n{record.sequence}\n")
    return len(retained) + 1


def run_mafft(input_path: Path, output_path: Path, *, mafft: str, threads: int = 1) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [mafft, "--auto", "--thread", str(threads), str(input_path)],
            stdout=handle,
            check=True,
        )


def run_pipeline(
    *,
    reference_path: str | Path,
    sites_path: str | Path,
    output_directory: str | Path,
    email: str,
    blastp_program: str = "blastp",
    mafft_program: str = "mafft",
    database: str = "nr_clustered",
    maximum_hits: int = 5000,
    evalue: float = 1e-5,
    minimum_query_coverage: float = 0.70,
    minimum_identity: float = 0.20,
    maximum_identity: float = 1.00,
    threads: int = 1,
    force: bool = False,
):
    reference_path = Path(reference_path)
    sites_path = Path(sites_path)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    reference = read_reference(reference_path)
    sites = read_binding_sites(sites_path, reference)
    blastp = require_program(blastp_program)
    mafft = require_program(mafft_program)

    blast_table = output / "blast_hits.tsv"
    if force or not blast_table.exists():
        run_remote_blast(
            reference_path,
            blast_table,
            blastp=blastp,
            database=database,
            maximum_hits=maximum_hits,
            evalue=evalue,
        )
    hits = parse_blast_table(
        blast_table,
        minimum_query_coverage=minimum_query_coverage,
        minimum_identity=minimum_identity,
        maximum_identity=maximum_identity,
    )
    if not hits:
        raise RuntimeError("BLAST completed but no hits passed the configured filters")

    homologs = output / "homologs.fasta"
    if force or not homologs.exists():
        fetch_protein_fasta([hit.accession for hit in hits], homologs, email=email)
    alignment_input = output / "alignment_input.fasta"
    sequence_count = make_alignment_input(reference_path, homologs, alignment_input)
    msa_path = output / "alignment.fasta"
    if force or not msa_path.exists():
        run_mafft(alignment_input, msa_path, mafft=mafft, threads=threads)

    records = read_fasta(msa_path, aligned=True)
    results = analyze_alignment(reference, sites, records)
    metadata = {
        "database": database,
        "maximum_hits": maximum_hits,
        "evalue": evalue,
        "minimum_query_coverage": minimum_query_coverage,
        "minimum_identity": minimum_identity,
        "maximum_identity": maximum_identity,
        "threads": threads,
        "passing_blast_hits": len(hits),
        "alignment_sequences": sequence_count,
        "blastp": blastp,
        "mafft": mafft,
    }
    with (output / "pipeline.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")
    return results, metadata
