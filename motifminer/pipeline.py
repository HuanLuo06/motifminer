from __future__ import annotations

import csv
import json
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path

from motifminer.analysis import analyze_alignment, find_aligned_reference, reference_position_map
from motifminer.io import (BindingSite, SequenceRecord, locate_unique_subsequence,
                           map_sites_to_domain, read_binding_sites, read_fasta, read_reference)


@dataclass(frozen=True)
class BlastHit:
    accession: str
    identity: float
    query_coverage: float
    evalue: float
    bitscore: float
    subject_length: int
    query_start: int = 0
    query_end: int = 0
    subject_start: int = 0
    subject_end: int = 0
    subject_sequence: str = ""


LOCAL_BLAST_FIELDS = ("sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
                      "evalue", "bitscore", "qlen", "slen", "sseq")


def require_program(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required program {name!r} was not found. Install it or provide its full path.")
    return path


def parse_blast_table(path, *, minimum_query_coverage, minimum_identity, maximum_identity,
                      maximum_evalue=float("inf"), local_format=False) -> list[BlastHit]:
    best: dict[str, BlastHit] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            expected = 12 if local_format else 9
            if len(fields) != expected:
                raise ValueError(f"Unexpected BLAST table format on line {line_number}: expected {expected} columns")
            if local_format:
                accession, pident, length, qs, qe, ss, se, ev, bits, qlen, slen, sseq = fields
            else:
                accession, pident, length, qs, qe, ev, bits, qlen, slen = fields
                ss = se = "0"; sseq = ""
            coverage = (abs(int(qe) - int(qs)) + 1) / int(qlen)
            hit = BlastHit(accession, float(pident) / 100.0, coverage, float(ev), float(bits),
                           int(slen), int(qs), int(qe), int(ss), int(se), sseq)
            if hit.evalue > maximum_evalue or hit.query_coverage < minimum_query_coverage or not minimum_identity <= hit.identity <= maximum_identity:
                continue
            previous = best.get(accession)
            if previous is None or hit.bitscore > previous.bitscore:
                best[accession] = hit
    return sorted(best.values(), key=lambda hit: (-hit.bitscore, hit.accession))


def run_remote_blast(reference_path, output_path, *, blastp, database, maximum_hits, evalue):
    subprocess.run([blastp, "-remote", "-query", str(reference_path), "-db", database,
                    "-evalue", str(evalue), "-max_target_seqs", str(maximum_hits),
                    "-outfmt", "6 sacc pident length qstart qend evalue bitscore qlen slen",
                    "-out", str(output_path)], check=True)


def run_local_blast(query_path, output_path, *, blastp, database, maximum_hits, evalue, threads):
    # sseq makes this independent of makeblastdb -parse_seqids and avoids blastdbcmd lookup.
    subprocess.run([blastp, "-query", str(query_path), "-db", str(database), "-evalue", str(evalue),
                    "-max_target_seqs", str(maximum_hits), "-num_threads", str(threads),
                    "-outfmt", "6 " + " ".join(LOCAL_BLAST_FIELDS), "-out", str(output_path)], check=True)


def fetch_protein_fasta(accessions, output_path, *, email, batch_size=200):
    with Path(output_path).open("w", encoding="utf-8") as handle:
        for start in range(0, len(accessions), batch_size):
            batch = accessions[start:start + batch_size]
            payload = urllib.parse.urlencode({"db": "protein", "id": ",".join(batch), "rettype": "fasta",
                                               "retmode": "text", "email": email, "tool": "motifminer"}).encode()
            request = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                             data=payload, headers={"User-Agent": f"MotifMiner/0.2 ({email})"})
            with urllib.request.urlopen(request, timeout=120) as response:
                handle.write(response.read().decode("utf-8"))
            if start + batch_size < len(accessions):
                time.sleep(0.4)


def _safe_identifier(identifier: str, number: int) -> str:
    cleaned = "".join(c if c.isalnum() or c in "_.-" else "_" for c in identifier)
    return f"hit{number:06d}_{cleaned[:120]}"


def write_filtered_hits(hits: list[BlastHit], path: Path) -> None:
    fields = [field.name for field in BlastHit.__dataclass_fields__.values()]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(asdict(hit) for hit in hits)


def write_homolog_regions(hits: list[BlastHit], path: Path) -> int:
    seen: set[str] = set()
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        for number, hit in enumerate(hits, 1):
            # BLAST sseq is the aligned subject region. Removing pairwise gaps yields the original
            # subject residues from sstart..send and works even without parse_seqids.
            sequence = hit.subject_sequence.replace("-", "").upper()
            if not sequence or sequence in seen:
                continue
            seen.add(sequence); count += 1
            handle.write(f">{_safe_identifier(hit.accession, number)} source={hit.accession} subject={hit.subject_start}-{hit.subject_end}\n{sequence}\n")
    return count


def make_alignment_input(reference_path: Path, homolog_path: Path, output_path: Path) -> int:
    reference = read_reference(reference_path); homologs = read_fasta(homolog_path)
    seen = {reference.sequence}; retained = []
    for record in homologs:
        if record.sequence not in seen:
            retained.append(record); seen.add(record.sequence)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(f">{reference.identifier}\n{reference.sequence}\n")
        for record in retained:
            handle.write(f">{record.identifier}\n{record.sequence}\n")
    return len(retained) + 1


def run_mafft(input_path, output_path, *, mafft, threads=1):
    with Path(output_path).open("w", encoding="utf-8") as handle:
        subprocess.run([mafft, "--auto", "--thread", str(threads), str(input_path)], stdout=handle, check=True)


def write_mapping_qc(path: Path, full_sites: list[BindingSite], query_sites: list[BindingSite],
                     query_reference: SequenceRecord, records: list[SequenceRecord]) -> None:
    columns = reference_position_map(find_aligned_reference(query_reference, records).sequence)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, ["full_position", "full_residue", "query_position", "query_residue", "msa_column"], delimiter="\t")
        writer.writeheader()
        for full, query in zip(full_sites, query_sites):
            writer.writerow({"full_position": full.position, "full_residue": full.residue,
                             "query_position": query.position, "query_residue": query.residue,
                             "msa_column": columns[query.position]})


def run_local_pipeline(*, reference_path, sites_path, output_directory, query_mode, database,
                       domain_path=None, database_name=None, blastp_program="blastp", mafft_program="mafft",
                       maximum_hits=5000, evalue=1e-5, minimum_query_coverage=0.70,
                       minimum_identity=0.20, maximum_identity=1.0, threads=1, force=False):
    if query_mode not in {"domain", "full"}:
        raise ValueError("query_mode must be 'domain' or 'full'")
    output = Path(output_directory); output.mkdir(parents=True, exist_ok=True)
    full_reference = read_reference(reference_path)
    full_sites = read_binding_sites(sites_path, full_reference)
    if query_mode == "domain":
        if not domain_path:
            raise ValueError("Domain mode requires --domain")
        query_reference = read_reference(domain_path)
        domain_start = locate_unique_subsequence(full_reference, query_reference)
        query_sites = map_sites_to_domain(full_sites, domain_start, len(query_reference.sequence))
        query_path = Path(domain_path)
    else:
        if domain_path:
            raise ValueError("--domain is only valid with --query-mode domain")
        query_reference = full_reference; query_sites = full_sites; query_path = Path(reference_path); domain_start = 1
    blastp = require_program(blastp_program); mafft = require_program(mafft_program)
    raw = output / "blast_hits_raw.tsv"
    if force or not raw.exists():
        run_local_blast(query_path, raw, blastp=blastp, database=database, maximum_hits=maximum_hits, evalue=evalue, threads=threads)
    hits = parse_blast_table(raw, minimum_query_coverage=minimum_query_coverage,
                             minimum_identity=minimum_identity, maximum_identity=maximum_identity,
                             maximum_evalue=evalue, local_format=True)
    if not hits:
        raise RuntimeError("BLAST completed but no hits passed the configured filters")
    write_filtered_hits(hits, output / "blast_hits_filtered.tsv")
    homologs = output / "homolog_regions.fasta"
    region_count = write_homolog_regions(hits, homologs)
    if not region_count:
        raise RuntimeError("Passing BLAST hits contained no extractable subject sequences")
    alignment_input = output / "alignment_input.fasta"
    sequence_count = make_alignment_input(query_path, homologs, alignment_input)
    msa = output / "alignment.fasta"
    if force or not msa.exists():
        run_mafft(alignment_input, msa, mafft=mafft, threads=threads)
    records = read_fasta(msa, aligned=True)
    results = analyze_alignment(query_reference, query_sites, records)
    write_mapping_qc(output / "site_mapping_qc.tsv", full_sites, query_sites, query_reference, records)
    metadata = {"workflow": "local", "query_mode": query_mode, "database": str(database),
                "database_name": database_name, "domain_start_in_full_reference": domain_start,
                "maximum_hits": maximum_hits, "evalue": evalue,
                "minimum_query_coverage": minimum_query_coverage, "minimum_identity": minimum_identity,
                "maximum_identity": maximum_identity, "threads": threads, "passing_blast_hits": len(hits),
                "homolog_regions": region_count, "alignment_sequences": sequence_count,
                "blastp": blastp, "mafft": mafft}
    with (output / "pipeline_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2); handle.write("\n")
    return results, metadata


def run_pipeline(*, reference_path, sites_path, output_directory, email, blastp_program="blastp",
                 mafft_program="mafft", database="nr_clustered", maximum_hits=5000, evalue=1e-5,
                 minimum_query_coverage=0.70, minimum_identity=0.20, maximum_identity=1.0,
                 threads=1, force=False):
    """Legacy remote workflow, deliberately retained for backward compatibility."""
    output = Path(output_directory); output.mkdir(parents=True, exist_ok=True)
    reference = read_reference(reference_path); sites = read_binding_sites(sites_path, reference)
    blastp = require_program(blastp_program); mafft = require_program(mafft_program)
    table = output / "blast_hits.tsv"
    if force or not table.exists():
        run_remote_blast(reference_path, table, blastp=blastp, database=database, maximum_hits=maximum_hits, evalue=evalue)
    hits = parse_blast_table(table, minimum_query_coverage=minimum_query_coverage,
                             minimum_identity=minimum_identity, maximum_identity=maximum_identity)
    if not hits: raise RuntimeError("BLAST completed but no hits passed the configured filters")
    homologs = output / "homologs.fasta"
    if force or not homologs.exists(): fetch_protein_fasta([h.accession for h in hits], homologs, email=email)
    alignment_input = output / "alignment_input.fasta"; sequence_count = make_alignment_input(Path(reference_path), homologs, alignment_input)
    msa = output / "alignment.fasta"
    if force or not msa.exists(): run_mafft(alignment_input, msa, mafft=mafft, threads=threads)
    results = analyze_alignment(reference, sites, read_fasta(msa, aligned=True))
    metadata = {"database": database, "maximum_hits": maximum_hits, "evalue": evalue,
                "minimum_query_coverage": minimum_query_coverage, "minimum_identity": minimum_identity,
                "maximum_identity": maximum_identity, "threads": threads, "passing_blast_hits": len(hits),
                "alignment_sequences": sequence_count, "blastp": blastp, "mafft": mafft}
    with (output / "pipeline.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2); handle.write("\n")
    return results, metadata
