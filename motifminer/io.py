from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = frozenset("BJZ")
MSA_SYMBOLS = STANDARD_AA | AMBIGUOUS_AA | frozenset("-X?")


@dataclass(frozen=True)
class SequenceRecord:
    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class BindingSite:
    position: int
    residue: str
    accepted_residues: frozenset[str]


def read_fasta(path: str | Path, *, aligned: bool = False) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    header: str | None = None
    chunks: list[str] = []

    def finish_record() -> None:
        if header is None:
            return
        sequence = "".join(chunks).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"FASTA record {header!r} has no sequence")
        allowed = MSA_SYMBOLS if aligned else STANDARD_AA | frozenset("X?")
        invalid = sorted(set(sequence) - allowed)
        if invalid:
            raise ValueError(
                f"FASTA record {header!r} contains invalid symbols: {', '.join(invalid)}"
            )
        identifier = header.split()[0]
        records.append(SequenceRecord(identifier, header, sequence))

    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError("FASTA contains an empty header")
                chunks = []
            else:
                if header is None:
                    raise ValueError("FASTA sequence appears before its first header")
                chunks.append(line)
    finish_record()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    identifiers = [record.identifier for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("FASTA identifiers must be unique")
    return records


def read_reference(path: str | Path) -> SequenceRecord:
    records = read_fasta(path)
    if len(records) != 1:
        raise ValueError("Reference FASTA must contain exactly one protein sequence")
    return records[0]


def _parse_accepted(value: str, reference_residue: str) -> frozenset[str]:
    cleaned = value.upper().replace(",", "").replace("/", "").replace(" ", "")
    residues = set(cleaned) if cleaned else {reference_residue}
    residues.add(reference_residue)
    invalid = residues - STANDARD_AA
    if invalid:
        raise ValueError(f"Invalid accepted amino acids: {', '.join(sorted(invalid))}")
    return frozenset(residues)


def read_binding_sites(path: str | Path, reference: SequenceRecord) -> list[BindingSite]:
    sites: list[BindingSite] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"position", "residue"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Binding-site CSV requires position and residue columns")
        for row_number, row in enumerate(reader, start=2):
            try:
                position = int(row["position"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid position on CSV row {row_number}") from exc
            residue = (row["residue"] or "").strip().upper()
            if position < 1 or position > len(reference.sequence):
                raise ValueError(
                    f"Position {position} is outside the reference sequence (length {len(reference.sequence)})"
                )
            if residue not in STANDARD_AA:
                raise ValueError(f"Invalid residue {residue!r} at position {position}")
            observed = reference.sequence[position - 1]
            if observed != residue:
                raise ValueError(
                    f"Binding-site mismatch at position {position}: supplied {residue}, reference has {observed}"
                )
            accepted = _parse_accepted(row.get("accepted_residues", ""), residue)
            sites.append(BindingSite(position, residue, accepted))
    if not sites:
        raise ValueError("Binding-site CSV contains no sites")
    positions = [site.position for site in sites]
    if len(positions) != len(set(positions)):
        raise ValueError("Binding-site positions must be unique")
    return sorted(sites, key=lambda site: site.position)

