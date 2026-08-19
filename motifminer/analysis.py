from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass

from motifminer.exchange import infer_exchangeable
from motifminer.io import BindingSite, SequenceRecord, STANDARD_AA


@dataclass(frozen=True)
class SiteResult:
    position: int
    reference_residue: str
    alignment_column: int
    total_n: int
    nongap_n: int
    gap_n: int
    unknown_n: int
    gap_fraction: float
    frequencies: dict[str, float]
    top_residues: list[dict[str, float | str]]
    exact_conservation: float
    automatic_residues: list[str]
    automatic_grouped_conservation: float
    extended_residues: list[str]
    extended_grouped_conservation: float
    user_residues: list[str]
    user_grouped_conservation: float
    alignment_quality: str
    decision: str

    def to_dict(self) -> dict:
        return asdict(self)


def validate_alignment(records: list[SequenceRecord]) -> int:
    lengths = {len(record.sequence) for record in records}
    if len(lengths) != 1:
        raise ValueError("All MSA sequences must have the same aligned length")
    return lengths.pop()


def find_aligned_reference(
    reference: SequenceRecord, records: list[SequenceRecord]
) -> SequenceRecord:
    id_matches = [record for record in records if record.identifier == reference.identifier]
    if id_matches:
        candidate = id_matches[0]
        if candidate.sequence.replace("-", "") != reference.sequence:
            raise ValueError(
                "MSA record matching the reference identifier does not match the reference sequence"
            )
        return candidate
    sequence_matches = [
        record for record in records if record.sequence.replace("-", "") == reference.sequence
    ]
    if len(sequence_matches) == 1:
        return sequence_matches[0]
    if not sequence_matches:
        raise ValueError("The reference protein could not be found in the MSA")
    raise ValueError("Multiple MSA records match the reference; use the same unique identifier")


def reference_position_map(aligned_reference: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    reference_position = 0
    for column_index, residue in enumerate(aligned_reference, start=1):
        if residue != "-":
            reference_position += 1
            mapping[reference_position] = column_index
    return mapping


def _grouped_frequency(residues: set[str] | frozenset[str], frequencies: dict[str, float]) -> float:
    return sum(frequencies.get(residue, 0.0) for residue in residues)


def analyze_site(
    site: BindingSite,
    alignment_column: int,
    records: list[SequenceRecord],
    *,
    top_n: int = 4,
    minimum_exchange_frequency: float = 0.05,
    conservation_threshold: float = 0.80,
    maximum_gap_fraction: float = 0.20,
) -> SiteResult:
    column = [record.sequence[alignment_column - 1] for record in records]
    counts = Counter(residue for residue in column if residue in STANDARD_AA)
    gap_n = column.count("-")
    unknown_n = len(column) - gap_n - sum(counts.values())
    nongap_n = sum(counts.values())
    if nongap_n == 0:
        raise ValueError(f"No amino acids are present at reference position {site.position}")
    frequencies = {
        residue: count / nongap_n
        for residue, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    }
    top_residues = [
        {"residue": residue, "frequency": frequency}
        for residue, frequency in list(frequencies.items())[:top_n]
    ]
    automatic = infer_exchangeable(
        site.residue,
        frequencies,
        minimum_frequency=minimum_exchange_frequency,
    )
    extended = infer_exchangeable(
        site.residue,
        frequencies,
        minimum_frequency=minimum_exchange_frequency,
        extended=True,
    )
    gap_fraction = gap_n / len(records)
    quality = "pass" if gap_fraction <= maximum_gap_fraction else "warning_high_gap"
    automatic_conservation = _grouped_frequency(automatic, frequencies)
    decision = (
        "candidate_motif_position"
        if quality == "pass" and automatic_conservation >= conservation_threshold
        else "review"
    )
    return SiteResult(
        position=site.position,
        reference_residue=site.residue,
        alignment_column=alignment_column,
        total_n=len(records),
        nongap_n=nongap_n,
        gap_n=gap_n,
        unknown_n=unknown_n,
        gap_fraction=gap_fraction,
        frequencies=frequencies,
        top_residues=top_residues,
        exact_conservation=frequencies.get(site.residue, 0.0),
        automatic_residues=sorted(automatic),
        automatic_grouped_conservation=automatic_conservation,
        extended_residues=sorted(extended),
        extended_grouped_conservation=_grouped_frequency(extended, frequencies),
        user_residues=sorted(site.accepted_residues),
        user_grouped_conservation=_grouped_frequency(site.accepted_residues, frequencies),
        alignment_quality=quality,
        decision=decision,
    )


def analyze_alignment(
    reference: SequenceRecord,
    sites: list[BindingSite],
    records: list[SequenceRecord],
    **kwargs,
) -> list[SiteResult]:
    validate_alignment(records)
    aligned_reference = find_aligned_reference(reference, records)
    mapping = reference_position_map(aligned_reference.sequence)
    return [
        analyze_site(site, mapping[site.position], records, **kwargs)
        for site in sites
    ]

