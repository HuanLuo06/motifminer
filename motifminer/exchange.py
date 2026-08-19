from __future__ import annotations


# Conservative groups are intentionally narrow. Extended groups are displayed
# separately because similarity does not guarantee equivalent ligand binding.
HIGH_CONFIDENCE_GROUPS = (
    frozenset("FY"),
    frozenset("KR"),
    frozenset("DE"),
    frozenset("ILV"),
    frozenset("ST"),
    frozenset("NQ"),
)

EXTENDED_GROUPS = HIGH_CONFIDENCE_GROUPS + (
    frozenset("FYW"),
    frozenset("AG"),
)


def compatible_with(reference: str, candidate: str, *, extended: bool = False) -> bool:
    if reference == candidate:
        return True
    groups = EXTENDED_GROUPS if extended else HIGH_CONFIDENCE_GROUPS
    return any(reference in group and candidate in group for group in groups)


def infer_exchangeable(
    reference: str,
    frequencies: dict[str, float],
    *,
    minimum_frequency: float = 0.05,
    extended: bool = False,
) -> frozenset[str]:
    accepted = {reference}
    for residue, frequency in frequencies.items():
        if frequency >= minimum_frequency and compatible_with(
            reference, residue, extended=extended
        ):
            accepted.add(residue)
    return frozenset(accepted)

