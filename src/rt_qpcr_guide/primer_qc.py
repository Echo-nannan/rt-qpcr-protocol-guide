"""Simple primer quality checks for qPCR planning."""

from __future__ import annotations

from dataclasses import dataclass


_VALID_BASES = {"A", "T", "G", "C"}


@dataclass(frozen=True)
class PrimerQC:
    sequence: str
    length: int
    gc_percent: float
    wallace_tm_c: int
    warnings: list[str]

    @property
    def passes_basic_rules(self) -> bool:
        return not self.warnings


@dataclass(frozen=True)
class PrimerPair:
    forward: str
    reverse: str
    product_size_bp: int


@dataclass(frozen=True)
class PrimerPairQC:
    forward: PrimerQC
    reverse: PrimerQC
    product_size_bp: int
    product_size_ok: bool
    tm_difference_ok: bool
    warnings: list[str]

    @property
    def passes(self) -> bool:
        return self.forward.passes_basic_rules and self.reverse.passes_basic_rules and not self.warnings


def normalize_sequence(sequence: str) -> str:
    seq = sequence.replace(" ", "").replace("\n", "").upper()
    if not seq:
        raise ValueError("primer sequence is empty")
    invalid = sorted(set(seq) - _VALID_BASES)
    if invalid:
        raise ValueError(f"invalid bases in primer sequence: {', '.join(invalid)}")
    return seq


def check_primer(sequence: str) -> PrimerQC:
    seq = normalize_sequence(sequence)
    length = len(seq)
    gc_count = seq.count("G") + seq.count("C")
    gc_percent = round(gc_count / length * 100, 2)
    wallace_tm = 2 * (seq.count("A") + seq.count("T")) + 4 * gc_count

    warnings: list[str] = []
    if not 18 <= length <= 24:
        warnings.append("Length outside 18-24 bp")
    if not 40 <= gc_percent <= 60:
        warnings.append("GC content outside 40-60%")
    if not 58 <= wallace_tm <= 64:
        warnings.append("Wallace Tm outside 58-64 degC")
    if "GGGG" in seq or "CCCC" in seq:
        warnings.append("Contains >=4 consecutive G/C bases")

    return PrimerQC(
        sequence=seq,
        length=length,
        gc_percent=gc_percent,
        wallace_tm_c=wallace_tm,
        warnings=warnings,
    )


def check_primer_pair(pair: PrimerPair) -> PrimerPairQC:
    forward = check_primer(pair.forward)
    reverse = check_primer(pair.reverse)
    product_size_ok = 80 <= pair.product_size_bp <= 200
    tm_difference_ok = abs(forward.wallace_tm_c - reverse.wallace_tm_c) <= 2

    warnings: list[str] = []
    if not product_size_ok:
        warnings.append("Product size outside 80-200 bp")
    if not tm_difference_ok:
        warnings.append("Forward/reverse Tm difference > 2 degC")

    return PrimerPairQC(
        forward=forward,
        reverse=reverse,
        product_size_bp=pair.product_size_bp,
        product_size_ok=product_size_ok,
        tm_difference_ok=tm_difference_ok,
        warnings=warnings,
    )

