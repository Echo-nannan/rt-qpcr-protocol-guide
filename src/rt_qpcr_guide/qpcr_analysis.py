"""Delta Delta Ct analysis using only the Python standard library."""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class DetailRow:
    sample_id: str
    group: str
    gene: str
    mean_ct: float
    ct_sd: float
    delta_ct: float
    delta_delta_ct: float
    fold_change: float
    log2fc: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SummaryRow:
    gene: str
    group: str
    n_samples: int
    mean_delta_ct: float
    mean_delta_delta_ct: float
    fold_change: float
    log2fc: float
    regulation: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class DDCTResult:
    detail_rows: list[DetailRow]
    summary_rows: list[SummaryRow]


def _mean(values: list[float]) -> float:
    return statistics.mean(values)


def _sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def _round(value: float) -> float:
    return round(value, 4)


def _sample_groups(sample_rows: list[dict[str, str]]) -> dict[str, str]:
    groups: dict[str, str] = {}
    for row in sample_rows:
        sample_id = row.get("sample_id", "").strip()
        group = row.get("group", "").strip()
        if sample_id and group:
            groups[sample_id] = group
    return groups


def _mean_ct_by_sample_gene(ct_rows: list[dict[str, str]]) -> dict[tuple[str, str], tuple[float, float]]:
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in ct_rows:
        sample_id = row.get("sample_id", "").strip()
        gene = row.get("gene", "").strip()
        ct_text = row.get("ct", "").strip()
        if not sample_id or not gene or not ct_text:
            continue
        ct_value = float(ct_text)
        if math.isfinite(ct_value):
            values[(sample_id, gene)].append(ct_value)

    return {key: (_mean(cts), _sd(cts)) for key, cts in values.items()}


def analyze_ddct(
    ct_rows: list[dict[str, str]],
    sample_rows: list[dict[str, str]],
    reference_gene: str,
    control_group: str,
) -> DDCTResult:
    groups = _sample_groups(sample_rows)
    mean_ct = _mean_ct_by_sample_gene(ct_rows)
    reference_gene = reference_gene.strip()
    control_group = control_group.strip()

    delta_ct: dict[tuple[str, str], float] = {}
    ct_sd: dict[tuple[str, str], float] = {}
    for (sample_id, gene), (target_mean, target_sd) in mean_ct.items():
        if gene == reference_gene:
            continue
        ref_key = (sample_id, reference_gene)
        if ref_key not in mean_ct:
            continue
        delta_ct[(sample_id, gene)] = target_mean - mean_ct[ref_key][0]
        ct_sd[(sample_id, gene)] = target_sd

    control_delta_by_gene: dict[str, list[float]] = defaultdict(list)
    for (sample_id, gene), value in delta_ct.items():
        if groups.get(sample_id) == control_group:
            control_delta_by_gene[gene].append(value)

    control_mean_by_gene = {
        gene: _mean(values)
        for gene, values in control_delta_by_gene.items()
        if values
    }

    detail_rows: list[DetailRow] = []
    for (sample_id, gene), value in sorted(delta_ct.items()):
        if gene not in control_mean_by_gene or sample_id not in groups:
            continue
        ddct = value - control_mean_by_gene[gene]
        fold = 2 ** (-ddct)
        detail_rows.append(
            DetailRow(
                sample_id=sample_id,
                group=groups[sample_id],
                gene=gene,
                mean_ct=_round(mean_ct[(sample_id, gene)][0]),
                ct_sd=_round(ct_sd[(sample_id, gene)]),
                delta_ct=_round(value),
                delta_delta_ct=_round(ddct),
                fold_change=_round(fold),
                log2fc=_round(-ddct),
            )
        )

    by_gene_group: dict[tuple[str, str], list[DetailRow]] = defaultdict(list)
    for row in detail_rows:
        by_gene_group[(row.gene, row.group)].append(row)

    summary_rows: list[SummaryRow] = []
    for (gene, group), rows in sorted(by_gene_group.items()):
        mean_delta_ct = _mean([row.delta_ct for row in rows])
        mean_ddct = _mean([row.delta_delta_ct for row in rows])
        fold = 2 ** (-mean_ddct)
        regulation = "up" if fold > 1.2 else "down" if fold < 0.8 else "stable"
        summary_rows.append(
            SummaryRow(
                gene=gene,
                group=group,
                n_samples=len(rows),
                mean_delta_ct=_round(mean_delta_ct),
                mean_delta_delta_ct=_round(mean_ddct),
                fold_change=_round(fold),
                log2fc=_round(-mean_ddct),
                regulation=regulation,
            )
        )

    return DDCTResult(detail_rows=detail_rows, summary_rows=summary_rows)


def write_detail_csv(rows: list[DetailRow], path: str | Path) -> None:
    _write_dataclass_csv(rows, path, list(DetailRow.__dataclass_fields__.keys()))


def write_summary_csv(rows: list[SummaryRow], path: str | Path) -> None:
    _write_dataclass_csv(rows, path, list(SummaryRow.__dataclass_fields__.keys()))


def _write_dataclass_csv(rows: list[DetailRow] | list[SummaryRow], path: str | Path, fieldnames: list[str]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.to_dict())

