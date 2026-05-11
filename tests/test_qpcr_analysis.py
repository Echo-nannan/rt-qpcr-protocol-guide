import csv
from pathlib import Path

from rt_qpcr_guide.qpcr_analysis import analyze_ddct, write_summary_csv


def test_analyze_ddct_computes_fold_change():
    ct_rows = [
        {"sample_id": "C1", "gene": "GAPDH", "ct": "18.0", "replicate": "1"},
        {"sample_id": "C1", "gene": "GeneA", "ct": "23.0", "replicate": "1"},
        {"sample_id": "C2", "gene": "GAPDH", "ct": "18.0", "replicate": "1"},
        {"sample_id": "C2", "gene": "GeneA", "ct": "23.0", "replicate": "1"},
        {"sample_id": "T1", "gene": "GAPDH", "ct": "18.0", "replicate": "1"},
        {"sample_id": "T1", "gene": "GeneA", "ct": "21.0", "replicate": "1"},
        {"sample_id": "T2", "gene": "GAPDH", "ct": "18.0", "replicate": "1"},
        {"sample_id": "T2", "gene": "GeneA", "ct": "21.0", "replicate": "1"},
    ]
    sample_rows = [
        {"sample_id": "C1", "group": "Control"},
        {"sample_id": "C2", "group": "Control"},
        {"sample_id": "T1", "group": "Treatment"},
        {"sample_id": "T2", "group": "Treatment"},
    ]

    result = analyze_ddct(
        ct_rows,
        sample_rows,
        reference_gene="GAPDH",
        control_group="Control",
    )

    gene_summary = {row.gene: row for row in result.summary_rows}
    assert gene_summary["GeneA"].mean_delta_delta_ct == -2.0
    assert gene_summary["GeneA"].fold_change == 4.0
    assert gene_summary["GeneA"].regulation == "up"


def test_write_summary_csv(tmp_path: Path):
    result = analyze_ddct(
        [
            {"sample_id": "C1", "gene": "GAPDH", "ct": "18"},
            {"sample_id": "C1", "gene": "GeneA", "ct": "23"},
            {"sample_id": "T1", "gene": "GAPDH", "ct": "18"},
            {"sample_id": "T1", "gene": "GeneA", "ct": "22"},
        ],
        [
            {"sample_id": "C1", "group": "Control"},
            {"sample_id": "T1", "group": "Treatment"},
        ],
        reference_gene="GAPDH",
        control_group="Control",
    )
    out = tmp_path / "summary.csv"

    write_summary_csv(result.summary_rows, out)

    with out.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    treatment = next(row for row in rows if row["group"] == "Treatment")
    assert treatment["gene"] == "GeneA"
    assert treatment["fold_change"] == "2.0"
