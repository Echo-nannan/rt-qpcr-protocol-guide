"""Smoke test for Pfaffl method (v2.5).

Verifies:
    * efficiency normalisation handles base / pct / 0~1 inputs
    * standard ΔΔCt fold-change with ΔΔCt = -1 → 2.0
    * Pfaffl with E=1.8 and same ΔΔCt = -1 → 1.8
    * default fall-back: gene without explicit eff → 2.0
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd
from complete_gui import DeltaCtCalculator, _normalise_efficiency


def make_df() -> pd.DataFrame:
    np.random.seed(0)
    rows = []
    for sample, group, dct_offset in [
        ("S1", "C", 0.0), ("S2", "C", 0.0),
        ("S3", "T", -1.0), ("S4", "T", -1.0),
    ]:
        for gene, base in [("GAPDH", 18.0), ("GeneA", 24.0), ("GeneB", 26.0)]:
            ct = base + (dct_offset if gene in {"GeneA", "GeneB"} else 0.0)
            for _ in range(2):
                rows.append({
                    "Sample": sample, "Group": group, "Gene": gene,
                    "Ct": float(ct + np.random.normal(0, 0.0005)),
                })
    return pd.DataFrame(rows)


def main() -> None:
    print("=== efficiency normalisation ===")
    print("base 1.95:", _normalise_efficiency(1.95))
    print("pct 95:", _normalise_efficiency(95))
    print("pct 100:", _normalise_efficiency(100))
    print("delta 0.85:", _normalise_efficiency(0.85))
    print("string '1.95':", _normalise_efficiency("1.95"))
    print("invalid 'abc':", _normalise_efficiency("abc"))

    df = make_df()

    print("\n=== ΔΔCt method (default) ===")
    calc = DeltaCtCalculator()
    calc.raw_data = df
    calc.calculate(
        sample_col="Sample", gene_col="Gene", ct_col="Ct",
        ref_gene="GAPDH", ctrl_sample="C", group_col="Group",
    )
    fc_col = "2^-ΔΔCt"

    def _pick(summary: pd.DataFrame, group: str, gene: str) -> float:
        mask = (summary["Group"].astype(str) == group) & (summary["Gene"].astype(str) == gene)
        vals = summary.loc[mask, fc_col].dropna()
        return float(vals.mean()) if not vals.empty else float("nan")

    fc = _pick(calc.summary, "T", "GeneA")
    print(f"GeneA T fold = {fc:.4f}  (expected ≈ 2.0)")
    assert abs(fc - 2.0) < 0.01

    print("\n=== Pfaffl method (GeneA E=1.8, GeneB defaults to 2.0) ===")
    calc2 = DeltaCtCalculator()
    calc2.raw_data = df
    calc2.calculate(
        sample_col="Sample", gene_col="Gene", ct_col="Ct",
        ref_gene="GAPDH", ctrl_sample="C", group_col="Group",
        method="pfaffl",
        gene_efficiencies={"GeneA": 1.8},
    )
    fc_a = _pick(calc2.summary, "T", "GeneA")
    fc_b = _pick(calc2.summary, "T", "GeneB")
    print(f"GeneA T fold = {fc_a:.4f}  (expected ≈ 1.8)")
    print(f"GeneB T fold = {fc_b:.4f}  (expected ≈ 2.0, default base)")
    assert abs(fc_a - 1.8) < 0.01
    assert abs(fc_b - 2.0) < 0.01
    print(f"calc2.last_method={calc2.last_method!r} eff={calc2.last_gene_efficiencies}")

    print("\n=== Pfaffl method (efficiency from percent: GeneA=95) ===")
    calc3 = DeltaCtCalculator()
    calc3.raw_data = df
    calc3.calculate(
        sample_col="Sample", gene_col="Gene", ct_col="Ct",
        ref_gene="GAPDH", ctrl_sample="C", group_col="Group",
        method="pfaffl",
        gene_efficiencies={"GeneA": 95},  # 95% → base 1.95
    )
    fc95 = _pick(calc3.summary, "T", "GeneA")
    print(f"GeneA T fold = {fc95:.4f}  (expected ≈ 1.95)")
    assert abs(fc95 - 1.95) < 0.01

    print("\n[OK] Pfaffl smoke test passed")


if __name__ == "__main__":
    main()
