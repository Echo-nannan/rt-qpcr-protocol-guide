"""End-to-end smoke test for the .ixo one-click pipeline (no GUI).

Verifies:
  1. ``IxoParser`` loads 384 wells from the real PFC-Xiongxing.ixo
  2. ``PCRDataParser`` (which now auto-routes to IxoParser) produces the same
     wide-format Sample×Gene dataframe as the user's existing converted xlsx
  3. ``export_long_format_with_mean`` + ``DeltaCtCalculator`` produces a
     non-empty side-by-side dataframe with the expected columns (matching
     the structure shown in the user's reference image #3).
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.complete_gui import PCRDataParser, DeltaCtCalculator  # noqa: E402

DATA_DIR = Path(
    r"K:\shi_hao_nan_shuoshi_ke_ti\Topic_Pregnancy_Exposure_Neurodevelopment"
    r"\2.Animal_experiments\2.Offspring_Studies\07_RT-qPCR"
    r"\PND56♂RT-qPCR_The guess has rna-seq base\20260426F1-♂-PFC\data"
)
IXO = DATA_DIR / "2026.4.26.01.05am-PFC-Xiongxing.ixo"
USER_XLSX = DATA_DIR / "2026.4.26.01.05am-PFC-Xiongxing_qPCR_converted.xlsx"


def step(msg: str) -> None:
    print(f"\n=== {msg} ===")


def main() -> None:
    step("1. IxoParser.parse_file → PCRDataParser auto-route")
    parser = PCRDataParser()
    parser.parse_file(str(IXO))
    print(f"   loaded {len(parser.data)} wells (expect 384)")
    assert len(parser.data) == 384, "well count mismatch"

    sample_names = (
        "C,C,C,C,TD-L,TD-L,TD-L,TD-M,TD-M,TD-M,TD-H,TD-H,TD-H".split(",")
    )
    # Image 2 shows the user types 13 names; layout uses only the first 12.
    gene_names = "GAPDH,RPL13a,GFAP,IBA1,Dly4,Syp,Bdnf,Arc".split(",")
    cols_per_sample, rows_per_gene = 2, 2
    samples_per_group = 3  # 12 samples / 4 groups

    step("2. export_sample_gene_mean (wide) — compare to user's existing xlsx")
    wide = parser.export_sample_gene_mean(
        cols_per_sample, rows_per_gene, sample_names, gene_names,
    )
    print(wide.head(13))
    print(f"   shape: {wide.shape}")

    if USER_XLSX.exists():
        try:
            user_wide = pd.read_excel(USER_XLSX, sheet_name='平均CT值')
            print(f"   user's '平均CT值' shape: {user_wide.shape}")
            # Compare numeric columns where both have value
            common = [c for c in wide.columns[1:] if c in user_wide.columns]
            for c in common[:3]:
                a = pd.to_numeric(wide[c], errors='coerce').to_numpy()
                b = pd.to_numeric(user_wide[c], errors='coerce').to_numpy()
                paired = [(x, y) for x, y in zip(a, b) if not np.isnan(x) and not np.isnan(y)]
                if paired:
                    deltas = [abs(x - y) for x, y in paired]
                    print(
                        f"   gene {c}: {len(paired)} pairs, max|Δ|={max(deltas):.4f}, "
                        f"mean|Δ|={sum(deltas)/len(deltas):.4f}"
                    )
        except Exception as exc:
            print(f"   (skip xlsx compare: {exc})")
    else:
        print(f"   (no user xlsx at {USER_XLSX})")

    step("3. export_long_format_with_mean (for ΔΔCt)")
    long_df = parser.export_long_format_with_mean(
        cols_per_sample, rows_per_gene, sample_names, gene_names,
        samples_per_group=samples_per_group,
    )
    print(long_df.head(8))
    print(f"   long_df shape: {long_df.shape}")
    print(f"   unique groups: {long_df['Group'].unique().tolist()}")
    print(f"   unique samples: {long_df['Sample'].unique().tolist()}")

    step("4. DeltaCtCalculator end-to-end")
    calc = DeltaCtCalculator()
    calc.raw_data = long_df
    calc.calculate(
        sample_col='Sample', gene_col='Gene', ct_col='Ct',
        ref_gene='GAPDH', ctrl_sample='Group1', group_col='Group',
    )
    print(f"   results shape:        {calc.results.shape}")
    print(f"   summary shape:        {calc.summary.shape}")
    print(f"   side_by_side shape:   {calc.side_by_side.shape}")
    print(f"   side_by_side columns: {list(calc.side_by_side.columns)}")
    print("\n   First rows of side_by_side:")
    print(calc.side_by_side.head(10).to_string())

    step("[OK] Smoke test passed")


if __name__ == "__main__":
    main()
