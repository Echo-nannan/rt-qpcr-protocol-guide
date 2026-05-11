"""Smoke test for the new wide→ΔΔCt tab logic, using the user's real xlsx.

Directly runs the same pipeline that the GUI uses internally:
    user xlsx ('原始数据' sheet) → unique sample IDs → long DF → DeltaCtCalculator
                                                              → side_by_side
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402

from src.complete_gui import DeltaCtCalculator  # noqa: E402

USER_XLSX = Path(
    r"K:\shi_hao_nan_shuoshi_ke_ti\Topic_Pregnancy_Exposure_Neurodevelopment"
    r"\2.Animal_experiments\2.Offspring_Studies\07_RT-qPCR"
    r"\PND56♂RT-qPCR_The guess has rna-seq base\20260426F1-♂-PFC\data"
    r"\2026.4.26.01.05am-PFC-Xiongxing_qPCR_converted.xlsx"
)


def step(msg: str) -> None:
    print(f"\n=== {msg} ===")


def build_long_from_wide(
    wide: pd.DataFrame, sample_col: str, ref_gene: str, target_genes: list[str]
) -> pd.DataFrame:
    """Replicates the GUI's _build_wide_long_df logic without Tk."""
    keep = [ref_gene] + [g for g in target_genes if g != ref_gene]
    long_rows = []
    seen: dict[str, int] = {}
    for _, row in wide.iterrows():
        if pd.isna(row[sample_col]):
            continue
        group = str(row[sample_col]).strip()
        if not group:
            continue
        seen[group] = seen.get(group, 0) + 1
        sample_id = f"{group}_{seen[group]:02d}"
        for gene in keep:
            if gene not in wide.columns:
                continue
            v = pd.to_numeric(row[gene], errors='coerce')
            if pd.notna(v):
                long_rows.append({
                    'Sample': sample_id, 'Group': group, 'Gene': gene, 'Ct': float(v),
                })
    return pd.DataFrame(long_rows)


def main() -> None:
    step("1. Load user's wide-format xlsx")
    xl = pd.ExcelFile(USER_XLSX)
    print(f"   sheets: {xl.sheet_names}")
    wide = pd.read_excel(USER_XLSX, sheet_name='原始数据')
    print(f"   '原始数据' shape: {wide.shape}")
    print(f"   columns: {list(wide.columns)}")
    print(f"   first 3 rows:\n{wide.head(3)}")

    sample_col = 'Sample'
    ref_gene = 'GAPDH'
    target_genes = [c for c in wide.columns if c not in {sample_col, ref_gene}]
    print(f"   ref={ref_gene}, targets={target_genes}")

    step("2. Build long DF (unique sample IDs per group)")
    long_df = build_long_from_wide(wide, sample_col, ref_gene, target_genes)
    print(f"   long_df shape: {long_df.shape}")
    print(f"   groups: {long_df['Group'].unique().tolist()}")
    print(f"   sample IDs (first 6): {long_df['Sample'].unique().tolist()[:6]}")

    step("3. ΔΔCt with control group = 'C'")
    calc = DeltaCtCalculator()
    calc.raw_data = long_df
    calc.calculate(
        sample_col='Sample', gene_col='Gene', ct_col='Ct',
        ref_gene=ref_gene, ctrl_sample='C', group_col='Group',
    )
    print(f"   detailed results: {calc.results.shape}")
    print(f"   summary: {calc.summary.shape}")
    print(f"   side_by_side: {calc.side_by_side.shape}")

    print("\n   First side_by_side rows (showing Bdnf section if present):")
    sbs = calc.side_by_side
    bdnf_idx = sbs.index[sbs['Target Name.1'] == 'Bdnf'].tolist()
    if bdnf_idx:
        print(sbs.iloc[bdnf_idx[0]: bdnf_idx[0] + 5].to_string())
    else:
        print(sbs.head(5).to_string())

    print("\n   Summary preview (Bdnf rows):")
    summ = calc.summary
    if 'Gene' in summ.columns:
        bdnf_summ = summ[summ['Gene'] == 'Bdnf'].head(8)
        if not bdnf_summ.empty:
            print(bdnf_summ.to_string())

    step("[OK] Wide-tab smoke test passed")


if __name__ == "__main__":
    main()
