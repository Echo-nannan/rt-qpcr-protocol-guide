"""Run a small end-to-end demo with the bundled CSV templates."""

from __future__ import annotations

from pathlib import Path

from rt_qpcr_guide.cli import main


ROOT = Path(__file__).resolve().parents[1]


def run() -> None:
    main(
        [
            "rt-plan",
            "--rna",
            str(ROOT / "templates" / "rna_concentration_template.csv"),
            "--output",
            str(ROOT / "results" / "demo_rt_plan.csv"),
        ]
    )
    main(
        [
            "analyze",
            "--ct",
            str(ROOT / "templates" / "ct_values_template.csv"),
            "--samples",
            str(ROOT / "templates" / "sample_info_template.csv"),
            "--reference",
            "GAPDH",
            "--control",
            "Control",
            "--outdir",
            str(ROOT / "results" / "demo_ddct"),
        ]
    )


if __name__ == "__main__":
    run()

