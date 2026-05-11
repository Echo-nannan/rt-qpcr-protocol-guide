"""Command-line interface for the RT-qPCR guide project."""

from __future__ import annotations

import argparse
from pathlib import Path

from rt_qpcr_guide.io import read_csv, write_dicts
from rt_qpcr_guide.qpcr_analysis import analyze_ddct, write_detail_csv, write_summary_csv
from rt_qpcr_guide.rt_calculator import plan_reactions


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rt-qpcr-guide")
    sub = parser.add_subparsers(dest="command", required=True)

    rt = sub.add_parser("rt-plan", help="Create a reverse-transcription volume plan.")
    rt.add_argument("--rna", required=True, help="RNA concentration CSV.")
    rt.add_argument("--target-ng", type=float, default=1000, help="RNA input per sample.")
    rt.add_argument("--output", required=True, help="Output CSV path.")
    rt.set_defaults(func=_cmd_rt_plan)

    analyze = sub.add_parser("analyze", help="Run Delta Delta Ct analysis.")
    analyze.add_argument("--ct", required=True, help="Ct values CSV.")
    analyze.add_argument("--samples", required=True, help="Sample information CSV.")
    analyze.add_argument("--reference", required=True, help="Reference gene, e.g. GAPDH.")
    analyze.add_argument("--control", required=True, help="Control group name.")
    analyze.add_argument("--outdir", required=True, help="Output directory.")
    analyze.set_defaults(func=_cmd_analyze)

    return parser


def _cmd_rt_plan(args: argparse.Namespace) -> None:
    reactions = plan_reactions(read_csv(args.rna), target_rna_ng=args.target_ng)
    write_dicts(
        args.output,
        [reaction.to_dict() for reaction in reactions],
        [
            "sample_id",
            "concentration_ng_ul",
            "target_rna_ng",
            "rna_volume_ul",
            "gdna_wiper_mix_ul",
            "water_to_16_ul",
            "rt_supermix_ul",
            "total_volume_ul",
            "status",
        ],
    )
    print(f"Wrote RT plan: {args.output}")


def _cmd_analyze(args: argparse.Namespace) -> None:
    result = analyze_ddct(
        read_csv(args.ct),
        read_csv(args.samples),
        reference_gene=args.reference,
        control_group=args.control,
    )
    outdir = Path(args.outdir)
    write_detail_csv(result.detail_rows, outdir / "detail.csv")
    write_summary_csv(result.summary_rows, outdir / "summary.csv")
    print(f"Wrote analysis results: {outdir}")


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

