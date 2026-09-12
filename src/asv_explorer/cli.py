"""Command-line interface for taxonomic abundance visualization."""

from __future__ import annotations

import argparse
from pathlib import Path

from .analysis import load_and_summarize, plot_nested_donut


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize a taxonomic abundance CSV and create a nested-donut chart."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/example_taxonomic_abundance.csv"),
        help="Input CSV (default: data/example_taxonomic_abundance.csv).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for the summary and chart (default: results).",
    )
    parser.add_argument("--parent-column", default="parent_taxon")
    parser.add_argument("--taxon-column", default="taxon")
    parser.add_argument("--abundance-column", default="abundance")
    parser.add_argument(
        "--min-label-percent",
        type=float,
        default=1.0,
        help="Hide labels and percentages below this threshold (default: 1.0).",
    )
    parser.add_argument(
        "--title", default="Relative abundance of transcripts by taxonomy"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_label_percent < 0:
        raise SystemExit("--min-label-percent must be non-negative")

    summary = load_and_summarize(
        args.input,
        parent_column=args.parent_column,
        taxon_column=args.taxon_column,
        abundance_column=args.abundance_column,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "taxonomic_abundance_summary.csv"
    chart_path = args.output_dir / "taxonomic_abundance_donut.png"
    summary.to_csv(summary_path, index=False)
    plot_nested_donut(
        summary,
        chart_path,
        title=args.title,
        min_label_percent=args.min_label_percent,
    )
    print(f"Summary: {summary_path}")
    print(f"Chart:   {chart_path}")
    return 0
