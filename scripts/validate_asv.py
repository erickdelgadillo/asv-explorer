#!/usr/bin/env python3

"""Command-line interface for ASV input validation."""

import argparse

from asv_explorer.io.report import generate_validation_report


def main():
    parser = argparse.ArgumentParser(
        description="Validate ASV counts and taxonomy tables."
    )

    parser.add_argument(
        "--counts",
        required=True,
        help="Path to the ASV count table.",
    )

    parser.add_argument(
        "--taxonomy",
        required=True,
        help="Path to the taxonomy table.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory for validation outputs.",
    )

    args = parser.parse_args()

    report = generate_validation_report(
        counts_path=args.counts,
        taxonomy_path=args.taxonomy,
        outdir=args.outdir,
    )

    print("ASV validation completed successfully.")
    print(f"ASVs: {report['dataset']['n_asvs']}")
    print(f"Samples: {report['dataset']['n_samples']}")
    print(f"Total reads: {report['dataset']['total_reads']}")


if __name__ == "__main__":
    main()