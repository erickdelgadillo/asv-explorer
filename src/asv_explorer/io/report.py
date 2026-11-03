"""Generate validation reports for ASV datasets."""

import json
from pathlib import Path

import pandas as pd

from asv_explorer.io.dada2 import (
    load_dada2_counts,
    load_dada2_taxonomy,
    find_zero_count_samples,
    validate_asv_alignment,
)


def generate_validation_report(
    counts_path: str | Path,
    taxonomy_path: str | Path,
    outdir: str | Path,
) -> dict:
    """Validate ASV inputs and generate summary reports."""

    counts_path = Path(counts_path)
    taxonomy_path = Path(taxonomy_path)
    outdir = Path(outdir)

    # 1. Load and validate inputs
    counts = load_dada2_counts(counts_path)
    taxonomy = load_dada2_taxonomy(taxonomy_path)

    validate_asv_alignment(counts, taxonomy)

    # 2. Extract abundance matrix
    abundance = counts.drop(columns="ASV_ID")

    n_asvs = len(counts)
    n_samples = len(abundance.columns)
    total_reads = int(abundance.to_numpy().sum())

    # 3. Generate sample summary
    sample_summary = pd.DataFrame({
        "sample_id": abundance.columns,
        "total_reads": abundance.sum(axis=0).to_numpy(),
        "observed_asvs": (abundance > 0).sum(axis=0).to_numpy(),
    })

    # 4. Generate feature summary
    feature_summary = pd.DataFrame({
        "ASV_ID": counts["ASV_ID"],
        "total_reads": abundance.sum(axis=1),
        "prevalence": (abundance > 0).sum(axis=1),
    })

    # 5. Verify read conservation
    sample_total = int(sample_summary["total_reads"].sum())
    feature_total = int(feature_summary["total_reads"].sum())

    if not (total_reads == sample_total == feature_total):
        raise ValueError("Read conservation check failed")

    # 6. Build validation report
    report = {
        "dataset": {
            "n_asvs": n_asvs,
            "n_samples": n_samples,
            "total_reads": total_reads,
        },
        "validation": {
            "taxonomy_alignment": True,
            "read_conservation": True,
        },
        "samples": {
            "zero_count": find_zero_count_samples(counts),
        },
        "provenance": {
            "counts_file": str(counts_path.resolve()),
            "taxonomy_file": str(taxonomy_path.resolve()),
        },
    }

    # 7. Write outputs
    outdir.mkdir(parents=True, exist_ok=True)

    sample_summary.to_csv(
        outdir / "sample_summary.tsv",
        sep="\t",
        index=False,
    )

    feature_summary.to_csv(
        outdir / "feature_summary.tsv",
        sep="\t",
        index=False,
    )

    with open(outdir / "validation.json", "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=4)

    return report