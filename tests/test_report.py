"""Tests for ASV validation report generation."""

import json

import pandas as pd
import pytest

from asv_explorer.io.report import generate_validation_report


@pytest.fixture
def example_dataset(tmp_path):
    """Create a small, valid ASV dataset."""

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2", "ASV_3"],
        "Sample_A": [10, 20, 0],
        "Sample_B": [5, 0, 15],
        "Sample_Zero": [0, 0, 0],
    })

    taxonomy = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2", "ASV_3"],
        "Kingdom": ["Bacteria", "Bacteria", "Bacteria"],
        "Phylum": ["Pseudomonadota", "Bacteroidota", "Cyanobacteriota"],
    })

    counts_path = tmp_path / "counts.tsv"
    taxonomy_path = tmp_path / "taxonomy.tsv"
    outdir = tmp_path / "results"

    counts.to_csv(counts_path, sep="\t", index=False)
    taxonomy.to_csv(taxonomy_path, sep="\t", index=False)

    return counts_path, taxonomy_path, outdir


def test_report_generation(example_dataset):
    """Verify that all expected output files are generated."""

    counts_path, taxonomy_path, outdir = example_dataset

    generate_validation_report(counts_path, taxonomy_path, outdir)

    assert (outdir / "validation.json").exists()
    assert (outdir / "sample_summary.tsv").exists()
    assert (outdir / "feature_summary.tsv").exists()


def test_report_metadata(example_dataset):
    """Verify dataset dimensions and read counts."""

    counts_path, taxonomy_path, outdir = example_dataset

    report = generate_validation_report(
        counts_path, taxonomy_path, outdir
    )

    assert report["dataset"]["n_asvs"] == 3
    assert report["dataset"]["n_samples"] == 3
    assert report["dataset"]["total_reads"] == 50
    assert report["validation"]["taxonomy_alignment"] is True
    assert report["validation"]["read_conservation"] is True
    assert report["samples"]["zero_count"] == ["Sample_Zero"]


def test_sample_summary(example_dataset):
    """Verify read totals and observed ASVs per sample."""

    counts_path, taxonomy_path, outdir = example_dataset

    generate_validation_report(counts_path, taxonomy_path, outdir)

    summary = pd.read_csv(
        outdir / "sample_summary.tsv",
        sep="\t",
    ).set_index("sample_id")

    assert summary.loc["Sample_A", "total_reads"] == 30
    assert summary.loc["Sample_A", "observed_asvs"] == 2

    assert summary.loc["Sample_B", "total_reads"] == 20
    assert summary.loc["Sample_B", "observed_asvs"] == 2

    assert summary.loc["Sample_Zero", "total_reads"] == 0
    assert summary.loc["Sample_Zero", "observed_asvs"] == 0


def test_feature_summary(example_dataset):
    """Verify abundance and prevalence per ASV."""

    counts_path, taxonomy_path, outdir = example_dataset

    generate_validation_report(counts_path, taxonomy_path, outdir)

    summary = pd.read_csv(
        outdir / "feature_summary.tsv",
        sep="\t",
    ).set_index("ASV_ID")

    assert summary.loc["ASV_1", "total_reads"] == 15
    assert summary.loc["ASV_1", "prevalence"] == 2

    assert summary.loc["ASV_2", "total_reads"] == 20
    assert summary.loc["ASV_2", "prevalence"] == 1

    assert summary.loc["ASV_3", "total_reads"] == 15
    assert summary.loc["ASV_3", "prevalence"] == 1


def test_json_output(example_dataset):
    """Verify that the JSON output is readable and consistent."""

    counts_path, taxonomy_path, outdir = example_dataset

    report = generate_validation_report(
        counts_path, taxonomy_path, outdir
    )

    with open(
        outdir / "validation.json",
        encoding="utf-8",
    ) as handle:
        saved_report = json.load(handle)

    assert saved_report == report


def test_report_rejects_unmatched_taxonomy(example_dataset):
    """Reject datasets with missing taxonomic assignments."""

    counts_path, taxonomy_path, outdir = example_dataset

    taxonomy = pd.read_csv(taxonomy_path, sep="\t")
    taxonomy = taxonomy[taxonomy["ASV_ID"] != "ASV_3"]
    taxonomy.to_csv(taxonomy_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="Missing taxonomy"):
        generate_validation_report(
            counts_path, taxonomy_path, outdir
        )

    assert not (outdir / "validation.json").exists()