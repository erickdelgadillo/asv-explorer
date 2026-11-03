"""Integration tests for the Nextflow ASV validation workflow."""

import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]

def test_nextflow_validation(tmp_path):
    """Run Nextflow and verify validation outputs."""

    assert shutil.which("nextflow") is not None, (
        "Nextflow is required for integration testing"
    )

    assert shutil.which("docker") is not None, (
        "Docker is required for integration testing"
    )

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2", "ASV_3"],
        "Sample_A": [10, 20, 0],
        "Sample_B": [5, 0, 15],
        "Sample_Zero": [0, 0, 0],
    })

    taxonomy = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2", "ASV_3"],
        "Kingdom": ["Bacteria", "Bacteria", "Bacteria"],
    })

    counts_path = tmp_path / "counts.tsv"
    taxonomy_path = tmp_path / "taxonomy.tsv"
    outdir = tmp_path / "results"

    counts.to_csv(counts_path, sep="\t", index=False)
    taxonomy.to_csv(taxonomy_path, sep="\t", index=False)

    command = [
        "nextflow",
        "run",
        str(PROJECT_ROOT / "main.nf"),
        "--counts",
        str(counts_path),
        "--taxonomy",
        str(taxonomy_path),
        "--outdir",
        str(outdir),
    "-with-docker",
    ]

    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )

    assert completed.returncode == 0, (
        f"Nextflow failed:\n"
        f"STDOUT:\n{completed.stdout}\n"
        f"STDERR:\n{completed.stderr}"
    )

    assert (outdir / "validation.json").exists()
    assert (outdir / "sample_summary.tsv").exists()
    assert (outdir / "feature_summary.tsv").exists()

    with open(
        outdir / "validation.json",
        encoding="utf-8",
    ) as handle:
        report = json.load(handle)

    assert report["dataset"]["n_asvs"] == 3
    assert report["dataset"]["n_samples"] == 3
    assert report["dataset"]["total_reads"] == 50
    assert report["samples"]["zero_count"] == ["Sample_Zero"]
    assert report["validation"]["read_conservation"] is True