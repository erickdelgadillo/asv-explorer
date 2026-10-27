"""Regression tests for validation, taxonomy grouping, and the installed CLI."""

from pathlib import Path

import pandas as pd
import pytest

from asv_explorer.analysis import load_and_summarize, plot_nested_donut
from asv_explorer.cli import main


def write_csv(tmp_path: Path, rows: list[dict]) -> Path:
    path = tmp_path / "input.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


@pytest.mark.parametrize(
    ("abundance", "message"),
    [
        (-0.1, "non-negative"),
        (None, "missing or non-numeric"),
        ("not-a-number", "missing or non-numeric"),
        (float("inf"), "infinite"),
        (float("-inf"), "infinite"),
        (0, "no positive values"),
    ],
)
def test_rejects_invalid_abundance(tmp_path: Path, abundance, message: str):
    path = write_csv(tmp_path, [{"parent_taxon": "A", "taxon": "a", "abundance": abundance}])
    with pytest.raises(ValueError, match=message):
        load_and_summarize(path)


def test_missing_columns_reported(tmp_path: Path):
    path = write_csv(tmp_path, [{"parent_taxon": "A", "taxon": "a"}])
    with pytest.raises(ValueError, match="Missing required column.*abundance"):
        load_and_summarize(path)


def test_custom_columns_and_missing_taxonomy(tmp_path: Path):
    path = write_csv(
        tmp_path,
        [
            {"Class": "A", "Taxonomy": "x", "TPM": 10},
            {"Class": "A", "Taxonomy": "x", "TPM": 20},
            {"Class": None, "Taxonomy": None, "TPM": 30},
            {"Class": "B", "Taxonomy": "zero", "TPM": 0},
        ],
    )
    result = load_and_summarize(
        path, parent_column="Class", taxon_column="Taxonomy", abundance_column="TPM"
    )
    assert len(result) == 2
    assert set(zip(result["parent_taxon"], result["taxon"], result["abundance"])) == {
        ("A", "x", 30),
        ("Unclassified", "Unclassified", 30),
    }
    assert result["percentage"].sum() == pytest.approx(100.0)


def test_cli_end_to_end_creates_expected_artifacts(tmp_path: Path):
    source = write_csv(
        tmp_path,
        [{"parent_taxon": "A", "taxon": "x", "abundance": 7},
         {"parent_taxon": "B", "taxon": "y", "abundance": 3}],
    )
    output = tmp_path / "nested" / "results"
    assert main(["--input", str(source), "--output-dir", str(output)]) == 0
    csv_file = output / "taxonomic_abundance_summary.csv"
    png_file = output / "taxonomic_abundance_donut.png"
    assert csv_file.is_file()
    assert png_file.is_file() and png_file.stat().st_size > 0
    summary = pd.read_csv(csv_file)
    assert summary["percentage"].sum() == pytest.approx(100.0)


def test_cli_rejects_negative_label_threshold(tmp_path: Path):
    with pytest.raises(SystemExit, match="non-negative"):
        main(["--min-label-percent", "-1", "--output-dir", str(tmp_path)])


def test_nested_donut_returns_output_path(tmp_path: Path):
    summary = pd.DataFrame(
        {"parent_taxon": ["A", "A", "B"], "taxon": ["x", "y", "z"],
         "abundance": [2.0, 3.0, 5.0], "percentage": [20.0, 30.0, 50.0]}
    )
    target = tmp_path / "plots" / "donut.png"
    assert plot_nested_donut(summary, target, min_label_percent=2) == target
    assert target.stat().st_size > 0
