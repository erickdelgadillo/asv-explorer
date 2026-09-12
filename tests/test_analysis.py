from pathlib import Path

import pandas as pd
import pytest

from asv_explorer.analysis import load_and_summarize, plot_nested_donut


def test_load_and_summarize_groups_rows_and_calculates_percentages(tmp_path: Path):
    input_path = tmp_path / "abundance.csv"
    pd.DataFrame(
        {
            "parent_taxon": ["A", "A", "B"],
            "taxon": ["x", "x", "y"],
            "abundance": [10, 20, 30],
        }
    ).to_csv(input_path, index=False)

    result = load_and_summarize(input_path)

    assert result["abundance"].tolist() == [30, 30]
    assert result["percentage"].sum() == pytest.approx(100)


def test_load_and_summarize_rejects_negative_abundance(tmp_path: Path):
    input_path = tmp_path / "abundance.csv"
    pd.DataFrame(
        {"parent_taxon": ["A"], "taxon": ["x"], "abundance": [-1]}
    ).to_csv(input_path, index=False)

    with pytest.raises(ValueError, match="non-negative"):
        load_and_summarize(input_path)


def test_plot_nested_donut_creates_png(tmp_path: Path):
    summary = pd.DataFrame(
        {
            "parent_taxon": ["A", "A", "B"],
            "taxon": ["x", "z", "y"],
            "abundance": [30.0, 10.0, 60.0],
            "percentage": [30.0, 10.0, 60.0],
        }
    )
    output_path = tmp_path / "chart.png"

    plot_nested_donut(summary, output_path)

    assert output_path.exists()
    assert output_path.stat().st_size > 0
