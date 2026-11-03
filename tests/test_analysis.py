from pathlib import Path
import math

import pandas as pd
import pytest
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.text import Annotation

from asv_explorer.analysis import (
    REFERENCE_TAXON_COLORS,
    load_and_summarize,
    plot_nested_donut,
)


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


@pytest.mark.parametrize("order", [(0, 1, 2, 3), (3, 2, 1, 0), (0, 2, 1, 3)])
@pytest.mark.parametrize("categorical", [False, True])
def test_nested_donut_aligns_child_wedges_with_parents(
    tmp_path: Path, monkeypatch, order, categorical
):
    summary = pd.DataFrame(
        {
            "parent_taxon": ["Dinoflagellata", "Ochrophyta", "Dinoflagellata", "Ochrophyta"],
            "taxon": ["Gonyaulacales", "Thalassiosira", "Peridiniales", "Pseudo-nitzschia"],
            "abundance": [10.0, 20.0, 30.0, 40.0],
            "percentage": [10.0, 20.0, 30.0, 40.0],
        }
    ).iloc[list(order)]
    if categorical:
        summary["parent_taxon"] = pd.Categorical(
            summary["parent_taxon"],
            categories=["Ochrophyta", "Dinoflagellata", "Unused"],
            ordered=True,
        )
    original = summary.copy(deep=True)
    figures = []
    # Inspect real Matplotlib artists; skip only file serialization.
    monkeypatch.setattr(Figure, "savefig", lambda self, *args, **kwargs: figures.append(self))

    plot_nested_donut(summary, tmp_path / "interleaved.png")

    pd.testing.assert_frame_equal(summary, original)
    ax = figures[0].axes[0]
    parents = {w.get_label().split("\n")[0]: w for w in ax.patches if w.r == 0.75}
    children = [w for w in ax.patches if w.r == 1.12]
    annotations = [text for text in ax.texts if isinstance(text, Annotation)]
    assert list(parents) == list(summary["parent_taxon"].unique())
    assert len(children) == len(annotations) == len(summary)
    intervals = {name: [] for name in parents}
    for child in children:
        angle = math.radians((child.theta1 + child.theta2) / 2)
        anchor = (1.12 * math.cos(angle), 1.12 * math.sin(angle))
        matches = [a for a in annotations if a.xy == pytest.approx(anchor)]
        assert len(matches) == 1
        taxon = matches[0].get_text().split("  ")[0]
        row = summary.loc[summary["taxon"] == taxon].iloc[0]
        parent = parents[row["parent_taxon"]]
        assert child.theta1 >= parent.theta1 - 1e-4
        assert child.theta2 <= parent.theta2 + 1e-4
        assert child.theta2 - child.theta1 == pytest.approx(3.6 * row["percentage"], abs=1e-4)
        assert child.get_facecolor() == pytest.approx(to_rgba(REFERENCE_TAXON_COLORS[taxon]))
        intervals[row["parent_taxon"]].append((child.theta1, child.theta2))
    for name, spans in intervals.items():
        spans.sort()
        assert spans[0][0] == pytest.approx(parents[name].theta1, abs=1e-4)
        assert spans[-1][1] == pytest.approx(parents[name].theta2, abs=1e-4)
        for left, right in zip(spans, spans[1:]):
            assert left[1] == pytest.approx(right[0], abs=1e-4)
