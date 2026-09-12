"""Load, validate, summarize, and visualize taxonomic abundance data."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

_matplotlib_cache = Path(tempfile.gettempdir()) / "asv-explorer-matplotlib"
_matplotlib_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_matplotlib_cache))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb


REQUIRED_OUTPUT_COLUMNS = ["parent_taxon", "taxon", "abundance", "percentage"]


def load_and_summarize(
    input_path: str | Path,
    *,
    parent_column: str = "parent_taxon",
    taxon_column: str = "taxon",
    abundance_column: str = "abundance",
) -> pd.DataFrame:
    """Read a CSV abundance table and return one validated row per taxon."""
    input_path = Path(input_path)
    table = pd.read_csv(input_path)

    required = {parent_column, taxon_column, abundance_column}
    missing = sorted(required.difference(table.columns))
    if missing:
        raise ValueError(
            f"Missing required column(s): {', '.join(missing)}. "
            f"Available columns: {', '.join(table.columns)}"
        )

    selected = table[[parent_column, taxon_column, abundance_column]].copy()
    selected.columns = ["parent_taxon", "taxon", "abundance"]
    selected["parent_taxon"] = selected["parent_taxon"].fillna("Unclassified").astype(str)
    selected["taxon"] = selected["taxon"].fillna("Unclassified").astype(str)
    selected["abundance"] = pd.to_numeric(selected["abundance"], errors="coerce")

    if selected["abundance"].isna().any():
        raise ValueError("The abundance column contains missing or non-numeric values.")
    if not np.isfinite(selected["abundance"]).all():
        raise ValueError("The abundance column contains infinite values.")
    if (selected["abundance"] < 0).any():
        raise ValueError("Abundance values must be non-negative.")

    summary = (
        selected.groupby(["parent_taxon", "taxon"], as_index=False, sort=False)["abundance"]
        .sum()
        .loc[lambda frame: frame["abundance"] > 0]
    )
    total = summary["abundance"].sum()
    if total <= 0:
        raise ValueError("The abundance table has no positive values to plot.")

    summary["percentage"] = summary["abundance"] / total * 100
    parent_order = (
        summary.groupby("parent_taxon")["abundance"]
        .sum()
        .sort_values(ascending=False)
        .index
    )
    summary["parent_taxon"] = pd.Categorical(
        summary["parent_taxon"], categories=parent_order, ordered=True
    )
    summary = summary.sort_values(
        ["parent_taxon", "abundance"], ascending=[True, False], kind="stable"
    ).reset_index(drop=True)
    summary["parent_taxon"] = summary["parent_taxon"].astype(str)
    return summary[REQUIRED_OUTPUT_COLUMNS]


def _shade(color: tuple[float, float, float], position: float) -> tuple[float, float, float]:
    """Blend a base color toward white for distinguishable child taxa."""
    rgb = np.asarray(to_rgb(color))
    amount = 0.08 + 0.52 * position
    return tuple(rgb + (1 - rgb) * amount)


def _autopct(minimum: float):
    def format_percentage(value: float) -> str:
        return f"{value:.1f}%" if value >= minimum else ""

    return format_percentage


def plot_nested_donut(
    summary: pd.DataFrame,
    output_path: str | Path,
    *,
    title: str = "Taxonomic relative abundance",
    min_label_percent: float = 1.0,
) -> Path:
    """Plot parent taxa in the inner ring and child taxa in the outer ring."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    parent = (
        summary.groupby("parent_taxon", sort=False, as_index=False)["abundance"]
        .sum()
    )
    parent_percentages = parent["abundance"] / parent["abundance"].sum() * 100
    parent_labels = [
        name if percentage >= max(min_label_percent, 2.0) else ""
        for name, percentage in zip(parent["parent_taxon"], parent_percentages)
    ]
    palette = plt.get_cmap("tab20")
    parent_colors = [palette(i % 20)[:3] for i in range(len(parent))]
    color_by_parent = dict(zip(parent["parent_taxon"], parent_colors))

    child_colors: list[tuple[float, float, float]] = []
    for parent_name, group in summary.groupby("parent_taxon", sort=False):
        base = color_by_parent[parent_name]
        count = len(group)
        positions = np.linspace(0, 1, count, endpoint=False) if count > 1 else [0.25]
        child_colors.extend(_shade(base, float(position)) for position in positions)

    labels = [
        taxon if percentage >= min_label_percent else ""
        for taxon, percentage in zip(summary["taxon"], summary["percentage"])
    ]

    fig, ax = plt.subplots(figsize=(13, 10), constrained_layout=True)
    ax.pie(
        parent["abundance"],
        radius=0.78,
        labels=parent_labels,
        colors=parent_colors,
        labeldistance=0.52,
        startangle=90,
        counterclock=False,
        textprops={"fontsize": 9},
        wedgeprops={"width": 0.36, "edgecolor": "white"},
    )
    ax.pie(
        summary["abundance"],
        radius=1.15,
        labels=labels,
        colors=child_colors,
        labeldistance=1.05,
        pctdistance=0.86,
        autopct=_autopct(min_label_percent),
        startangle=90,
        counterclock=False,
        rotatelabels=True,
        textprops={"fontsize": 8},
        wedgeprops={"width": 0.36, "edgecolor": "white"},
    )
    ax.set_title(title, fontsize=15, pad=22)
    ax.set(aspect="equal")
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output_path
