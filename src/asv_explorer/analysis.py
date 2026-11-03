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

# Palette used in bacterial-suppression/results/figures/main/4.2.Metatranscriptomes.png.
# Hex values are the exact equivalents of the R colour names in my.colors_All.
REFERENCE_TAXON_COLORS = {
    "Other Spirotrichea": "#473C8B",
    "Spirotrichea": "#473C8B",
    "Choreotrichida": "#6A5ACD",
    "Tintinnida": "#836FFF",
    "Other Dinoflagellata": "#000080",
    "Dinoflagellata": "#000080",
    "Dinophyceae": "#000080",
    "Gonyaulacales": "#00C5CD",
    "Alexandrium": "#00F5FF",
    "Gymnodiniales": "#87CEFA",
    "Peridiniales": "#4169E1",
    "Suessiales": "#6CA6CD",
    "Symbiodinium": "#87CEFF",
    "Syndiniales": "#CAE1FF",
    "Other Stramenopiles": "#006400",
    "Stramenopiles": "#006400",
    "Other Ochrophyta": "#008B45",
    "Ochrophyta": "#008B45",
    "Other Bacillariophyta": "#90EE90",
    "Bacillariophyta": "#90EE90",
    "Thalassiosira": "#00CD00",
    "Chaetoceros": "#CAFF70",
    "Skeletonema": "#008B45",
    "Pseudo-nitzschia": "#00FF00",
    "Glaucophyta": "#CD3700",
    "Other Chlorophyta": "#295345",
    "Chlorophyta": "#295345",
    "Micromonas": "#7FFFD4",
    "Ostreococcus": "#66CDAA",
    "Amoebozoa": "#8B5F65",
    "Hacrobia": "#CDC673",
    "Haptophyta": "#CDC673",
    "Rhizaria": "#CDB5CD",
    "Other eukaryotes": "#545454",
    "Unclassified": "#A7ADB4",
}


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


def _format_abundance(value: float) -> str:
    """Format a total abundance compactly for the centre of the chart."""
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:.0f}"


def _spread_label_positions(
    labels: list[dict[str, float | int]],
    *,
    minimum_gap: float = 0.15,
    lower: float = -1.16,
    upper: float = 1.16,
) -> list[dict[str, float | int]]:
    """Spread labels vertically while preserving their circular order."""
    labels = sorted(labels, key=lambda item: float(item["target_y"]))
    if not labels:
        return labels

    positions = [float(item["target_y"]) for item in labels]
    for index in range(1, len(positions)):
        positions[index] = max(positions[index], positions[index - 1] + minimum_gap)

    if positions[-1] > upper:
        shift = positions[-1] - upper
        positions = [position - shift for position in positions]
    for index in range(len(positions) - 2, -1, -1):
        positions[index] = min(positions[index], positions[index + 1] - minimum_gap)
    if positions[0] < lower:
        shift = lower - positions[0]
        positions = [position + shift for position in positions]

    for item, position in zip(labels, positions):
        item["label_y"] = position
    return labels


def plot_nested_donut(
    summary: pd.DataFrame,
    output_path: str | Path,
    *,
    title: str = "Taxonomic relative abundance",
    min_label_percent: float = 1.0,
) -> Path:
    """Plot two taxonomic rings with collision-reduced labels on the outer ring."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Keep each parent's children contiguous in both rings, colors and labels.
    # Preserve encounter order within and between groups without mutating input.
    summary = pd.concat(
        [group for _, group in summary.groupby("parent_taxon", sort=False, observed=True)],
        ignore_index=True,
    )
    parent = (
        summary.groupby("parent_taxon", sort=False, observed=True, as_index=False)["abundance"]
        .sum()
    )
    parent_percentages = parent["abundance"] / parent["abundance"].sum() * 100
    parent_labels = [
        f"{name}\n{percentage:.1f}%"
        if percentage >= max(min_label_percent, 3.0)
        else ""
        for name, percentage in zip(parent["parent_taxon"], parent_percentages)
    ]
    fallback_palette = plt.get_cmap("tab20")
    parent_colors = [
        REFERENCE_TAXON_COLORS.get(name, fallback_palette(index % 20)[:3])
        for index, name in enumerate(parent["parent_taxon"])
    ]
    color_by_parent = dict(zip(parent["parent_taxon"], parent_colors))

    child_colors: list[tuple[float, float, float]] = []
    for parent_name, group in summary.groupby("parent_taxon", sort=False, observed=True):
        base = color_by_parent[parent_name]
        count = len(group)
        positions = np.linspace(0, 1, count, endpoint=False) if count > 1 else [0.25]
        for (_, row), position in zip(group.iterrows(), positions):
            child_colors.append(
                REFERENCE_TAXON_COLORS.get(
                    row["taxon"], _shade(base, float(position))
                )
            )

    fig, donut_ax = plt.subplots(figsize=(13, 11))
    donut_ax.pie(
        parent["abundance"],
        radius=0.75,
        labels=parent_labels,
        colors=parent_colors,
        labeldistance=0.55,
        startangle=90,
        counterclock=False,
        textprops={"fontsize": 9, "fontweight": "bold"},
        wedgeprops={"width": 0.34, "edgecolor": "white", "linewidth": 1.5},
    )
    outer_wedges, _ = donut_ax.pie(
        summary["abundance"],
        radius=1.12,
        labels=None,
        colors=child_colors,
        startangle=90,
        counterclock=False,
        wedgeprops={"width": 0.35, "edgecolor": "white", "linewidth": 1.0},
    )

    left_labels: list[dict[str, float | int]] = []
    right_labels: list[dict[str, float | int]] = []
    for index, (wedge, row) in enumerate(zip(outer_wedges, summary.itertuples())):
        if row.percentage < min_label_percent:
            continue
        angle = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)
        item: dict[str, float | int] = {
            "index": index,
            "angle": float(angle),
            "target_y": float(np.sin(angle) * 1.24),
        }
        (right_labels if np.cos(angle) >= 0 else left_labels).append(item)

    for side, labels in ((1, right_labels), (-1, left_labels)):
        for item in _spread_label_positions(labels):
            index = int(item["index"])
            angle = float(item["angle"])
            row = summary.iloc[index]
            duplicate_taxon = (summary["taxon"] == row["taxon"]).sum() > 1
            label_name = (
                f"{row['parent_taxon']} › {row['taxon']}"
                if duplicate_taxon
                else row["taxon"]
            )
            donut_ax.annotate(
                f"{label_name}  {row['percentage']:.1f}%",
                xy=(np.cos(angle) * 1.12, np.sin(angle) * 1.12),
                xytext=(side * 1.48, float(item["label_y"])),
                ha="left" if side > 0 else "right",
                va="center",
                fontsize=9,
                color="#263442",
                arrowprops={
                    "arrowstyle": "-",
                    "color": child_colors[index],
                    "linewidth": 1.1,
                    "connectionstyle": "arc3,rad=0",
                },
            )

    total_abundance = float(summary["abundance"].sum())
    donut_ax.text(
        0,
        0.03,
        _format_abundance(total_abundance),
        ha="center",
        va="center",
        fontsize=24,
        fontweight="bold",
        color="#243447",
    )
    donut_ax.text(
        0,
        -0.12,
        f"total abundance\n{len(summary)} taxa",
        ha="center",
        va="center",
        fontsize=10,
        color="#607080",
        linespacing=1.35,
    )
    donut_ax.set_title("Taxonomic composition", fontsize=13, fontweight="bold", pad=15)
    donut_ax.set(aspect="equal")
    fig.suptitle(title, fontsize=19, fontweight="bold", y=0.97, color="#17212b")
    fig.text(
        0.5,
        0.925,
        "Inner ring: higher taxonomic groups  •  Outer ring: individual taxa",
        ha="center",
        fontsize=10,
        color="#607080",
    )
    fig.text(
        0.5,
        0.045,
        f"Taxa below {min_label_percent:g}% remain visible in the outer ring but are not labelled.",
        ha="center",
        fontsize=9,
        color="#6f7d89",
    )
    fig.subplots_adjust(left=0.08, right=0.92, top=0.88, bottom=0.09)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output_path
