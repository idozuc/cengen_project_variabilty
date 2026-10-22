#!/usr/bin/env python3
"""Figures 10-11: marker expression on the ASE and AWC positive-control UMAPs.

Cells are coloured by log1p(raw UMI) of each known marker; outlines enclose
the central 95% of each GMM cluster. The UMAP is the same as in Figure 9.

Inputs: OUTPUT_DIR/data/positive_ASE.csv, data/positive_AWC.csv
        (02_export_positive_control_pcs.R)
        OUTPUT_DIR/data/expression_histogram_data.csv.gz (01_export_expression_data.R)
"""

from __future__ import annotations

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

from common import calculate_umap, cluster_palette, data_file, ordered_labels, save, set_style


DATASETS = {
    "ASE": {
        "figure": "fig10",
        "coordinate_file": "positive_ASE.csv",
        "display_name": "Pooled ASEL + ASER",
        "markers": ["gcy-6", "gcy-5", "gcy-7", "gcy-22"],
        "layout": (2, 2),
    },
    "AWC": {
        "figure": "fig11",
        "coordinate_file": "positive_AWC.csv",
        "display_name": "Pooled AWC_ON + AWC_OFF",
        "markers": ["str-2", "srsx-3"],
        "layout": (1, 2),
    },
}

EXPRESSION_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "marker_expression",
    ["#D9D9D9", "#FEE6CE", "#F16913", "#7F0000"],
)
BOUNDARY_FRACTION = 0.95


def cluster_boundaries(coordinates: pd.DataFrame) -> dict[str, np.ndarray]:
    """Convex hull around the 95% of each cluster's cells closest to its median."""
    boundaries = {}
    for cluster in ordered_labels(coordinates["cluster"]):
        members = coordinates.loc[
            coordinates["cluster"].astype(str) == cluster, ["UMAP1", "UMAP2"]
        ].to_numpy(dtype=float)
        if len(members) < 3:
            continue
        distances = np.linalg.norm(members - np.median(members, axis=0), axis=1)
        keep = max(3, int(np.ceil(BOUNDARY_FRACTION * len(members))))
        central = members[np.argsort(distances)[:keep]]
        polygon = central[ConvexHull(central).vertices]
        boundaries[cluster] = np.vstack([polygon, polygon[0]])
    return boundaries


def marker_frame(coordinates: pd.DataFrame, expression: pd.DataFrame, dataset: str, gene: str) -> pd.DataFrame:
    marker = expression.loc[
        (expression["dataset"] == dataset) & (expression["gene_name"] == gene),
        ["cell", "identity", "raw_count", "log_count"],
    ]
    if set(marker["cell"]) != set(coordinates["cell"]) or marker["cell"].duplicated().any():
        raise ValueError(f"{dataset}/{gene}: expression and UMAP cells differ")

    joined = coordinates.merge(marker, on="cell", how="left", validate="one_to_one")
    if not (joined["truth"] == joined["identity"]).all():
        raise ValueError(f"{dataset}/{gene}: identity labels do not match")
    # Draw zero-expression cells first and the strongest expressers on top.
    return joined.sort_values("log_count", kind="stable")


def marker_panel(ax, data: pd.DataFrame, gene: str, boundaries, palette) -> None:
    norm = mcolors.Normalize(vmin=0, vmax=max(float(data["log_count"].max()), 1e-12))
    points = ax.scatter(
        data["UMAP1"],
        data["UMAP2"],
        c=data["log_count"],
        cmap=EXPRESSION_CMAP,
        norm=norm,
        s=10,
        alpha=0.90,
        edgecolors="none",
        rasterized=True,
    )
    for cluster, polygon in boundaries.items():
        # A white underlay keeps the hull legible over strongly expressing cells.
        ax.plot(polygon[:, 0], polygon[:, 1], color="white", linewidth=2.4, zorder=4)
        ax.plot(polygon[:, 0], polygon[:, 1], color=palette[cluster], linewidth=1.25, zorder=5)

    positive = int((data["raw_count"] > 0).sum())
    total = len(data)
    ax.set_title(
        f"{gene}\n{positive:,}/{total:,} cells detected ({100 * positive / total:.1f}%)",
        fontsize=9.5,
        fontstyle="italic",
        pad=7,
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    sns.despine(ax=ax, left=True, bottom=True)

    colorbar = ax.figure.colorbar(points, ax=ax, fraction=0.048, pad=0.025)
    colorbar.set_label("log1p(raw UMI)", fontsize=7.5)
    colorbar.ax.tick_params(labelsize=7, length=2)
    colorbar.outline.set_linewidth(0.5)


def plot_dataset(dataset: str, spec: dict, expression: pd.DataFrame) -> None:
    coordinates = pd.read_csv(
        data_file(spec["coordinate_file"]),
        dtype={"cell": str, "cluster": str, "truth": str},
    )
    embedding = calculate_umap(coordinates)
    coordinates["UMAP1"] = embedding[:, 0]
    coordinates["UMAP2"] = embedding[:, 1]
    palette = cluster_palette(coordinates["cluster"])
    boundaries = cluster_boundaries(coordinates)

    rows, columns = spec["layout"]
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(3.7 * columns, 3.8 if rows == 1 else 3.35 * rows),
        squeeze=False,
    )
    for ax, gene in zip(axes.flat, spec["markers"]):
        marker_panel(ax, marker_frame(coordinates, expression, dataset, gene), gene, boundaries, palette)

    fig.legend(
        handles=[
            Line2D([0], [0], color=colour, linewidth=1.5, label=f"Cluster {cluster}")
            for cluster, colour in palette.items()
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925 if rows == 1 else 0.945),
        ncol=len(palette),
        frameon=False,
        handlelength=2.2,
        columnspacing=1.2,
        title="GMM cluster boundary (central 95% of cells)",
        fontsize=8,
        title_fontsize=8,
    )
    fig.suptitle(
        f"Marker expression in {spec['display_name']}",
        fontsize=12,
        fontweight="bold",
        y=0.985,
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.96,
        bottom=0.08,
        top=0.70 if rows == 1 else 0.86,
        wspace=0.34,
        hspace=0.34,
    )
    save(
        fig,
        f"{spec['figure']}_{dataset.lower()}_marker_expression_umaps.svg",
        dpi=300,
        facecolor="white",
    )


def main() -> None:
    expression = pd.read_csv(data_file("expression_histogram_data.csv.gz"), dtype={"cell": str})
    set_style({"font.size": 8, "axes.linewidth": 0.7})
    for dataset, spec in DATASETS.items():
        plot_dataset(dataset, spec, expression)


if __name__ == "__main__":
    main()
