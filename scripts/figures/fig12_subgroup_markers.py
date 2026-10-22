#!/usr/bin/env python3
"""Figure 12: marker detection per cluster in the HSN, SMD and NSM partitions.

Dumbbell plot of the fraction of cells with >= 1 UMI for each gene in each
cluster of the accepted solution.

Input: OUTPUT_DIR/data/subgroup_markers.csv (03_export_subgroup_markers.R).
"""

import matplotlib.pyplot as plt
import pandas as pd

from common import data_file, save


COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
MARKERS = ["o", "s", "D"]  # second encoding so clusters never rely on colour alone
INK, MUTED, GRID = "#2b2b2b", "#6b6b6b", "#e4e4e0"
CELL_TYPES = ["HSN", "SMD", "NSM"]


def main() -> None:
    d = pd.read_csv(data_file("subgroup_markers.csv"))

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": INK,
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0))

    for ax, letter, cell_type in zip(axes, "ABC", CELL_TYPES):
        x = d[d.cell_type == cell_type]
        genes = list(dict.fromkeys(x.gene))[::-1]
        y = {gene: i for i, gene in enumerate(genes)}
        for gene in genes:
            v = x[x.gene == gene].detection_fraction
            ax.plot([v.min(), v.max()], [y[gene], y[gene]], color=GRID, lw=2, zorder=1,
                    solid_capstyle="round")
        n_clusters = x.cluster.nunique()
        for k, (cluster, sub) in enumerate(x.groupby("cluster")):
            dodge = (k - (n_clusters - 1) / 2) * 0.22  # keep equal values visible
            ax.scatter(sub.detection_fraction, sub.gene.map(y) + dodge, s=30, marker=MARKERS[k],
                       color=COLORS[k], edgecolor="white", linewidth=0.8, zorder=2,
                       label=f"cluster {cluster} (n={sub.cells.iloc[0]:,})")
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontstyle="italic")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.7, len(genes) - 0.3)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xlabel("Fraction of cells detected")
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_title(cell_type, color=INK, fontsize=9, loc="left")
        ax.text(-0.42, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", color=INK)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, fontsize=7,
                  handletextpad=0.2, labelcolor=INK)

    fig.tight_layout(w_pad=2.5)
    save(fig, "fig12_subgroup_markers.pdf")


if __name__ == "__main__":
    main()
