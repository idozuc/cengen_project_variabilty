#!/usr/bin/env python3
"""Figures 13-15: GMM clusters and one gene's expression on the stored UMAP.

Point colour = GMM cluster, point size = log1p(raw UMI) of the gene.
Figure 13: HSN, tph-1.  Figure 14: SMD, flp-5.  Figure 15: NSM, del-3.

Input: OUTPUT_DIR/data/gene_umap_expression.csv (04_export_gene_umap_data.R). The UMAP
coordinates are the ones stored by the clustering run, not recomputed here.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

from common import cluster_palette, data_file, save, set_style


FIGURES = [("fig13", "HSN", "tph-1"), ("fig14", "SMD", "flp-5"), ("fig15", "NSM", "del-3")]
MIN_SIZE, SIZE_RANGE = 8, 62  # point area = 8 + 62 * expression / max expression


def plot_gene(frame: pd.DataFrame, figure: str, cell_type: str, gene: str) -> None:
    frame = frame.sort_values("log_count", kind="stable")
    palette = cluster_palette(frame["cluster"])
    maximum = float(frame["log_count"].max())

    def area(value):
        return MIN_SIZE + (SIZE_RANGE * value / maximum if maximum > 0 else 0)

    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    ax.scatter(
        frame["UMAP1"],
        frame["UMAP2"],
        c=frame["cluster"].map(palette),
        s=area(frame["log_count"]),
        alpha=0.75,
        edgecolors="none",
        rasterized=True,
        zorder=2,
    )

    positive = int((frame["raw_count"] > 0).sum())
    total = len(frame)
    ax.set_title(f"{positive:,}/{total:,} cells detected ({100 * positive / total:.1f}%)", fontsize=9.5, pad=8)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, left=True, bottom=True)

    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#666666",
                   markeredgecolor="none", markersize=np.sqrt(area(value)), alpha=0.75,
                   label=f"{value:.2g}")
            for value in np.unique(np.array([0.0, maximum / 2, maximum]))
        ],
        title="Gene expression\nlog1p(raw UMI)",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=7.5,
        title_fontsize=8,
        labelspacing=1.0,
    )
    cluster_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=colour,
               markeredgecolor="none", markersize=6, label=f"Cluster {cluster}")
        for cluster, colour in palette.items()
    ]
    fig.legend(
        handles=cluster_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.89),
        ncol=min(len(cluster_handles), 4),
        frameon=False,
        handlelength=2.2,
        columnspacing=1.2,
        title="GMM cluster",
        fontsize=8,
        title_fontsize=8,
    )
    fig.suptitle(
        f"{gene} expression in {cell_type}",
        fontsize=12,
        fontweight="bold",
        fontstyle="italic",
        y=0.985,
    )
    fig.subplots_adjust(left=0.11, right=0.91, bottom=0.10, top=0.74)
    save(
        fig,
        f"{figure}_{cell_type}_{gene}_expression_cluster_color_expression_size_umap.svg",
        dpi=300,
        facecolor="white",
    )


def main() -> None:
    data = pd.read_csv(data_file("gene_umap_expression.csv"), dtype={"cell": str, "cluster": str})
    set_style({"font.size": 8, "axes.linewidth": 0.7})
    for figure, cell_type, gene in FIGURES:
        frame = data[(data.cell_type == cell_type) & (data.gene == gene)]
        if frame.empty:
            raise ValueError(f"No exported rows for {cell_type}/{gene}")
        plot_gene(frame.copy(), figure, cell_type, gene)


if __name__ == "__main__":
    main()
