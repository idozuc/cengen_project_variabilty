#!/usr/bin/env python3
"""Figure 9: clustering of pooled populations with known neuronal identities.

PCA (PC1-PC2) and UMAP of ASEL + ASER and AWC_ON + AWC_OFF, coloured by
known identity or by GMM cluster. UMAP is computed here from the selected
PCs (seed 20260809), the same embedding as in Figures 10-11.

Input: OUTPUT_DIR/data/positive_ASE.csv, data/positive_AWC.csv and
       OUTPUT_DIR/data/positive_control_manifest.csv (02_export_positive_control_pcs.R).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import calculate_umap, cluster_palette, data_file, ordered_labels, save, set_style


IDENTITY_COLORS = {
    "ASEL": "#0072B2",
    "ASER": "#D55E00",
    "AWC_ON": "#009E73",
    "AWC_OFF": "#CC79A7",
}
ORDER = ["positive_ASE", "positive_AWC"]


def point_size(number_of_cells: int) -> float:
    for threshold, size in ((2200, 3.2), (1400, 4.0), (800, 5.2), (400, 6.5)):
        if number_of_cells > threshold:
            return size
    return 8.0


def scatter_categories(ax, x, y, labels: pd.Series, palette, *, show_legend, legend_prefix) -> None:
    labels = labels.astype(str)
    size = point_size(len(labels))
    for label in ordered_labels(labels):
        selected = (labels == label).to_numpy()
        ax.scatter(
            np.asarray(x)[selected],
            np.asarray(y)[selected],
            s=size,
            alpha=0.72,
            color=palette[label],
            edgecolors="none",
            rasterized=True,
            label=f"{legend_prefix} {label} (n={int(selected.sum()):,})",
        )
    if show_legend:
        ax.legend(
            frameon=False,
            fontsize=7,
            markerscale=1.8,
            handletextpad=0.35,
            borderaxespad=0.25,
            loc="best",
        )


def style_projection_axis(ax, xlabel: str, ylabel: str, *, umap_axis: bool = False) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal", adjustable="datalim")
    ax.tick_params(labelsize=7, length=2.5, width=0.6)
    if umap_axis:
        ax.set_xticks([])
        ax.set_yticks([])
    sns.despine(ax=ax)


def pc_label(component: str, value: float) -> str:
    return component if pd.isna(value) else f"{component} ({100 * value:.1f}%)"


def statistics(row: pd.Series) -> str:
    return (
        f"N={int(row.n_cells):,}   "
        f"G={int(row.primary_g)}   "
        f"ARI vs known identity={row.external_ari:.3f}   "
        f"ΔBIC vs G=1={row.bic_delta_from_one:.1f}   "
        f"median stability ARI={row.stability_median:.3f}"
    )


def main() -> None:
    manifest = pd.read_csv(data_file("positive_control_manifest.csv")).set_index("dataset_id")

    set_style({"font.size": 8, "axes.linewidth": 0.7})
    fig, axes = plt.subplots(2, 4, figsize=(12.4, 6.5))
    fig.suptitle(
        "Positive-control recovery of known neuronal identities",
        fontsize=12,
        fontweight="bold",
        y=0.992,
    )

    for row_index, dataset_id in enumerate(ORDER):
        row = manifest.loc[dataset_id]
        frame = pd.read_csv(
            data_file(row.data_file), dtype={"cell": str, "cluster": str, "truth": str}
        )
        if len(frame) != int(row.n_cells):
            raise ValueError(f"{dataset_id}: coordinate rows do not match manifest N")
        embedding = calculate_umap(frame)

        truth_palette = {label: IDENTITY_COLORS[label] for label in ordered_labels(frame["truth"])}
        gmm_palette = cluster_palette(frame["cluster"])
        panels = [
            (frame["PC1"], frame["PC2"], "truth", truth_palette, True, "", "PC1–PC2: known identity"),
            (frame["PC1"], frame["PC2"], "cluster", gmm_palette, True, "Cluster", "PC1–PC2: GMM cluster"),
            (embedding[:, 0], embedding[:, 1], "truth", truth_palette, False, "", "UMAP: known identity"),
            (embedding[:, 0], embedding[:, 1], "cluster", gmm_palette, False, "Cluster", "UMAP: GMM cluster"),
        ]
        for column, (x, y, label_column, palette, legend, prefix, title) in enumerate(panels):
            ax = axes[row_index, column]
            scatter_categories(
                ax, x, y, frame[label_column], palette,
                show_legend=legend, legend_prefix=prefix,
            )
            ax.set_title(title, fontsize=9)
            if column < 2:
                style_projection_axis(
                    ax,
                    pc_label("PC1", row.pc1_variance_fraction),
                    pc_label("PC2", row.pc2_variance_fraction),
                )
            else:
                style_projection_axis(ax, "UMAP 1", "UMAP 2", umap_axis=True)

    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.07, top=0.86, wspace=0.30, hspace=0.58)
    for position, dataset_id in zip([0.915, 0.485], ORDER):
        row = manifest.loc[dataset_id]
        fig.text(
            0.065,
            position,
            f"{row.display_name}   {statistics(row)}",
            ha="left",
            va="center",
            fontsize=8.5,
            fontweight="bold",
        )

    save(fig, "fig09_positive_controls_pca_umap.svg", dpi=300, facecolor="white")


if __name__ == "__main__":
    main()
