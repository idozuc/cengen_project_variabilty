#!/usr/bin/env python3
"""Figure 8: recurrent WormCat category-3 enrichment among HVGs.

Heatmap of fold enrichment for categories enriched in >= 3 cell types,
faceted by neuron functional group. Blank cells are not significant.

Inputs: BASICS_DIR/analysis/wormcat_hvg_cengen_t2.csv
        (scripts/analysis/05_wormcat_cengen_filtered.py)
        BASICS_DIR/analysis/neuron_functional_groups.csv
        (scripts/analysis/01_prepare_mappings.R)
"""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import basics_file, save, set_style


GROUP_ORDER = ["Sensory", "Interneuron", "Motor", "Pharyngeal"]
NEUROPEPTIDE = "Neuronal function: synaptic function: neuropeptide"
MINIMUM_CELL_TYPES = 3

SHORT_LABELS = {
    NEUROPEPTIDE: "Neuropeptide",
    "Unassigned: regulated by multiple stresses": "Regulated by multiple stresses",
    "Transcription: chromatin structure: histone": "Histone / chromatin",
    "Stress response: detoxification: GST": "GST detoxification",
    "Ribosome: subunit": "Ribosomal subunit",
    "Metabolism: insulin": "Insulin",
    "Stress response: detoxification: ugt": "UGT detoxification",
    "Non-coding RNA: snRNA": "snRNA",
    "Stress response: unassigned": "Stress response, unassigned",
    "Transmembrane protein: seven transmembrane receptor: str": "7TM receptor: str",
    "Transmembrane transport: potassium channel": "Potassium channel",
    "Muscle function": "Muscle function",
    "Stress response: ER": "ER stress",
    "Globin": "Globin",
    "Extracellular material: galectin": "Galectin",
}

# The published mapping uses individual or alternate class names for these
# pooled CeNGEN labels. These assignments keep the corresponding class-level
# samples in the same functional groups.
GROUP_OVERRIDES = {
    "AS": "Motor",
    "AVB": "Interneuron",
    "AVD": "Interneuron",
    "AVE": "Interneuron",
    "DA": "Motor",
    "DB": "Motor",
    "PVC": "Interneuron",
    "RIV_stressed": "Motor",
    "RMD_LR": "Motor",
    "SMD_stressed": "Motor",
    "URA": "Motor",
    "VA": "Motor",
    "VB": "Motor",
    "VB02": "Motor",
}


def main() -> None:
    enrichment = pd.read_csv(basics_file("analysis", "wormcat_hvg_cengen_t2.csv"))
    enrichment = enrichment.loc[
        (enrichment["category_level"] == 3)
        & (enrichment["FDR"] <= 0.05)
        & (enrichment["fold_enrichment"] > 1)
    ].copy()
    if enrichment.duplicated(["wormcat_family", "cell_type"]).any():
        raise ValueError("Duplicate category/cell-type enrichment rows found")

    functional_groups = pd.read_csv(
        basics_file("analysis", "neuron_functional_groups.csv")
    )
    group_map = dict(zip(functional_groups["neuron"], functional_groups["neuron_functional_group"]))
    group_map.update(GROUP_OVERRIDES)
    enrichment["functional_group"] = enrichment["cell_type"].map(group_map)
    missing = sorted(enrichment.loc[enrichment["functional_group"].isna(), "cell_type"].unique())
    if missing:
        raise ValueError(f"Missing functional groups for: {missing}")

    ordered_cells = [
        cell
        for group in GROUP_ORDER
        for cell in sorted(enrichment.loc[enrichment["functional_group"] == group, "cell_type"].unique())
    ]
    fold_matrix = enrichment.pivot(
        index="wormcat_family", columns="cell_type", values="fold_enrichment"
    ).reindex(columns=ordered_cells)

    # Neuropeptide first, then by number of enriched cell types.
    frequency = enrichment.groupby("wormcat_family")["cell_type"].nunique()
    recurrent = frequency.loc[frequency >= MINIMUM_CELL_TYPES].index
    categories = [NEUROPEPTIDE] + sorted(
        [category for category in recurrent if category != NEUROPEPTIDE],
        key=lambda category: (-frequency[category], category),
    )

    color_map = plt.get_cmap("viridis").copy()
    color_map.set_bad("#F2F2F2")
    # Colour scale spans all significant enrichments, not only the ones shown.
    norm = mcolors.LogNorm(vmin=1, vmax=float(np.nanmax(fold_matrix.to_numpy())))

    set_style({"font.size": 8})
    fig = plt.figure(figsize=(12.0, 9.2))
    outer_grid = fig.add_gridspec(2, 1, hspace=0.48)
    top_grid = outer_grid[0].subgridspec(1, 2, width_ratios=[24, 20], wspace=0.08)
    bottom_grid = outer_grid[1].subgridspec(1, 2, width_ratios=[24, 6], wspace=0.08)
    axes = [
        fig.add_subplot(top_grid[0, 0]),
        fig.add_subplot(top_grid[0, 1]),
        fig.add_subplot(bottom_grid[0, 0]),
        fig.add_subplot(bottom_grid[0, 1]),
    ]

    for panel_index, (ax, group) in enumerate(zip(axes, GROUP_ORDER)):
        group_cells = [cell for cell in ordered_cells if group_map[cell] == group]
        panel = fold_matrix.loc[categories, group_cells].copy()
        panel.index = [SHORT_LABELS[category] for category in categories]

        sns.heatmap(
            panel,
            ax=ax,
            cmap=color_map,
            norm=norm,
            mask=panel.isna(),
            linewidths=0.4,
            linecolor="white",
            cbar=False,
        )
        ax.set_title(f"{group} (n={len(group_cells)} cell types)", fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=7.5)
        if panel_index % 2 == 0:
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
        else:
            ax.set_yticklabels([])

    fig.suptitle("Recurrent WormCat category-3 enrichment among HVGs", fontsize=12, y=0.985)
    fig.text(
        0.5,
        0.955,
        "Categories enriched in at least three cell types. "
        "Color: fold enrichment; blank: not significant.",
        ha="center",
        va="top",
        fontsize=9,
    )
    fig.supxlabel("Cell type", fontsize=9, y=0.025)
    fig.supylabel("WormCat category 3", fontsize=9, x=0.01)

    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=color_map),
        cax=fig.add_axes([0.925, 0.32, 0.014, 0.36]),
        ticks=[1, 2, 4, 8, 16],
    )
    colorbar.set_label("Fold enrichment")
    colorbar.ax.set_yticklabels(["1", "2", "4", "8", "16"])
    colorbar.outline.set_linewidth(0.5)

    fig.subplots_adjust(left=0.20, right=0.90, bottom=0.10, top=0.90)
    save(fig, "fig08_wormcat_category3_enrichment_heatmap.svg")


if __name__ == "__main__":
    main()
