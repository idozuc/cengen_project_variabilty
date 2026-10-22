#!/usr/bin/env python3
"""Figure 7: gene-family enrichment in HVGs.

A: neuropeptide share of all eligible genes versus HVGs in RMD_DV (donuts).
B: number of cell types with significant enrichment for each WormCat
   category-3 term enriched in >= 3 cell types.

Input: BASICS_DIR/analysis/wormcat_hvg_cengen_t2.csv
       (scripts/analysis/05_wormcat_cengen_filtered.py; HVGs and background
       restricted to genes CeNGEN calls expressed at threshold 2).
"""

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import basics_file, save, set_style


CELL_TYPE = "RMD_DV"
NEUROPEPTIDE = "Neuronal function: synaptic function: neuropeptide"
NEUROPEPTIDE_COLOR = "#D55E00"
OTHER_COLOR = "#D9D9D9"
MINIMUM_CELL_TYPES = 3


def format_scientific(value: float) -> str:
    exponent = int(np.floor(np.log10(value)))
    return rf"${value / 10**exponent:.2f} \times 10^{{{exponent}}}$"


def draw_donut(ax, neuropeptide_count, total_count, title, center_label):
    ax.pie(
        [neuropeptide_count, total_count - neuropeptide_count],
        colors=[NEUROPEPTIDE_COLOR, OTHER_COLOR],
        startangle=90,
        counterclock=False,
        radius=1,
        wedgeprops={"width": 0.36, "edgecolor": "white", "linewidth": 1.0},
    )
    ax.text(0, 0, center_label, ha="center", va="center", fontsize=10, linespacing=1.25)
    ax.set_title(title, fontsize=11, pad=10)
    ax.text(
        0,
        -1.20,
        f"Neuropeptide: {neuropeptide_count:,} "
        f"({100 * neuropeptide_count / total_count:.2f}%)",
        ha="center",
        va="top",
        fontsize=9,
        color=NEUROPEPTIDE_COLOR,
        fontweight="bold",
    )
    ax.set_aspect("equal")


def neuropeptide_donuts(enrichment: pd.DataFrame) -> None:
    selected = enrichment.loc[
        (enrichment["cell_type"] == CELL_TYPE)
        & (enrichment["category_level"] == 3)
        & (enrichment["wormcat_family"] == NEUROPEPTIDE)
    ]
    if len(selected) != 1:
        raise ValueError(f"Expected one {CELL_TYPE} neuropeptide row, found {len(selected)}")

    row = selected.iloc[0]
    # WormCat columns: AC/n_background = category/all background genes,
    # RGS/n_hvg = category/all HVGs.
    background_np, background_total = int(row["AC"]), int(row["n_background"])
    hvg_np, hvg_total = int(row["RGS"]), int(row["n_hvg"])
    fold_enrichment = (hvg_np / hvg_total) / (background_np / background_total)
    if not np.isclose(fold_enrichment, row["fold_enrichment"]):
        raise ValueError("Recalculated fold enrichment does not match the source table")

    set_style({"font.size": 9}, theme="white")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7))
    draw_donut(
        axes[0], background_np, background_total,
        "All analysis-eligible genes", f"{background_total:,}\neligible genes",
    )
    draw_donut(axes[1], hvg_np, hvg_total, "Highly variable genes", f"{hvg_total:,}\nHVGs")

    fig.suptitle(f"Neuropeptide enrichment in {CELL_TYPE}", fontsize=13, y=0.98)
    fig.text(
        0.5,
        0.55,
        f"{fold_enrichment:.2f}× enrichment\n\nFDR = {format_scientific(row['FDR'])}",
        ha="center",
        va="center",
        fontsize=10,
    )
    fig.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, color=NEUROPEPTIDE_COLOR),
            plt.Rectangle((0, 0), 1, 1, color=OTHER_COLOR),
        ],
        labels=["Neuropeptide category", "Other genes"],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.subplots_adjust(left=0.04, right=0.96, bottom=0.20, top=0.80, wspace=0.20)
    save(fig, "fig07A_neuropeptide_enrichment_RMD_DV_donuts.svg")


def enriched_cell_type_counts(enrichment: pd.DataFrame) -> None:
    significant = enrichment.loc[
        (enrichment["category_level"] == 3)
        & (enrichment["FDR"] <= 0.05)
        & (enrichment["fold_enrichment"] > 1)
    ]
    counts = (
        significant.groupby("wormcat_family")
        .agg(enriched_cell_types=("cell_type", "nunique"))
        .sort_values("enriched_cell_types")
        .reset_index()
    )
    counts = counts.loc[counts["enriched_cell_types"] >= MINIMUM_CELL_TYPES]
    labels = counts["wormcat_family"].map(lambda text: "\n".join(textwrap.wrap(text, width=45)))

    set_style({"font.size": 9})
    fig, ax = plt.subplots(figsize=(6.5, max(5, 0.32 * len(counts))))
    bars = ax.barh(labels, counts["enriched_cell_types"], color="#4C78A8", edgecolor="none")
    ax.bar_label(bars, padding=3, fontsize=8)
    ax.set_xlabel("Cell types with significant enrichment")
    ax.set_ylabel("WormCat category 3")
    ax.set_xlim(0, counts["enriched_cell_types"].max() * 1.1)
    sns.despine()
    fig.tight_layout()
    save(fig, "fig07B_wormcat_category3_enriched_cell_types.svg")


def main() -> None:
    enrichment = pd.read_csv(basics_file("analysis", "wormcat_hvg_cengen_t2.csv"))
    neuropeptide_donuts(enrichment)
    enriched_cell_type_counts(enrichment)


if __name__ == "__main__":
    main()
