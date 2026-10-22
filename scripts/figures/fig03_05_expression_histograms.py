#!/usr/bin/env python3
"""Expression histograms.

Figure 3 (A-F): known ASE/AWC markers in the pooled positive controls.
Figure 5A: matched HVG/LVG examples in URX.
Figure 5B: high-mean, low-zero HVGs in AVK and RIG.

Input: OUTPUT_DIR/data/expression_histogram_data.csv.gz and
       OUTPUT_DIR/data/expression_histogram_manifest.csv (01_export_expression_data.R).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import PercentFormatter

from common import data_file, save, set_style


HVG_COLOR = "#D55E00"
LVG_COLOR = "#0072B2"
IDENTITY_COLORS = {
    "ASEL": "#0072B2",
    "ASER": "#E69F00",
    "AWC_ON": "#009E73",
    "AWC_OFF": "#CC79A7",
}

# Figure 3 panel order follows the manuscript caption.
MARKER_PANELS = [
    ("A", "ASE", "gcy-5", ["ASEL", "ASER"], "ASEL + ASER"),
    ("B", "ASE", "gcy-7", ["ASEL", "ASER"], "ASEL + ASER"),
    ("C", "ASE", "gcy-22", ["ASEL", "ASER"], "ASEL + ASER"),
    ("D", "ASE", "gcy-6", ["ASEL", "ASER"], "ASEL + ASER"),
    ("E", "AWC", "srsx-3", ["AWC_ON", "AWC_OFF"], "AWC_ON + AWC_OFF"),
    ("F", "AWC", "str-2", ["AWC_ON", "AWC_OFF"], "AWC_ON + AWC_OFF"),
]


def log_bins(values: pd.Series) -> np.ndarray:
    """0.25-wide bins centred on multiples of 0.25 of log1p(count)."""
    return np.arange(-0.125, np.ceil(values.max() / 0.25) * 0.25 + 0.375, 0.25)


def percentage_histogram(
    axis: plt.Axes,
    values: np.ndarray,
    bins: np.ndarray,
    color: str,
    *,
    label: str | None = None,
    histtype: str = "bar",
    alpha: float = 0.82,
    linewidth: float = 0.4,
) -> None:
    axis.hist(
        values,
        bins=bins,
        weights=np.full(values.size, 1.0 / values.size),
        color=color,
        label=label,
        histtype=histtype,
        alpha=alpha,
        linewidth=linewidth,
        edgecolor="white" if histtype == "bar" else color,
    )


def basics_label(row) -> str:
    return (
        f"μ = {row.Mu:.2f}\n"
        f"ε = {row.Epsilon:.2f}\n"
        f"zeros = {row.pct_zero_master:.1f}%"
    )


def marker_histograms(data: pd.DataFrame) -> None:
    for panel, dataset, gene, identities, pooled_label in MARKER_PANELS:
        subset = data[data["dataset"].eq(dataset) & data["gene_name"].eq(gene)]
        bins = log_bins(subset["log_count"])

        fig, axis = plt.subplots(figsize=(3.6, 2.8))
        for identity in identities:
            identity_data = subset[subset["identity"].eq(identity)]
            zero_fraction = 100 * identity_data["raw_count"].eq(0).mean()
            percentage_histogram(
                axis,
                identity_data["log_count"].to_numpy(),
                bins,
                IDENTITY_COLORS[identity],
                label=(
                    f"{identity} (n={len(identity_data):,}; "
                    f"zeros={zero_fraction:.1f}%)"
                ),
                histtype="step",
                alpha=1.0,
                linewidth=1.6,
            )

        axis.set_title(f"{gene} in {pooled_label}", pad=7)
        axis.set_xlabel("log(1 + raw UMI count)")
        axis.set_ylabel("Cells within identity")
        axis.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        axis.legend(frameon=False, loc="upper right", handlelength=1.8)
        sns.despine(fig=fig)
        fig.tight_layout()
        save(fig, f"fig03{panel}_{dataset.lower()}_{gene}_marker_histogram.svg")


def urx_matched_histograms(data: pd.DataFrame, manifest: pd.DataFrame) -> None:
    rows = manifest[manifest["dataset"].eq("URX")].sort_values("panel_order")
    urx = data[data["dataset"].eq("URX")]
    bins = log_bins(urx["log_count"])

    fig, axes = plt.subplots(1, 3, figsize=(8.2, 2.85), sharex=True, sharey=True)
    for axis, row in zip(axes, rows.itertuples(index=False)):
        subset = urx[urx["panel_order"].eq(row.panel_order)]
        color = HVG_COLOR if row.expected_call == "HVG" else LVG_COLOR
        percentage_histogram(axis, subset["log_count"].to_numpy(), bins, color)
        axis.set_title(f"{row.gene_name} ({row.expected_call})", pad=7)
        axis.text(
            0.98, 0.96, basics_label(row),
            transform=axis.transAxes, ha="right", va="top", fontsize=8,
        )
        axis.set_xlabel("log(1 + raw UMI count)")
        axis.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    axes[0].set_ylabel("Cells")
    fig.suptitle("Matched expression distributions in URX", fontsize=11, y=1.02)
    sns.despine(fig=fig)
    fig.tight_layout()
    save(fig, "fig05A_urx_matched_hvg_lvg_histograms.svg")


def high_mean_low_zero_histograms(data: pd.DataFrame, manifest: pd.DataFrame) -> None:
    rows = manifest[manifest["dataset"].isin(["AVK", "RIG"])].sort_values("panel_order")
    selected = data[data["dataset"].isin(["AVK", "RIG"])]
    bins = log_bins(selected["log_count"])

    fig, axes = plt.subplots(1, 2, figsize=(5.7, 2.85), sharex=True, sharey=True)
    for axis, row in zip(axes, rows.itertuples(index=False)):
        subset = selected[selected["panel_order"].eq(row.panel_order)]
        percentage_histogram(axis, subset["log_count"].to_numpy(), bins, HVG_COLOR)
        axis.set_title(f"{row.gene_name} in {row.dataset} (HVG)", pad=7)
        axis.text(
            0.98, 0.96, basics_label(row),
            transform=axis.transAxes, ha="right", va="top", fontsize=8,
        )
        axis.set_xlabel("log(1 + raw UMI count)")
        axis.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    axes[0].set_ylabel("Cells")
    fig.suptitle("High-mean, low-zero HVGs", fontsize=11, y=1.02)
    sns.despine(fig=fig)
    fig.tight_layout()
    save(fig, "fig05B_high_mean_low_zero_hvg_histograms.svg")


def main() -> None:
    data = pd.read_csv(data_file("expression_histogram_data.csv.gz"))
    manifest = pd.read_csv(data_file("expression_histogram_manifest.csv"))

    set_style({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })
    marker_histograms(data)
    urx_matched_histograms(data, manifest)
    high_mean_low_zero_histograms(data, manifest)


if __name__ == "__main__":
    main()
