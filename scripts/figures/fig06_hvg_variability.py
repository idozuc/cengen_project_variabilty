#!/usr/bin/env python3
"""Figure 6: BASiCS residual variability and detection frequency.

A: distribution of HVG calls across 5-point zero-expression bins.
B: fraction of eligible gene/cell-type pairs called HVG in each bin.
C: mean expression (Mu) versus residual variability (Epsilon).

Input: BASICS_DIR/basics_master_with_detection.csv.gz
       (one row per gene/cell-type pair; scripts/analysis/04_add_detection_columns.R).
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import PercentFormatter

from common import basics_file, save, set_style


HVG_COLOR = "#D55E00"
LVG_COLOR = "#0072B2"
BIN_EDGES = np.arange(0, 105, 5)


def load_master() -> pd.DataFrame:
    df = pd.read_csv(
        basics_file("basics_master_with_detection.csv.gz")
    )
    df["experiment_detection_fraction"] = df["detected_experiments"] / df["n_experiments"]

    # Analysis HVGs must also be detected in more than half of the experiments.
    df["analysis_lvg"] = df["LVG"] & df["analysis_eligible"]
    df["analysis_hvg"] = (
        df["HVG"]
        & df["analysis_eligible"]
        & (df["experiment_detection_fraction"] > 0.5)
    )
    df["analysis_other"] = ~(df["analysis_lvg"] | df["analysis_hvg"])
    assert not (df["analysis_lvg"] & df["analysis_hvg"]).any()
    print(df[["analysis_other", "analysis_lvg", "analysis_hvg"]].sum())
    return df


def binned(frame: pd.DataFrame) -> pd.Series:
    return pd.cut(frame["pct_zero"], bins=BIN_EDGES, include_lowest=True)


def line_plot(x: pd.Series, y: pd.Series, ylabel: str, formatter_max: float, name: str) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    ax.plot(x, y, color=HVG_COLOR, linewidth=1.5, marker="o", markersize=4)
    ax.set_xlabel("Cells with zero expression (%)")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(PercentFormatter(formatter_max))
    sns.despine()
    fig.tight_layout()
    save(fig, name)


def hvg_distribution(df: pd.DataFrame) -> None:
    hvgs = df.loc[df["analysis_hvg"], ["pct_zero"]].copy()
    hvgs["pct_zero_bin"] = binned(hvgs)
    counts = hvgs.groupby("pct_zero_bin", observed=True).size()
    midpoints = counts.index.map(lambda interval: interval.mid).astype(float)
    line_plot(
        midpoints,
        100 * counts / counts.sum(),
        "HVGs in bin (%)",
        100,
        "fig06A_hvg_distribution_by_pct_zero.svg",
    )


def hvg_fraction(df: pd.DataFrame) -> None:
    eligible = df.loc[df["analysis_eligible"], ["pct_zero", "analysis_hvg"]].copy()
    eligible["pct_zero_bin"] = binned(eligible)
    fraction = eligible.groupby("pct_zero_bin", observed=True)["analysis_hvg"].mean()
    midpoints = fraction.index.map(lambda interval: interval.mid).astype(float)
    line_plot(
        midpoints,
        fraction,
        "Fraction of genes classified as HVG",
        1,
        "fig06B_hvg_fraction_vs_pct_zero.svg",
    )


def mu_vs_epsilon(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    groups = [
        ("analysis_other", "Other genes", "#BDBDBD", 2, 0.20),
        ("analysis_lvg", "LVG", LVG_COLOR, 5, 0.65),
        ("analysis_hvg", "HVG", HVG_COLOR, 5, 0.65),
    ]
    for zorder, (column, label, color, size, alpha) in enumerate(groups, start=1):
        mask = df[column]
        ax.scatter(
            df.loc[mask, "Mu"],
            df.loc[mask, "Epsilon"],
            s=size,
            color=color,
            alpha=alpha,
            edgecolors="none",
            label=f"{label} (n={mask.sum():,})",
            zorder=zorder,
            rasterized=True,
        )

    ax.set_xscale("log")  # Mu spans several orders of magnitude
    ax.axhline(0, color="#444444", linewidth=0.6, linestyle="--")
    ax.set_xlabel("Mean expression (Mu)")
    ax.set_ylabel("Residual variability (Epsilon)")
    ax.legend(frameon=False, markerscale=2.5, loc="best")
    sns.despine()
    fig.tight_layout()
    save(fig, "fig06C_mu_vs_epsilon.svg")


def main() -> None:
    df = load_master()
    set_style({"font.size": 9})
    hvg_distribution(df)
    hvg_fraction(df)
    mu_vs_epsilon(df)


if __name__ == "__main__":
    main()
