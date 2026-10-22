#!/usr/bin/env python3
"""Figure 4: canonical-marker evidence within the known ASE and AWC identities.

Each cell is classed by which tested markers it expresses (>= 1 UMI):
the markers expected for its identity, the opposing identity's markers,
both, or none.

Input: OUTPUT_DIR/data/expression_histogram_data.csv.gz (01_export_expression_data.R).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter

from common import data_file, save, set_style


CATEGORY_ORDER = [
    "Expected marker(s)",
    "Both marker groups",
    "Opposing marker(s) only",
    "No tested marker",
]

CATEGORY_COLORS = {
    "Expected marker(s)": "#009E73",
    "Both marker groups": "#CC79A7",
    "Opposing marker(s) only": "#D55E00",
    "No tested marker": "#BDBDBD",
}

IDENTITY_COLORS = {
    "ASEL": "#0072B2",
    "ASER": "#E69F00",
    "AWC_ON": "#009E73",
    "AWC_OFF": "#CC79A7",
}

PANELS = [
    ("A", "AWC", {"AWC_ON": ["str-2"], "AWC_OFF": ["srsx-3"]}),
    ("B", "ASE", {"ASEL": ["gcy-6", "gcy-7"], "ASER": ["gcy-5", "gcy-22"]}),
]


def marker_matrix(data: pd.DataFrame, dataset: str, markers: list[str]) -> pd.DataFrame:
    subset = data[data["dataset"].eq(dataset) & data["gene_name"].isin(markers)]
    wide = subset.pivot(
        index=["cell", "identity"],
        columns="gene_name",
        values="raw_count",
    ).reset_index()

    missing = sorted(set(markers) - set(wide.columns))
    if missing:
        raise ValueError(f"Missing {dataset} markers: {', '.join(missing)}")
    if wide[markers].isna().any().any():
        raise ValueError(f"Incomplete marker-by-cell matrix for {dataset}.")
    return wide


def classify_evidence(wide: pd.DataFrame, marker_groups: dict[str, list[str]]) -> pd.DataFrame:
    records = []
    for identity, expected_markers in marker_groups.items():
        identity_data = wide[wide["identity"].eq(identity)]
        opposing_markers = [
            marker
            for other, markers in marker_groups.items()
            if other != identity
            for marker in markers
        ]
        expected = identity_data[expected_markers].gt(0).any(axis=1)
        opposing = identity_data[opposing_markers].gt(0).any(axis=1)

        category = np.select(
            [expected & ~opposing, expected & opposing, ~expected & opposing],
            CATEGORY_ORDER[:3],
            default="No tested marker",
        )
        counts = pd.Series(category).value_counts()
        total = len(identity_data)
        for status in CATEGORY_ORDER:
            count = int(counts.get(status, 0))
            records.append({
                "identity": identity,
                "status": status,
                "count": count,
                "fraction": count / total,
            })
    return pd.DataFrame.from_records(records)


def add_segment_label(axis: plt.Axes, left: float, width: float, y: float, count: int) -> None:
    if count == 0:
        return
    label = f"{count:,}\n({width:.1%})"
    if width >= 0.075:
        axis.text(
            left + width / 2, y, label,
            ha="center", va="center", fontsize=8, color="#222222",
        )
    else:
        # Narrow segment: label above the bar with a leader line.
        axis.annotate(
            label.replace("\n", " "),
            xy=(left + width / 2, y + 0.27),
            xytext=(left + width / 2, y + 0.47),
            ha="center",
            va="bottom",
            fontsize=7.5,
            color="#333333",
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.6},
        )


def evidence_figure(data: pd.DataFrame, panel: str, dataset: str, marker_groups: dict) -> None:
    identities = list(marker_groups)
    markers = [marker for group in marker_groups.values() for marker in group]
    summary = classify_evidence(marker_matrix(data, dataset, markers), marker_groups)

    fig, axis = plt.subplots(figsize=(5.8, 3.25))
    y_positions = np.arange(len(identities))[::-1]

    for y, identity in zip(y_positions, identities):
        left = 0.0
        identity_summary = summary[summary["identity"].eq(identity)].set_index("status")
        for status in CATEGORY_ORDER:
            width = float(identity_summary.at[status, "fraction"])
            axis.barh(
                y, width, left=left, height=0.55,
                color=CATEGORY_COLORS[status], edgecolor="white", linewidth=0.8,
            )
            add_segment_label(axis, left, width, y, int(identity_summary.at[status, "count"]))
            left += width

    axis.set_yticks(y_positions, labels=identities)
    for tick, identity in zip(axis.get_yticklabels(), identities):
        tick.set_color(IDENTITY_COLORS[identity])
        tick.set_fontweight("bold")

    axis.set_xlim(0, 1)
    axis.set_ylim(-0.65, len(identities) - 0.25)
    axis.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    axis.set_xlabel("Cells within assigned identity")
    axis.set_title(f"Canonical marker evidence among assigned {dataset} identities", pad=12)

    axis.legend(
        handles=[Patch(facecolor=CATEGORY_COLORS[s], label=s) for s in CATEGORY_ORDER],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.26),
        ncol=2,
        columnspacing=1.2,
        handlelength=1.2,
    )
    marker_description = "   |   ".join(
        f"{identity}: {', '.join(markers)}" for identity, markers in marker_groups.items()
    )
    fig.text(
        0.5, 0.015, f"Expected markers — {marker_description}",
        ha="center", va="bottom", fontsize=7.5, color="#555555",
    )

    sns.despine(ax=axis, left=True)
    axis.tick_params(axis="y", length=0)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    save(fig, f"fig04{panel}_{dataset.lower()}_identity_marker_evidence.svg")


def main() -> None:
    data = pd.read_csv(data_file("expression_histogram_data.csv.gz"))
    set_style({
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 9,
    })
    for panel, dataset, marker_groups in PANELS:
        evidence_figure(data, panel, dataset, marker_groups)


if __name__ == "__main__":
    main()
