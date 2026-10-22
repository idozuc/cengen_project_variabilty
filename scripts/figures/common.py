"""Shared configuration, plot style and helpers for scripts/figures/fig*.py.

Every script takes --config with a figures YAML (config/figures.example.yaml).
Paths in it are relative to the repository root, where the scripts are run.
"""

from __future__ import annotations

import argparse
import sys
from functools import cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml


SEED = 20260809

CLUSTER_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#000000",
]


@cache
def config() -> dict[str, Path]:
    parser = argparse.ArgumentParser(
        description=sys.modules["__main__"].__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path, required=True, help="Figures YAML")
    values = yaml.safe_load(parser.parse_args().config.read_text())
    output_dir = Path(values["output_dir"])
    return {
        "basics_dir": Path(values["basics_dir"]),
        "data_dir": output_dir / "data",  # tables written by the R exports
        "panel_dir": output_dir / "panels",
    }


def existing(path: Path, hint: str = "") -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing input: {path}{hint}")
    return path


def basics_file(*parts: str) -> Path:
    return existing(config()["basics_dir"].joinpath(*parts))


def data_file(name: str) -> Path:
    return existing(
        config()["data_dir"] / name,
        ". Run the scripts/figures/0*_export_*.R scripts first.",
    )


def set_style(rc: dict, theme: str = "ticks") -> None:
    sns.set_theme(style=theme, context="paper")
    plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", **rc})


def save(fig: plt.Figure, name: str, **savefig_kwargs) -> Path:
    panel_dir = config()["panel_dir"]
    panel_dir.mkdir(parents=True, exist_ok=True)
    path = panel_dir / name
    fig.savefig(path, bbox_inches="tight", **savefig_kwargs)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


def ordered_labels(values: pd.Series) -> list[str]:
    """Unique labels, numeric ones first in numeric order."""
    labels = values.dropna().astype(str).unique().tolist()

    def key(value: str) -> tuple[int, float | str]:
        try:
            return (0, float(value))
        except ValueError:
            return (1, value)

    return sorted(labels, key=key)


def cluster_palette(values: pd.Series) -> dict[str, str]:
    labels = ordered_labels(values)
    if len(labels) > len(CLUSTER_COLORS):
        colors = sns.color_palette("colorblind", len(labels)).as_hex()
    else:
        colors = CLUSTER_COLORS[: len(labels)]
    return dict(zip(labels, colors))


def pc_columns(frame: pd.DataFrame) -> list[str]:
    columns = [
        column
        for column in frame.columns
        if column.startswith("PC") and column[2:].isdigit()
    ]
    return sorted(columns, key=lambda column: int(column[2:]))


def calculate_umap(frame: pd.DataFrame) -> np.ndarray:
    """2-D UMAP of the PC scores. Seeded and single-threaded, so reproducible."""
    import umap

    values = frame[pc_columns(frame)].to_numpy(dtype=float)
    if values.shape[1] < 2 or not np.isfinite(values).all():
        raise ValueError("UMAP needs at least two finite PC columns")

    reducer = umap.UMAP(
        n_neighbors=min(15, len(frame) - 1),
        min_dist=0.30,
        metric="euclidean",
        n_components=2,
        random_state=SEED,
        transform_seed=SEED,
        n_jobs=1,
    )
    return reducer.fit_transform(values)
