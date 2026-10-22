#!/usr/bin/env python3
"""WormCat HVG enrichment restricted to genes CeNGEN calls expressed.

Guards against ambient RNA: a gene/cell-type pair enters neither the HVG list
nor the background unless the CeNGEN call table (e.g. threshold 2) calls the
gene expressed in that cell type. Cell types without CeNGEN calls are dropped;
genes absent from the call table are kept (no call either way).

The test follows the WormCat tool: one-sided hypergeometric test of annotated
HVGs against the annotated background, Benjamini-Hochberg within each cell
type and category level over the categories that contain at least one HVG.
Fold enrichment uses all HVGs and all background genes. A row is significant
at FDR <= 0.05 and fold enrichment > 1.

HVGs are BASiCS HVG calls among analysis-eligible pairs that are also detected
in more than half of the cell type's experiments.

Output: the rows with FDR <= 0.05, one per cell type x level x category.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


COLUMNS = [
    "cell_type", "category_level", "wormcat_family", "RGS", "AC",
    "n_hvg", "n_background", "fold_enrichment", "PValue", "FDR",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--master", type=Path, required=True,
                        help="basics_master_with_detection.csv.gz (scripts/analysis/04_add_detection_columns.R)")
    parser.add_argument("--wormcat", type=Path, required=True, help="WormCat whole-genome CSV")
    parser.add_argument("--cengen-calls", type=Path, required=True,
                        help="CeNGEN gene x cell-type call table, e.g. 021821_medium_threshold2.csv")
    parser.add_argument("--output", type=Path, required=True, help="Output CSV")
    return parser.parse_args()


def bh(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted P values, as R's p.adjust(method = "BH")."""
    order = np.argsort(p.values)[::-1]
    q = np.minimum.accumulate(p.values[order] * len(p) / np.arange(len(p), 0, -1))
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(q, 1)
    return pd.Series(adjusted, index=p.index)


def load_pairs(master_path: Path, wormcat_path: Path, calls_path: Path) -> pd.DataFrame:
    """Eligible gene/cell-type pairs with final HVG flag, WormCat categories and CeNGEN call."""
    m = pd.read_csv(master_path)
    m = m[m.analysis_eligible].copy()
    detected = m.experiment_detection.str.split("/", expand=True).astype(int)
    m["final_hvg"] = m.HVG & (detected[0] > detected[1] / 2)

    wormcat = pd.read_csv(wormcat_path, encoding="utf-8-sig")
    wormcat = wormcat.drop_duplicates("Wormbase ID").set_index("Wormbase ID")
    for level in (1, 2, 3):
        m[f"cat{level}"] = m.stable_id.map(wormcat[f"Category {level}"])

    calls = pd.read_csv(calls_path).drop(columns=["Unnamed: 0", "gene_name"], errors="ignore")
    calls = calls.melt(id_vars="Wormbase_ID", var_name="cell_type", value_name="tpm")
    calls["expressed"] = calls.tpm > 0
    return m.merge(
        calls[["Wormbase_ID", "cell_type", "expressed"]],
        how="left",
        left_on=["stable_id", "cell_type"],
        right_on=["Wormbase_ID", "cell_type"],
    ).drop(columns="Wormbase_ID")


def cengen_filter(pairs: pd.DataFrame) -> pd.DataFrame:
    covered = pairs.groupby("cell_type").expressed.apply(lambda s: s.notna().any())
    keep = covered.reindex(pairs.cell_type).values & (pairs.expressed != False)  # noqa: E712 (NaN = no call, kept)
    return pairs[keep]


def wormcat_enrichment(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell_type, group in pairs.groupby("cell_type"):
        n_all, N_all = int(group.final_hvg.sum()), len(group)
        for level in (1, 2, 3):
            annotated = group.dropna(subset=[f"cat{level}"])
            n, N = int(annotated.final_hvg.sum()), len(annotated)
            hvg_counts = annotated[annotated.final_hvg][f"cat{level}"].value_counts()
            background_counts = annotated[f"cat{level}"].value_counts()
            for category, k in hvg_counts.items():
                K = background_counts[category]
                rows.append(dict(
                    cell_type=cell_type, category_level=level, wormcat_family=category,
                    RGS=k, AC=K, n_hvg=n_all, n_background=N_all,
                    fold_enrichment=(k / n_all) / (K / N_all),
                    PValue=hypergeom.sf(k - 1, N, K, n),
                ))
    result = pd.DataFrame(rows)
    result["FDR"] = result.groupby(["cell_type", "category_level"]).PValue.transform(bh)
    return result


def main() -> None:
    args = parse_args()
    filtered = cengen_filter(load_pairs(args.master, args.wormcat, args.cengen_calls))
    result = wormcat_enrichment(filtered)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.loc[result.FDR <= 0.05, COLUMNS].to_csv(args.output, index=False)

    hvgs = filtered[filtered.final_hvg]
    significant = result[(result.FDR <= 0.05) & (result.fold_enrichment > 1)]
    print(f"{len(hvgs)} final HVG pairs in {filtered.cell_type.nunique()} cell types; "
          f"{len(significant)} significant enrichments")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
