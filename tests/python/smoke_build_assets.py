#!/usr/bin/env python3
"""Dependency-light executable smoke test for explorer asset generation."""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    requested = os.environ.get("CENGEN_SMOKE_OUTPUT")
    context = tempfile.TemporaryDirectory() if not requested else None
    try:
        work = Path(requested) if requested else Path(context.name)
        work.mkdir(parents=True, exist_ok=True)
        master = pd.DataFrame({
            "cell_type": ["SMD", "SMD", "RIA", "RIA"],
            "stable_id": ["WBGene1", "WBGene2", "WBGene1", "WBGene2"],
            "gene_name": ["flp-1", "npr-22", "flp-1", "npr-22"],
            "Mu": [1.0] * 4, "Delta": [1.0] * 4,
            "Epsilon": [2.0, 1.5, 1.8, 1.2], "Prob": [0.99, 0.98, 0.97, 0.96],
            "HVG": [True] * 4, "LVG": [False] * 4,
            "n_cells": [30, 30, 25, 25], "n_experiments": [3] * 4,
            "stress_gene": [False] * 4, "analysis_eligible": [True] * 4,
        })
        wormcat = pd.DataFrame({
            "condition": ["supported"], "cell_type": ["SMD"], "term_id": ["category_1::test process"],
            "term": ["test process"], "ontology": ["category_1"],
            "p_value": [0.001], "q_value": [0.01], "significant": [True],
            "term_selected_genes": [3],
        })
        master_path, wormcat_path, output = work / "master.csv.gz", work / "wormcat.csv.gz", work / "assets"
        master.to_csv(master_path, index=False)
        wormcat.to_csv(wormcat_path, index=False)
        subprocess.run([
            sys.executable, str(root / "explorer" / "build_assets.py"),
            "--master", str(master_path), "--wormcat", str(wormcat_path), "--output-dir", str(output),
        ], check=True)
        manifest = json.loads((output / "schema_manifest.json").read_text())
        assert manifest["schema_version"] == 2
        assert len(pd.read_parquet(output / "basics.parquet")) == 4
    finally:
        if context is not None:
            context.cleanup()
    print("Explorer asset smoke test passed")


if __name__ == "__main__":
    main()
