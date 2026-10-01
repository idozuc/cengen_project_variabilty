import json
from pathlib import Path

import pandas as pd

from explorer import build_assets


def test_build_assets_from_current_schema(tmp_path: Path, monkeypatch):
    master = pd.DataFrame(
        {
            "cell_type": ["SMD", "SMD", "RIA", "RIA"],
            "stable_id": ["WBGene1", "WBGene2", "WBGene1", "WBGene2"],
            "gene_name": ["flp-1", "npr-22", "flp-1", "npr-22"],
            "Mu": [1.0, 1.0, 1.0, 1.0],
            "Delta": [1.0, 1.0, 1.0, 1.0],
            "Epsilon": [2.0, 1.5, 1.8, 1.2],
            "Prob": [0.99, 0.98, 0.97, 0.96],
            "HVG": [True, True, True, True],
            "LVG": [False, False, False, False],
            "n_cells": [30, 30, 25, 25],
            "n_experiments": [3, 3, 3, 3],
            "stress_gene": [False, False, False, False],
            "analysis_eligible": [True, True, True, True],
        }
    )
    go = pd.DataFrame(
        {
            "condition": ["supported"], "cell_type": ["SMD"],
            "term_id": ["GO:1"], "term": ["test process"],
            "ontology": ["biological_process"], "p_value": [0.001],
            "q_value": [0.01], "significant": [True], "term_selected_genes": [3],
        }
    )
    master_path = tmp_path / "master.csv.gz"
    go_path = tmp_path / "go.csv.gz"
    output = tmp_path / "assets"
    master.to_csv(master_path, index=False)
    go.to_csv(go_path, index=False)
    monkeypatch.setattr(
        "sys.argv",
        ["build_assets.py", "--master", str(master_path), "--go", str(go_path), "--output-dir", str(output)],
    )
    build_assets.main()

    manifest = json.loads((output / "schema_manifest.json").read_text())
    assert manifest["schema_version"] == 1
    assert set(manifest["files"]) == {
        "basics.parquet", "basics_summary.parquet", "go_supported.parquet",
        "neuropeptide_family_signatures.csv", "gpcr_ligand_receptor_pairs.csv",
    }
    result = pd.read_parquet(output / "basics.parquet")
    assert not result.duplicated(["cell_type", "stable_id"]).any()
    assert set(result["nuisance"]) == {"biological"}


def test_duplicate_master_keys_are_rejected(tmp_path: Path):
    frame = pd.DataFrame({column: [None] for column in build_assets.MASTER_COLUMNS})
    duplicate = pd.concat([frame, frame], ignore_index=True)
    key = duplicate.cell_type.astype(str) + "::" + duplicate.stable_id.astype(str)
    assert key.duplicated().any()
