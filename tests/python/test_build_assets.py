import json
from pathlib import Path

import pandas as pd
import pytest

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
    wormcat = pd.DataFrame(
        {
            "condition": ["supported"], "cell_type": ["SMD"],
            "term_id": ["category_1::test process"], "term": ["test process"],
            "ontology": ["category_1"], "p_value": [0.001],
            "q_value": [0.01], "significant": [True], "term_selected_genes": [3],
        }
    )
    master_path = tmp_path / "master.csv.gz"
    wormcat_path = tmp_path / "wormcat.csv.gz"
    output = tmp_path / "assets"
    master.to_csv(master_path, index=False)
    wormcat.to_csv(wormcat_path, index=False)
    monkeypatch.setattr(
        "sys.argv",
        ["build_assets.py", "--master", str(master_path), "--wormcat", str(wormcat_path), "--output-dir", str(output)],
    )
    build_assets.main()

    manifest = json.loads((output / "schema_manifest.json").read_text())
    assert manifest["schema_version"] == 2
    assert set(manifest["files"]) == {
        "basics.parquet", "basics_summary.parquet", "wormcat_supported.parquet",
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


def wormcat_fixture():
    return pd.DataFrame({
        "condition": ["supported"] * 3, "cell_type": ["SMD"] * 3,
        "term_id": [f"category_{i}::test" for i in (1, 2, 3)],
        "term": ["test"] * 3,
        "ontology": [f"category_{i}" for i in (1, 2, 3)],
        "p_value": [0.001] * 3, "q_value": [0.01] * 3,
        "significant": [True, False, True], "term_selected_genes": [3] * 3,
    })


def test_wormcat_levels_and_empty_results():
    frame = wormcat_fixture()
    build_assets.validate_wormcat(frame, {"SMD"})
    build_assets.validate_wormcat(frame.iloc[:0], {"SMD"})


@pytest.mark.parametrize("failure", ["missing", "duplicate", "cell_type", "condition", "level", "id"])
def test_invalid_wormcat_is_rejected(failure):
    frame = wormcat_fixture()
    if failure == "missing":
        frame = frame.drop(columns="q_value")
    elif failure == "duplicate":
        frame = pd.concat([frame, frame.iloc[:1]])
    elif failure == "cell_type":
        frame.loc[0, "cell_type"] = "UNKNOWN"
    elif failure == "condition":
        frame.loc[0, "condition"] = "full"
    elif failure == "level":
        frame.loc[0, "ontology"] = "biological_process"
    else:
        frame.loc[0, "term_id"] = "wrong"
    with pytest.raises(ValueError):
        build_assets.validate_wormcat(frame, {"SMD"})
