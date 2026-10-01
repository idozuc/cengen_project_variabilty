#!/usr/bin/env python3
"""Build versioned Streamlit assets from current BASiCS and WormCat outputs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd
import yaml

try:
    from .domain import CELL_TYPE_TO_FAMILY, RECEPTOR_TO_LIGANDS
except ImportError:  # Direct script execution.
    from domain import CELL_TYPE_TO_FAMILY, RECEPTOR_TO_LIGANDS

SCHEMA_VERSION = 2
MASTER_COLUMNS = {
    "cell_type", "stable_id", "gene_name", "Mu", "Delta", "Epsilon", "Prob",
    "HVG", "LVG", "n_cells", "n_experiments", "stress_gene", "analysis_eligible",
}
WORMCAT_COLUMNS = {
    "condition", "cell_type", "term_id", "term", "ontology", "p_value", "q_value",
    "significant", "term_selected_genes",
}


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--config", type=Path, help="Optional YAML with master, wormcat, and output_dir")
    result.add_argument("--master", type=Path, help="Current BASiCS master CSV or CSV.GZ")
    result.add_argument("--wormcat", type=Path, help="Current supported WormCat CSV or CSV.GZ")
    result.add_argument("--output-dir", type=Path, help="Explorer output directory")
    return result


def parse_args() -> argparse.Namespace:
    args = parser().parse_args()
    configured = {}
    if args.config:
        configured = yaml.safe_load(args.config.read_text()) or {}
    for key in ("master", "wormcat", "output_dir"):
        if getattr(args, key) is None and configured.get(key):
            setattr(args, key, Path(configured[key]))
        if getattr(args, key) is None:
            parser().error(f"--{key.replace('_', '-')} is required")
    return args


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} lacks required columns: {', '.join(missing)}")


def validate_wormcat(frame: pd.DataFrame, cell_types: set[str]) -> None:
    require_columns(frame, WORMCAT_COLUMNS, "WormCat table")
    if frame.condition.isna().any() or not frame.condition.eq("supported").all():
        raise ValueError("WormCat input must contain only the supported condition")
    if frame[["cell_type", "term_id", "term", "ontology"]].isna().any().any():
        raise ValueError("WormCat input contains missing category identifiers")
    if frame.duplicated(["cell_type", "term_id"]).any():
        raise ValueError("WormCat input contains duplicate cell-type/category keys")
    if not set(frame.cell_type).issubset(cell_types):
        raise ValueError("WormCat cell types do not match the BASiCS master")
    if not frame.ontology.isin(["category_1", "category_2", "category_3"]).all():
        raise ValueError("WormCat category levels must be category_1, category_2, or category_3")
    expected = frame.ontology + "::" + frame.term
    if not frame.term_id.eq(expected).all():
        raise ValueError("WormCat term IDs do not match their category level and name")
    if not frame.significant.isin([True, False]).all():
        raise ValueError("WormCat significant must contain boolean values")


def nuisance_class(row: pd.Series) -> str:
    gene = str(row.gene_name)
    patterns = (
        ("mitochondrial", r"^(nduo|ctb|ctc|cox|atp|cyc)-"),
        ("ribosomal", r"^(rps|rpl|mrps|mrpl)-"),
        ("histone", r"^his-"),
        ("cytoskeletal", r"^(act|tba|tbb|mua|pat|mup|ifb|ifc|ife|ifp)-"),
        ("transgene", r"^(unc-119|rol-6|lin-15[AB]|pha-1|myo-2|myo-3)$"),
    )
    if bool(row.stress_gene):
        return "stress"
    for label, pattern in patterns:
        if re.search(pattern, gene):
            return label
    return "biological"


def build_summary(master: pd.DataFrame) -> pd.DataFrame:
    summary = master.groupby("cell_type", as_index=False).agg(
        n_cells=("n_cells", "first"),
        n_experiments=("n_experiments", "first"),
        n_genes_tested=("stable_id", "size"),
        n_hvg=("HVG", "sum"),
    )
    biological = master[(master.HVG) & (master.analysis_eligible) & (master.nuisance == "biological")]
    counts = biological.groupby("cell_type").size().rename("n_hvg_biological")
    summary = summary.join(counts, on="cell_type")
    summary["n_hvg_biological"] = summary.n_hvg_biological.fillna(0).astype(int)
    return summary


def build_signatures(master: pd.DataFrame) -> pd.DataFrame:
    work = master[(master.HVG) & (master.analysis_eligible) & master.gene_name.str.match(r"^(flp|nlp|ins|pdf|ntr)-", na=False)].copy()
    work["family"] = work.cell_type.map(CELL_TYPE_TO_FAMILY).fillna("other")
    family_size = master.assign(family=master.cell_type.map(CELL_TYPE_TO_FAMILY).fillna("other")).groupby("family").cell_type.nunique()
    result = work.groupby(["family", "gene_name"]).cell_type.nunique().rename("n_hvg_types").reset_index()
    result["n_types_in_family"] = result.family.map(family_size)
    result["frac_family"] = result.n_hvg_types / result.n_types_in_family
    return result[result.frac_family >= 0.5].sort_values(["family", "frac_family"], ascending=[True, False])


def gpcr_class(gene: str) -> str | None:
    if re.match(r"^(npr|ckr|frpr|pdfr|nmur|lat)-", str(gene)):
        return "neuropeptide_receptor"
    if re.match(r"^(dop|ser|tyra|octr|mod)-", str(gene)):
        return "monoamine_receptor"
    return None


def build_ligand_receptor(master: pd.DataFrame) -> pd.DataFrame:
    eligible = master[(master.HVG) & (master.analysis_eligible)].copy()
    peptides = eligible[eligible.gene_name.str.match(r"^(flp|nlp|ins|pdf|ntr)-", na=False)]
    receptors = eligible[eligible.gene_name.map(gpcr_class).notna()]
    rows = []
    for receptor in receptors.itertuples(index=False):
        ligands = RECEPTOR_TO_LIGANDS.get(receptor.gene_name, [])
        matches = peptides[peptides.gene_name.isin(ligands)]
        if matches.empty:
            rows.append({
                "receptor": receptor.gene_name, "receptor_ct": receptor.cell_type,
                "gpcr_class": gpcr_class(receptor.gene_name), "receptor_epsilon": receptor.Epsilon,
                "ligand": ",".join(ligands) if ligands else None, "ligand_ct": None,
                "ligand_epsilon": None, "same_ct": None, "ligand_mapped": bool(ligands),
            })
        else:
            for ligand in matches.itertuples(index=False):
                rows.append({
                    "receptor": receptor.gene_name, "receptor_ct": receptor.cell_type,
                    "gpcr_class": gpcr_class(receptor.gene_name), "receptor_epsilon": receptor.Epsilon,
                    "ligand": ligand.gene_name, "ligand_ct": ligand.cell_type,
                    "ligand_epsilon": ligand.Epsilon, "same_ct": ligand.cell_type == receptor.cell_type,
                    "ligand_mapped": True,
                })
    return pd.DataFrame(rows, columns=[
        "receptor", "receptor_ct", "gpcr_class", "receptor_epsilon", "ligand",
        "ligand_ct", "ligand_epsilon", "same_ct", "ligand_mapped",
    ])


def main() -> None:
    args = parse_args()
    master = pd.read_csv(args.master)
    wormcat = pd.read_csv(args.wormcat)
    require_columns(master, MASTER_COLUMNS, "BASiCS master")
    key = master.cell_type.astype(str) + "::" + master.stable_id.astype(str)
    if key.duplicated().any():
        raise ValueError("BASiCS master contains duplicate cell-type/gene keys")
    validate_wormcat(wormcat, set(master.cell_type))

    master["nuisance"] = master.apply(nuisance_class, axis=1)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    master.to_parquet(output / "basics.parquet", index=False)
    build_summary(master).to_parquet(output / "basics_summary.parquet", index=False)
    wormcat.to_parquet(output / "wormcat_supported.parquet", index=False)
    build_signatures(master).to_csv(output / "neuropeptide_family_signatures.csv", index=False)
    build_ligand_receptor(master).to_csv(output / "gpcr_ligand_receptor_pairs.csv", index=False)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "master_source": str(args.master.resolve()),
        "wormcat_source": str(args.wormcat.resolve()),
        "files": [
            "basics.parquet", "basics_summary.parquet", "wormcat_supported.parquet",
            "neuropeptide_family_signatures.csv", "gpcr_ligand_receptor_pairs.csv",
        ],
    }
    (output / "schema_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
