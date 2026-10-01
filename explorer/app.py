"""
BASiCS within-type HVG explorer — Streamlit app.
Run: streamlit run explorer/app.py
"""

import json
import os
import re
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DIR = Path(os.environ.get("CENGEN_OUTPUT_DIR", PROJECT_ROOT / "outputs" / "explorer"))
CLUSTER_DIR = DIR / "clustering"
SCHEMA_VERSION = 2
CLUSTER_SCHEMA_VERSION = 1

st.set_page_config(page_title="BASiCS HVG Explorer", layout="wide")

# ── family map & gene-class regex (from scripts/12_explorer.py) ──────────────
FAMILIES = {
    "motor_body":  ["DA","DB","VA","VB","AS","VD_DD","VC","VC_4_5"],
    "motor_head":  ["SMB","SMD","SMD_stressed","RMD_DV","RMD_LR","RME_DV","RME_LR","RMF","RMG","RMH"],
    "sensory":     ["AFD","AWA","AWB","AWC_ON","AWC_OFF","ASEL","ASER","ASH","ASI","ASJ","ASK",
                    "ADF","ADL","AQR","BAG","FLP","PDE","URX","URY","PHA","PHB","PHC",
                    "OLQ","CEP","IL1","IL2_DV","IL2_LR"],
    "interneuron": ["RIA","RIB","RIC","RID","RIF","RIG","RIH","RIM","RIP","RIR","RIS","RIV",
                    "RIV_stressed","AIB","AIA","AIM","AIN","AIY","AVA","AVB","AVD","AVE","AVF",
                    "AVG","AVH","AVJ","AVK","AVL","DVA","DVB","DVC","LUA","PVC","PVD","PVM",
                    "PVP","PVQ","PVR","PVW","PLM","PLN","BDU","ADA","ADE","AUA","SAA","SAB",
                    "SDQ","URA","URB","AVM","ALM","PQR"],
    "pharyngeal":  ["I1","I4","I5","I6","M1","M2","M3","M5","MI","NSM","MC"],
    "other":       ["HSN","SIA","SIB"],
}
CT_TO_FAMILY = {ct: f for f, cts in FAMILIES.items() for ct in cts}

GENE_CLASS_PATTERNS = [
    ("neuropeptide",          r"^(flp|nlp|ins|pdf|ntr)-"),
    ("neuropeptide_receptor", r"^(npr|ckr|frpr|pdfr|nmur|lat)-"),
    ("monoamine_receptor",    r"^(dop|ser|tyra|octr|mod)-"),
    ("ribosomal",             r"^(rps|rpl|mrps|mrpl)-"),
    ("mitochondrial",         r"^(nduo|ctb|ctc|cox|atp|cyc)-"),
    ("histone",               r"^his-"),
    ("stress",                r"^hsp-|^xbp-1$"),
    ("cytoskeletal",          r"^(act|tba|tbb|mua|pat|mup|ifb|ifc|ife|ifp)-"),
    ("transgene",             r"^(unc-119|rol-6|lin-15[AB]|pha-1|myo-2|myo-3)$"),
]
NUISANCE_CLASSES = {"ribosomal","mitochondrial","histone","stress","cytoskeletal","transgene"}

def classify_gene(gene):
    for cls, pat in GENE_CLASS_PATTERNS:
        if re.search(pat, str(gene)):
            return cls
    return "other"

FAMILY_COLORS = {
    "sensory": "#4C9BE8", "interneuron": "#E8954C", "motor_body": "#6BBF6B",
    "motor_head": "#A67BC8", "pharyngeal": "#E86B6B", "other": "#AAAAAA",
}

# ── data loader ──────────────────────────────────────────────────────────────
@st.cache_data
def load():
    manifest_path = DIR / "schema_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing explorer schema manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"Explorer schema {manifest.get('schema_version')} is unsupported; expected {SCHEMA_VERSION}"
        )
    df = pd.read_parquet(DIR / "basics.parquet")
    required = {
        "cell_type", "stable_id", "gene_name", "Mu", "Delta", "Epsilon", "Prob",
        "HVG", "n_cells", "n_experiments", "analysis_eligible", "nuisance",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"basics.parquet lacks required columns: {', '.join(missing)}")
    df["neuron_family"] = df["cell_type"].map(CT_TO_FAMILY).fillna("other")
    df["gene_class"] = df["gene_name"].apply(classify_gene)
    for c in ["Mu","Delta","Epsilon","Prob"]:
        df[c] = df[c].round(3)
    summary = pd.read_parquet(DIR / "basics_summary.parquet")
    wormcat_df = pd.read_parquet(DIR / "wormcat_supported.parquet")
    try:
        from explorer.build_assets import validate_wormcat
    except ImportError:
        from build_assets import validate_wormcat
    validate_wormcat(wormcat_df, set(df.cell_type))
    wormcat_df = wormcat_df[wormcat_df.significant].copy()
    lr = pd.read_csv(DIR / "gpcr_ligand_receptor_pairs.csv")
    sigs = pd.read_csv(DIR / "neuropeptide_family_signatures.csv")
    return df, summary, wormcat_df, lr, sigs

df, summary, wormcat_df, lr, sigs = load()

# ── nuisance-adjusted clustering assets ──────────────────────────────────────
@st.cache_data
def load_cluster_manifest():
    path = CLUSTER_DIR / "manifest.parquet"
    if not path.exists():
        return pd.DataFrame()
    manifest = pd.read_parquet(path)
    if "schema_version" not in manifest or not (manifest["schema_version"] == CLUSTER_SCHEMA_VERSION).all():
        raise ValueError(f"Unsupported clustering schema in {path}")
    return manifest.sort_values("cell_type").reset_index(drop=True)


@st.cache_data(max_entries=8)
def load_cluster_assets(cell_type):
    directory = CLUSTER_DIR / cell_type
    cells = pd.read_parquet(directory / "cells.parquet")
    cluster_summary = pd.read_parquet(directory / "summary.parquet")
    qc = pd.read_parquet(directory / "qc_by_resolution.parquet")
    markers = pd.read_parquet(directory / "markers.parquet")
    return cells, cluster_summary, qc, markers


@st.cache_data(max_entries=32)
def load_cluster_expression(cell_type, feature):
    path = CLUSTER_DIR / cell_type / "expression.parquet"
    return pd.read_parquet(path, columns=["cell_id", feature]).rename(
        columns={feature: "normalized_expression"}
    )


cluster_manifest = load_cluster_manifest()

# ── sidebar filters ──────────────────────────────────────────────────────────
st.sidebar.title("Filters")
st.sidebar.caption("Applied to Browse, Gene, Family tabs")

all_families = list(FAMILIES.keys())
all_classes  = sorted(df["gene_class"].unique())

sel_families = st.sidebar.multiselect("Neuron family", all_families, default=all_families)
sel_classes  = st.sidebar.multiselect("Gene class", all_classes, default=all_classes)
bio_only     = st.sidebar.toggle("Biological genes only (exclude nuisance)", value=True)
min_eps      = st.sidebar.slider("Min Epsilon", 0.0, 5.0, 0.0, 0.1)
min_prob     = st.sidebar.slider("Min Prob (HVG posterior)", 0.0, 1.0, 0.0, 0.05)

def apply_filters(d, use_hvg=True):
    m = d["neuron_family"].isin(sel_families) & d["gene_class"].isin(sel_classes)
    if use_hvg:
        m &= (d["HVG"] == True) & (d["analysis_eligible"] == True)
    if bio_only:
        m &= d["nuisance"] == "biological"
    m &= d["Epsilon"] >= min_eps
    m &= d["Prob"] >= min_prob
    return d[m]

filt = apply_filters(df)

# ── tabs ─────────────────────────────────────────────────────────────────────
st.title("BASiCS within-type variability — CeNGEN L4 neurons")
t1, t2, t3, t4, t5, t6, t7, t8 = st.tabs(
    ["📋 Browse", "🔲 Matrix", "🧬 Gene", "🔬 Cell type", "👪 Family",
     "🧪 WormCat enrichment", "🔗 Ligand–Receptor", "🫧 Clusters"]
)

# ── Tab 1: Browse ────────────────────────────────────────────────────────────
with t1:
    c1, c2, c3 = st.columns(3)
    c1.metric("Rows shown", f"{len(filt):,}")
    c2.metric("Unique genes", filt["gene_name"].nunique())
    c3.metric("Cell types", filt["cell_type"].nunique())
    cols = ["gene_name","gene_class","cell_type","neuron_family","Epsilon","Mu","Prob","n_cells","n_experiments"]
    st.dataframe(
        filt[cols].sort_values("Epsilon", ascending=False),
        width="stretch", height=600, hide_index=True,
    )
    st.download_button("Download filtered CSV",
                       filt[cols].to_csv(index=False).encode(),
                       "basics_filtered.csv", "text/csv")

# ── Tab 2: Matrix (gene × cell type) ─────────────────────────────────────────
with t2:
    st.caption("Genes (rows) × cell types (columns). Cell = Epsilon where the gene is HVG "
               "(blank/0 = not HVG). Filters below are independent of the sidebar.")

    mc1, mc2 = st.columns(2)
    with mc1:
        m_gclass = st.multiselect("Gene class", all_classes,
                                  default=["neuropeptide"] if "neuropeptide" in all_classes else all_classes[:1],
                                  key="mx_gclass")
        m_genes = st.multiselect("Specific genes (optional — overrides gene class)",
                                 sorted(df["gene_name"].unique()), key="mx_genes")
    with mc2:
        m_fam = st.multiselect("Neuron family", all_families, default=all_families, key="mx_fam")
        m_cts = st.multiselect("Specific cell types (optional — overrides family)",
                               sorted(df["cell_type"].unique()), key="mx_cts")

    m_val = st.radio("Cell value", ["Epsilon", "HVG (0/1)"], horizontal=True, key="mx_val")
    m_hvgonly = st.toggle("Only HVG cells (blank the rest)", value=True, key="mx_hvgonly")

    # row (gene) selection
    md = df.copy()
    if m_genes:
        md = md[md["gene_name"].isin(m_genes)]
    else:
        md = md[md["gene_class"].isin(m_gclass)]
    # column (cell type) selection
    if m_cts:
        md = md[md["cell_type"].isin(m_cts)]
    else:
        md = md[md["neuron_family"].isin(m_fam)]

    # keep only genes that are HVG somewhere in the current selection (avoids empty rows)
    hvg_genes = md.loc[md["HVG"] == True, "gene_name"].unique()
    md = md[md["gene_name"].isin(hvg_genes)]

    if md.empty:
        st.info("No HVG genes match the current selection.")
    else:
        value_col = "Epsilon" if m_val == "Epsilon" else "HVG"
        work = md.copy()
        if m_hvgonly:
            work.loc[work["HVG"] != True, "Epsilon"] = np.nan
        work["HVG"] = work["HVG"].astype(int)
        pivot = work.pivot_table(index="gene_name", columns="cell_type",
                                 values=value_col, aggfunc="first")
        # order rows by how many cell types they're HVG in, cols by n HVG
        row_order = (md[md["HVG"] == True].groupby("gene_name")["cell_type"].nunique()
                     .reindex(pivot.index).fillna(0).sort_values(ascending=False).index)
        col_order = (md[md["HVG"] == True].groupby("cell_type")["gene_name"].nunique()
                     .reindex(pivot.columns).fillna(0).sort_values(ascending=False).index)
        pivot = pivot.loc[row_order, col_order]

        st.write(f"**{pivot.shape[0]} genes × {pivot.shape[1]} cell types**")

        # heatmap (cap height/width to stay readable; scroll via plotly)
        fig = px.imshow(
            pivot, aspect="auto",
            color_continuous_scale="YlOrRd",
            labels=dict(x="Cell type", y="Gene", color=value_col),
        )
        fig.update_layout(
            height=max(400, min(pivot.shape[0]*16 + 150, 1400)),
            xaxis_tickangle=-45,
        )
        fig.update_xaxes(tickfont=dict(size=8))
        fig.update_yaxes(tickfont=dict(size=8))
        st.plotly_chart(fig, width="stretch")

        with st.expander("Show as table"):
            st.dataframe(pivot, width="stretch", height=500)
        st.download_button("Download matrix CSV", pivot.to_csv().encode(),
                           "basics_matrix.csv", "text/csv")

# ── Tab 3: Gene ──────────────────────────────────────────────────────────────
with t3:
    genes = sorted(df["gene_name"].unique())
    default_idx = genes.index("flp-9") if "flp-9" in genes else 0
    gene = st.selectbox("Select gene", genes, index=default_idx)

    gdf = df[df["gene_name"] == gene].copy()
    ghvg = gdf[gdf["HVG"] == True].sort_values("Epsilon", ascending=False)

    c1, c2, c3 = st.columns(3)
    c1.metric("HVG in # cell types", len(ghvg))
    c2.metric("Mean Epsilon (HVG)", f"{ghvg['Epsilon'].mean():.2f}" if len(ghvg) else "—")
    c3.metric("Gene class", classify_gene(gene))

    if len(ghvg):
        fig = px.bar(
            ghvg, x="Epsilon", y="cell_type", orientation="h",
            color="neuron_family", color_discrete_map=FAMILY_COLORS,
            hover_data=["Mu","Prob","n_cells"],
            title=f"{gene} — Epsilon across cell types where it is HVG",
        )
        fig.update_layout(height=max(400, len(ghvg)*18), yaxis={"categoryorder":"total ascending"})
        st.plotly_chart(fig, width="stretch")
    else:
        st.info(f"{gene} is not called HVG in any cell type.")

    st.subheader("Per-cell-type values")
    st.dataframe(
        gdf[["cell_type","neuron_family","HVG","Epsilon","Mu","Prob","n_cells"]]
        .sort_values("Epsilon", ascending=False),
        width="stretch", hide_index=True, height=300,
    )

# ── Tab 4: Cell type ─────────────────────────────────────────────────────────
with t4:
    cts = sorted(df["cell_type"].unique())
    default_ct = cts.index("RIA") if "RIA" in cts else 0
    ct = st.selectbox("Select cell type", cts, index=default_ct)
    topn = st.slider("Top N genes", 10, 100, 30, 5)

    cdf = df[df["cell_type"] == ct]
    chvg = cdf[cdf["HVG"] == True]
    if bio_only:
        chvg = chvg[chvg["nuisance"] == "biological"]

    s = summary[summary["cell_type"] == ct]
    c1, c2, c3, c4 = st.columns(4)
    if len(s):
        c1.metric("Cells", int(s["n_cells"].values[0]))
        c2.metric("Batches", int(s["n_experiments"].values[0]))
        c3.metric("HVG total", int(s["n_hvg"].values[0]))
        c4.metric("HVG biological", int(s["n_hvg_biological"].values[0]))

    top = chvg.sort_values("Epsilon", ascending=False).head(topn)
    fig = px.bar(
        top, x="Epsilon", y="gene_name", orientation="h",
        color="gene_class", hover_data=["Mu","Prob"],
        title=f"{ct} — top {topn} HVG genes by Epsilon",
    )
    fig.update_layout(height=max(400, topn*18), yaxis={"categoryorder":"total ascending"})
    st.plotly_chart(fig, width="stretch")

    st.subheader("WormCat enrichment")
    wc_level = st.selectbox("Category level", ["category_1", "category_2", "category_3"], key="wc_cell_level")
    cgo = wormcat_df[(wormcat_df["cell_type"] == ct) & (wormcat_df["ontology"] == wc_level)].copy()
    if len(cgo):
        cgo["-log10(q)"] = -np.log10(cgo["q_value"].clip(lower=1e-300))
        cgo = cgo.sort_values("-log10(q)", ascending=False)
        fig2 = px.bar(cgo.head(20), x="-log10(q)", y="term", orientation="h",
                      color="ontology", hover_data=["term_id","term_selected_genes"])
        fig2.update_layout(height=max(300, min(len(cgo),20)*22), yaxis={"categoryorder":"total ascending"})
        st.plotly_chart(fig2, width="stretch")
        st.dataframe(cgo[["ontology","term_id","term","p_value","q_value","term_selected_genes"]],
                     width="stretch", hide_index=True)
    else:
        st.info(f"No significant supported WormCat categories for {ct} at this level.")

# ── Tab 5: Family ────────────────────────────────────────────────────────────
with t5:
    st.caption("Dot = one (family, gene). Size = fraction of family that's HVG. Color = mean Epsilon.")
    cls = st.selectbox("Gene class", all_classes,
                       index=all_classes.index("neuropeptide") if "neuropeptide" in all_classes else 0)
    min_fam = st.slider("Show genes HVG in ≥ N families", 1, 6, 3)

    sub = df[(df["gene_class"] == cls)]
    if bio_only:
        sub = sub[sub["nuisance"] == "biological"]

    fam_size = df.drop_duplicates("cell_type").groupby("neuron_family")["cell_type"].count()
    agg = (sub.groupby(["neuron_family","gene_name"])
             .agg(n_hvg=("HVG","sum"),
                  mean_eps=("Epsilon", lambda x: x[sub.loc[x.index,"HVG"]].mean()))
             .reset_index())
    agg["frac"] = agg["n_hvg"] / agg["neuron_family"].map(fam_size)
    agg = agg[agg["n_hvg"] > 0]
    gene_fam_breadth = agg.groupby("gene_name")["neuron_family"].nunique()
    keep = gene_fam_breadth[gene_fam_breadth >= min_fam].index
    agg = agg[agg["gene_name"].isin(keep)]
    agg["mean_eps"] = agg["mean_eps"].fillna(0)

    if len(agg):
        gene_order = agg.groupby("gene_name")["n_hvg"].sum().sort_values(ascending=False).index.tolist()
        fig = px.scatter(
            agg, x="gene_name", y="neuron_family", size="frac", color="mean_eps",
            color_continuous_scale="YlOrRd", size_max=22,
            category_orders={"gene_name": gene_order, "neuron_family": all_families},
            hover_data=["n_hvg","frac","mean_eps"],
            title=f"{cls} variability by family (HVG in ≥{min_fam} families)",
        )
        fig.update_layout(height=400, xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")
    else:
        st.info("No genes match the current filters.")

    st.subheader("Neuropeptide family signatures (HVG in ≥50% of family members)")
    st.dataframe(sigs, width="stretch", hide_index=True)

# ── Tab 6: WormCat enrichment ────────────────────────────────────────────────
with t6:
    st.caption("Significant HVG enrichment after requiring detection support in at least two experiments.")
    level = st.selectbox("Category level", ["category_1", "category_2", "category_3"], key="wc_level")
    categories = wormcat_df[wormcat_df["ontology"] == level]
    view = st.radio("View", ["By cell type", "By category (recurrence)"], horizontal=True)
    if categories.empty:
        st.info("No significant supported WormCat categories at this level.")
    elif view == "By cell type":
        ct = st.selectbox("Cell type", sorted(categories["cell_type"].unique()), key="wc_ct")
        d = categories[categories["cell_type"] == ct].copy()
        d["-log10(q)"] = -np.log10(d["q_value"].clip(lower=1e-300))
        st.dataframe(d[["ontology","term_id","term","p_value","q_value","term_selected_genes"]]
                     .sort_values("q_value"), width="stretch", hide_index=True)
    else:
        rec = (categories.groupby(["ontology","term_id","term"])
                    .agg(n_types=("cell_type","nunique"),
                         cell_types=("cell_type", lambda x: ", ".join(sorted(x))))
                    .reset_index().sort_values("n_types", ascending=False))
        topk = st.slider("Show top N recurrent terms", 10, 60, 25)
        fig = px.bar(rec.head(topk), x="n_types", y="term", orientation="h", color="ontology",
                     title="WormCat categories enriched across the most cell types")
        fig.update_layout(height=max(400, topk*20), yaxis={"categoryorder":"total ascending"})
        st.plotly_chart(fig, width="stretch")
        st.dataframe(rec, width="stretch", hide_index=True)

# ── Tab 7: Ligand-Receptor ───────────────────────────────────────────────────
with t7:
    st.caption("Cognate neuropeptide-GPCR pairs where both ligand and receptor are HVG (in some cell type).")
    m = lr[(lr["ligand_mapped"] == True) & lr["ligand_ct"].notna()].copy()
    only_cross = st.toggle("Only cross-cell-type (sender ≠ receiver)", value=True)
    if only_cross:
        m = m[m["same_ct"] == False]
    rec_classes = sorted(m["gpcr_class"].dropna().unique())
    sel_rc = st.multiselect("Receptor class", rec_classes, default=rec_classes)
    m = m[m["gpcr_class"].isin(sel_rc)]
    min_le = st.slider("Min ligand Epsilon", 0.0, float(m["ligand_epsilon"].max() if len(m) else 1), 0.0, 0.5)
    m = m[m["ligand_epsilon"] >= min_le]

    st.metric("Pairs shown", len(m))
    if len(m):
        fig = px.scatter(
            m, x="ligand_epsilon", y="receptor_epsilon",
            color="gpcr_class", hover_data=["receptor","receptor_ct","ligand","ligand_ct"],
            title="Ligand vs receptor Epsilon (each point = one co-variable pair)",
        )
        fig.update_layout(height=500)
        st.plotly_chart(fig, width="stretch")
    st.dataframe(
        m[["receptor","receptor_ct","gpcr_class","receptor_epsilon","ligand","ligand_ct","ligand_epsilon"]]
        .sort_values("ligand_epsilon", ascending=False),
        width="stretch", hide_index=True, height=400,
    )

# ── Tab 8: accepted nuisance-adjusted clusters ───────────────────────────────
with t8:
    st.header("Nuisance-adjusted within-cell-type clusters")
    st.caption(
        "Only statistically accepted multi-cluster results are shown. The UMAP is "
        "computed from the same selected PC scores used by the GMM and is a display "
        "only: it does not determine the clusters or their acceptance."
    )

    if cluster_manifest.empty:
        st.info(
            "No clustering assets were found. Run the R clustering exporter to create "
            "outputs/explorer/clustering/manifest.parquet."
        )
    else:
        cluster_types = cluster_manifest["cell_type"].astype(str).tolist()
        default_cluster_type = cluster_types.index("SMD") if "SMD" in cluster_types else 0
        cluster_type = st.selectbox(
            "Cell type",
            cluster_types,
            index=default_cluster_type,
            key="cluster_cell_type",
        )

        cells, cluster_summary, qc, markers = load_cluster_assets(cluster_type)
        if len(cluster_summary) != 1:
            st.error(f"{cluster_type}: expected exactly one summary row.")
            st.stop()

        cluster_summary = cluster_summary.iloc[0]
        cells["cell_id"] = cells["cell_id"].astype(str)
        cells["cluster"] = cells["cluster"].astype(str)
        markers["feature"] = markers["feature"].astype(str)
        markers["target_cluster"] = markers["target_cluster"].astype(str)
        markers = markers.sort_values("score", ascending=False).reset_index(drop=True)

        # Analysis summary
        s1, s2 = st.columns(2)
        s1.metric("Accepted G", int(cluster_summary["primary_g"]))
        s2.metric("Cells", f"{int(cluster_summary['retained_cells']):,}")
        st.metric("Cluster sizes", str(cluster_summary["cluster_sizes"]))
        s4, s5 = st.columns(2)
        s4.metric("ΔBIC vs G=1", f"{cluster_summary['bic_delta_from_one']:.1f}")
        s5.metric("Stability ARI", f"{cluster_summary['stability_median']:.3f}")
        st.caption(
            f"Cheap-screen choice: G={int(cluster_summary['screen_g'])}; "
            f"selected PCs: {int(cluster_summary['selected_pc_count'])}; "
            f"supported markers: {int(cluster_summary['marker_count'])}."
        )

        # Marker table and filters. Row selection is placed before the gene
        # widget so a selected row can immediately update the UMAP on rerun.
        st.subheader("Marker genes")
        f1, f2 = st.columns(2)
        marker_types = sorted(markers["type"].dropna().unique().tolist())
        selected_marker_types = f1.multiselect(
            "Marker type",
            marker_types,
            default=marker_types,
            key=f"cluster_marker_types_{cluster_type}",
        )
        marker_clusters = sorted(
            markers["target_cluster"].dropna().unique().tolist(),
            key=lambda value: int(value) if str(value).isdigit() else str(value),
        )
        selected_marker_clusters = f2.multiselect(
            "Target cluster",
            marker_clusters,
            default=marker_clusters,
            key=f"cluster_marker_clusters_{cluster_type}",
        )
        f3, f4 = st.columns(2)
        pca_only = f3.toggle(
            "Entered PCA only",
            value=False,
            key=f"cluster_marker_pca_{cluster_type}",
        )
        maximum_q = f4.select_slider(
            "Maximum adjusted p-value",
            options=[0.001, 0.01, 0.05, 0.10, 1.0],
            value=0.10,
            key=f"cluster_marker_q_{cluster_type}",
        )
        marker_search = st.text_input(
            "Search marker name or WBGene ID",
            key=f"cluster_marker_search_{cluster_type}",
        ).strip().lower()

        filtered_markers = markers[
            markers["type"].isin(selected_marker_types)
            & markers["target_cluster"].isin(selected_marker_clusters)
        ].copy()
        filtered_markers["best_q"] = filtered_markers[
            ["detection_q", "magnitude_q"]
        ].min(axis=1, skipna=True)
        filtered_markers = filtered_markers[filtered_markers["best_q"] <= maximum_q]
        if pca_only:
            filtered_markers = filtered_markers[filtered_markers["entered_clustering"] == True]
        if marker_search:
            search_mask = (
                filtered_markers["display_name"].fillna("").str.lower().str.contains(
                    marker_search, regex=False
                )
            )
            filtered_markers = filtered_markers[search_mask]

        marker_columns = [
            "display_name", "type", "target_cluster", "detection_difference",
            "magnitude_difference", "detection_q", "magnitude_q",
            "entered_clustering", "score", "driver_rank", "centroid_contribution",
        ]
        marker_columns = [c for c in marker_columns if c in filtered_markers.columns]
        marker_display = filtered_markers[marker_columns].rename(columns={
            "display_name": "Marker",
            "type": "Type",
            "target_cluster": "Target cluster",
            "detection_difference": "Detection difference",
            "magnitude_difference": "Magnitude difference",
            "detection_q": "Detection q",
            "magnitude_q": "Magnitude q",
            "entered_clustering": "Entered PCA",
            "score": "Marker score",
            "driver_rank": "Driver rank",
            "centroid_contribution": "Centroid contribution",
        })

        marker_event = st.dataframe(
            marker_display,
            width="stretch",
            hide_index=True,
            height=330,
            column_config={
                "Detection difference": st.column_config.NumberColumn(format="%.3f"),
                "Magnitude difference": st.column_config.NumberColumn(format="%.3f"),
                "Detection q": st.column_config.NumberColumn(format="%.2e"),
                "Magnitude q": st.column_config.NumberColumn(format="%.2e"),
                "Marker score": st.column_config.NumberColumn(format="%.2f"),
                "Driver rank": st.column_config.NumberColumn(format="%d"),
                "Centroid contribution": st.column_config.NumberColumn(format="%.3f"),
            },
            on_select="rerun",
            selection_mode="single-row",
            key=f"cluster_marker_table_{cluster_type}",
        )
        selected_rows = marker_event.selection.rows
        selected_from_table = None
        if selected_rows:
            selected_position = selected_rows[0]
            if 0 <= selected_position < len(filtered_markers):
                selected_from_table = filtered_markers.iloc[selected_position]["feature"]

        st.download_button(
            "Download displayed markers",
            filtered_markers.to_csv(index=False).encode(),
            file_name=f"{cluster_type}_cluster_markers.csv",
            mime="text/csv",
            key=f"cluster_marker_download_{cluster_type}",
        )

        # Interactive UMAP
        st.subheader("Interactive clustering UMAP")
        all_features = markers["feature"].tolist()
        label_by_feature = markers.set_index("feature")["display_name"].to_dict()
        gene_key = f"cluster_gene_{cluster_type}"
        mode_key = f"cluster_colour_mode_{cluster_type}"
        if selected_from_table in all_features:
            st.session_state[gene_key] = selected_from_table
            st.session_state[mode_key] = "Gene expression"
        if st.session_state.get(gene_key) not in all_features:
            st.session_state[gene_key] = all_features[0]
        if st.session_state.get(mode_key) not in ["Cluster", "Gene expression"]:
            st.session_state[mode_key] = "Cluster"

        u1, u2, u3 = st.columns([1, 2.2, 1.2])
        colour_mode = u1.radio(
            "Colour cells by",
            ["Cluster", "Gene expression"],
            key=mode_key,
        )
        selected_feature = u2.selectbox(
            "Marker for expression view",
            all_features,
            format_func=lambda feature: label_by_feature.get(feature, feature),
            key=gene_key,
        )
        clip_colour = u3.toggle(
            "Clip positive colours at 99th percentile",
            value=True,
            key=f"cluster_clip_{cluster_type}",
        )

        hover_columns = [
            column for column in [
                "cell_id", "cluster", "Experiment", "Detection",
                "total_features_by_counts", "total_counts", "pct_counts_Mito",
            ] if column in cells.columns
        ]

        if colour_mode == "Cluster":
            plot_cells = cells.copy()
            plot_cells["cluster_label"] = "Cluster " + plot_cells["cluster"]
            cluster_order = sorted(
                plot_cells["cluster_label"].unique(),
                key=lambda value: int(value.split()[-1]) if value.split()[-1].isdigit() else value,
            )
            cluster_fig = px.scatter(
                plot_cells,
                x="UMAP_1",
                y="UMAP_2",
                color="cluster_label",
                category_orders={"cluster_label": cluster_order},
                hover_data=hover_columns,
                labels={"cluster_label": "Cluster"},
                title=f"{cluster_type}: accepted G={int(cluster_summary['primary_g'])}",
            )
            cluster_fig.update_traces(marker={"size": 6, "opacity": 0.82})
        else:
            expression = load_cluster_expression(cluster_type, selected_feature)
            plot_cells = cells.merge(
                expression,
                on="cell_id",
                how="left",
                validate="one_to_one",
            )
            if plot_cells["normalized_expression"].isna().any():
                st.error("Expression cells do not match the UMAP cells.")
                st.stop()

            positive = plot_cells["normalized_expression"] > 0
            positive_values = plot_cells.loc[positive, "normalized_expression"]
            colour_max = (
                positive_values.quantile(0.99)
                if clip_colour and len(positive_values)
                else positive_values.max() if len(positive_values) else 1.0
            )
            colour_max = max(float(colour_max), 1e-12)
            cluster_fig = go.Figure()

            def hover_matrix(frame):
                columns = hover_columns + ["normalized_expression"]
                return frame[columns].astype(object).to_numpy(), columns

            zero_cells = plot_cells.loc[~positive]
            zero_hover, zero_columns = hover_matrix(zero_cells)
            hover_lines = [f"{column}: %{{customdata[{index}]}}" for index, column in enumerate(zero_columns)]
            cluster_fig.add_trace(go.Scattergl(
                x=zero_cells["UMAP_1"],
                y=zero_cells["UMAP_2"],
                mode="markers",
                name="Zero UMI",
                marker={"color": "#c7c7c7", "size": 6, "opacity": 0.65},
                customdata=zero_hover,
                hovertemplate="<br>".join(hover_lines) + "<extra></extra>",
            ))

            positive_cells = plot_cells.loc[positive]
            positive_hover, positive_columns = hover_matrix(positive_cells)
            hover_lines = [f"{column}: %{{customdata[{index}]}}" for index, column in enumerate(positive_columns)]
            cluster_fig.add_trace(go.Scattergl(
                x=positive_cells["UMAP_1"],
                y=positive_cells["UMAP_2"],
                mode="markers",
                name="Positive",
                marker={
                    "color": positive_cells["normalized_expression"].clip(upper=colour_max),
                    "colorscale": "Viridis",
                    "cmin": 0,
                    "cmax": colour_max,
                    "size": 7,
                    "opacity": 0.88,
                    "colorbar": {"title": "log1p(UMI / SF)"},
                },
                customdata=positive_hover,
                hovertemplate="<br>".join(hover_lines) + "<extra></extra>",
            ))
            selected_label = label_by_feature.get(selected_feature, selected_feature)
            cluster_fig.update_layout(title=f"{cluster_type}: {selected_label}")

        cluster_fig.update_layout(
            height=650,
            legend_title_text="Cluster" if colour_mode == "Cluster" else "Detection",
            xaxis={"title": "Clustering UMAP 1", "showgrid": False, "zeroline": False},
            yaxis={"title": "Clustering UMAP 2", "showgrid": False, "zeroline": False},
        )
        st.plotly_chart(cluster_fig, width="stretch")
        st.caption(
            "Grey cells have zero observed UMI for the selected marker. Expression is "
            "log1p(raw UMI / supplied size factor). The original CeNGEN UMAP is not used."
        )

        # BASiCS values for the selected marker and cell type.
        selected_marker = markers.loc[markers["feature"] == selected_feature].iloc[0]
        selected_gene_name = selected_marker.get("gene_short_name")
        basics_match = df[
            (df["cell_type"] == cluster_type)
            & (df["gene_name"].astype(str) == str(selected_gene_name))
        ]
        st.subheader("Selected marker: clustering and BASiCS")
        b1, b2, b3, b4, b5 = st.columns(5)
        b1.metric("Marker type", selected_marker["type"])
        b2.metric("Target cluster", selected_marker["target_cluster"])
        b3.metric("Marker score", f"{selected_marker['score']:.2f}")
        if len(basics_match):
            basics_row = basics_match.sort_values("Prob", ascending=False).iloc[0]
            b4.metric("BASiCS Epsilon", f"{basics_row['Epsilon']:.3f}")
            b5.metric("BASiCS HVG", "Yes" if bool(basics_row["HVG"]) else "No")
            st.caption(
                f"BASiCS posterior probability: {basics_row['Prob']:.3f}; "
                f"nuisance annotation: {basics_row['nuisance']}."
            )
        else:
            b4.metric("BASiCS Epsilon", "—")
            b5.metric("BASiCS HVG", "Not matched")
            st.caption(
                "No BASiCS record matched this cell type and gene symbol. The clustering "
                "result remains available by WBGene ID."
            )

        # Resolution ladder and QC diagnostics.
        st.subheader("Resolution and QC diagnostics")
        st.caption(
            "Cramér's V measures association with categorical nuisances. Eta-squared and "
            "omega-squared measure continuous associations. Mitochondrial percentage is "
            "reported for biological review but is diagnostic-only and does not veto a result."
        )
        qc_columns = [
            "g", "is_primary", "status", "bic_delta_from_one",
            "bic_delta_from_previous", "minimum_cluster_size", "cluster_sizes",
            "stability_median", "structurally_supported", "stable",
            "core_nuisance_clear", "gating_qc_clear", "borderline",
            "cramers_v_experiment", "cramers_v_Detection", "eta_squared_stress",
            "omega_squared_stress", "eta_squared_log_size_factor",
            "omega_squared_log_size_factor", "eta_squared_total_features_by_counts",
            "omega_squared_total_features_by_counts", "eta_squared_pct_counts_Mito",
            "omega_squared_pct_counts_Mito", "mitochondrial_policy",
        ]
        qc_columns = [column for column in qc_columns if column in qc.columns]
        qc_view = qc[qc_columns].copy()

        def colour_qc_row(row):
            if bool(row.get("is_primary", False)):
                colour = "background-color: #d9edf7"
            elif not bool(row.get("structurally_supported", True)) or not bool(row.get("stable", True)):
                colour = "background-color: #f8d7da"
            elif not bool(row.get("core_nuisance_clear", True)):
                colour = "background-color: #fff3cd"
            elif bool(row.get("borderline", False)):
                colour = "background-color: #fff3cd"
            else:
                colour = "background-color: #e8f5e9"
            return [colour] * len(row)

        qc_styled = qc_view.style.apply(colour_qc_row, axis=1).format(precision=3, na_rep="—")
        st.dataframe(qc_styled, width="stretch", hide_index=True, height=420)
        st.download_button(
            "Download resolution/QC table",
            qc.to_csv(index=False).encode(),
            file_name=f"{cluster_type}_resolution_qc.csv",
            mime="text/csv",
            key=f"cluster_qc_download_{cluster_type}",
        )
