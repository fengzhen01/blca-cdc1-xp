from __future__ import annotations

import gzip
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

RAW = ROOT / "data" / "raw" / "GSE293189"
META_FILE = ROOT / "data" / "raw" / "GSE293189_sample_metadata.tsv"
OUT = ROOT / "results" / "single_cell"
FIG = ROOT / "figures" / "single_cell"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20261020
QC_MIN_GENES = 150
QC_MAX_GENES = 5000
QC_MIN_COUNTS = 300
QC_MAX_MITO = 0.25

GENE_SETS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "KRT17", "KRT20", "KRT7"],
    "T cell": ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "IL7R", "CCL5"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "KLRB1", "CCL5", "PRF1"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD37", "CD22", "CD19", "CD83"],
    "Plasma cell": ["JCHAIN", "MZB1", "DERL3", "SDC1", "IGKC"],
    "Myeloid/DC": ["LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "LILRB1", "LGALS3", "CD68"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "PLVAP"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "COL3A1", "COL6A1", "PDGFRA", "LUM"],
    "Pericyte": ["RGS5", "MCAM", "CSPG4", "NOTCH3", "ACTA2"],
    "Mast cell": ["TPSAB1", "TPSB2", "KIT", "MS4A2", "CPA3"],
    "MHC-II APC": ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74"],
    "Macrophage": ["C1QA", "C1QB", "C1QC", "APOE", "FCGR3A", "CTSD", "SPP1", "CD68"],
    "cDC1": ["CLEC9A", "XCR1", "BATF3", "CADM1", "IRF8", "WDFY4"],
    "cDC2": ["CD1C", "CLEC10A", "FCER1A", "CD1E", "CD1A"],
    "pDC": ["CLEC4C", "IL3RA", "GZMB", "JCHAIN", "TCF4", "IRF7"],
    "Mature DC": ["LAMP3", "CCR7", "FSCN1", "CD80", "CD86", "CD40", "IL12B", "CCL19"],
    "Cross-presentation": ["TAP1", "TAP2", "TAPBP", "SEC22B"],
    "IFN response": ["STAT1", "IRF1", "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"],
}

meta = pd.read_csv(META_FILE, sep="\t").set_index("gsm")


def savefig(fig: plt.Figure, stem: str) -> None:
    for ext, dpi in [("png", 400), ("pdf", 300), ("svg", 300)]:
        fig.savefig(FIG / f"{stem}.{ext}", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def sample_info(path: Path) -> dict[str, str]:
    gsm = re.match(r"(GSM\d+)", path.name).group(1)
    title = str(meta.loc[gsm, "title"])
    pm = re.search(r"Patient\s+([^ ]+)", title, re.I)
    patient = pm.group(1) if pm else gsm
    hm = re.search(r"\(([^)]+)\)", title)
    histology = hm.group(1) if hm else "unknown"
    specimen = "normal" if "normal" in title.lower() else "tumor"
    source = re.sub(r"^GSM\d+_", "", path.name).split("_")[0].split("-")[0]
    if source not in {"HY", "PA", "SG"}:
        source = "other"
    return {"gsm": gsm, "patient": patient, "histology": histology,
            "specimen": specimen, "batch_source": source, "title": title}


def first_pass(path: Path) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, list[str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        barcodes = header[1:]
        n = len(barcodes)
        total = np.zeros(n, dtype=np.float64)
        detected = np.zeros(n, dtype=np.int32)
        mito = np.zeros(n, dtype=np.float64)
        seen_genes = []
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != n + 1:
                continue
            values = np.asarray(parts[1:], dtype=np.float32)
            total += values
            detected += values > 0
            if parts[0].upper().startswith("MT-"):
                mito += values
            seen_genes.append(parts[0])
    return barcodes, total, detected, mito, seen_genes


def second_pass(path: Path, keep: np.ndarray, expected_genes: list[str] | None) -> tuple[sparse.csr_matrix, list[str]]:
    rows, cols, vals, genes = [], [], [], []
    gene_to_col = None if expected_genes is None else {g: i for i, g in enumerate(expected_genes)}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        gi = 0
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0]
            if gene_to_col is not None and gene not in gene_to_col:
                continue
            values = np.asarray(parts[1:], dtype=np.float32)[keep]
            nz = np.flatnonzero(values)
            if nz.size:
                rows.append(nz.astype(np.int32, copy=False))
                col = gi if gene_to_col is None else gene_to_col[gene]
                cols.append(np.full(nz.size, col, dtype=np.int32))
                vals.append(values[nz])
            genes.append(gene)
            gi += 1
    output_genes = genes if expected_genes is None else expected_genes
    if rows:
        rr, cc, vv = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
        mat = sparse.csr_matrix((vv, (rr, cc)), shape=(int(keep.sum()), len(output_genes)), dtype=np.float32)
    else:
        mat = sparse.csr_matrix((int(keep.sum()), len(output_genes)), dtype=np.float32)
    return mat, output_genes


def module_zscores(adata: sc.AnnData, sets: dict[str, list[str]]) -> pd.DataFrame:
    scores = pd.DataFrame(index=adata.obs_names)
    for name, genes in sets.items():
        use = [g for g in genes if g in adata.var_names]
        if not use:
            scores[name] = 0.0
            continue
        arr = adata[:, use].X
        arr = arr.toarray() if sparse.issparse(arr) else np.asarray(arr)
        mu = arr.mean(axis=0, keepdims=True)
        sd = arr.std(axis=0, ddof=1, keepdims=True)
        sd[sd == 0] = 1
        scores[name] = np.clip((arr - mu) / sd, -5, 5).mean(axis=1)
    return scores


def direct_expression(adata: sc.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=np.float32)
    x = adata[:, gene].X
    return np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()


print("Reading and quality-controlling GSE293189 ...", flush=True)
paths = sorted(RAW.glob("*.gz"))
prepass, gene_union = [], set()
for path in paths:
    info = sample_info(path)
    barcodes, total, detected, mito, sample_genes = first_pass(path)
    gene_union.update(sample_genes)
    prepass.append((path, info, barcodes, total, detected, mito))

genes = sorted(gene_union)
print(f"Union feature space: {len(genes)} genes", flush=True)
matrices, obs_parts, qc_rows = [], [], []
for path, info, barcodes, total, detected, mito in prepass:
    mito_frac = np.divide(mito, total, out=np.zeros_like(mito), where=total > 0)
    keep = ((detected >= QC_MIN_GENES) & (detected <= QC_MAX_GENES) &
            (total >= QC_MIN_COUNTS) & (mito_frac <= QC_MAX_MITO))
    mat, _ = second_pass(path, keep, genes)
    matrices.append(mat)
    ids = [f"{info['gsm']}:{barcodes[i]}" for i in np.flatnonzero(keep)]
    obs = pd.DataFrame(index=ids)
    for k, v in info.items():
        obs[k] = v
    obs["barcode"] = [barcodes[i] for i in np.flatnonzero(keep)]
    obs["total_counts"] = total[keep]
    obs["n_genes_by_counts"] = detected[keep]
    obs["pct_counts_mt"] = mito_frac[keep] * 100
    obs_parts.append(obs)
    qc_rows.append({**info, "raw_cells": len(barcodes), "qc_cells": int(keep.sum()),
                    "median_counts_raw": float(np.median(total)),
                    "median_genes_raw": float(np.median(detected)),
                    "median_mito_raw": float(np.median(mito_frac) * 100)})
    print(f"  {info['gsm']}: {len(barcodes)} -> {keep.sum()} cells", flush=True)

X = sparse.vstack(matrices, format="csr")
obs = pd.concat(obs_parts, axis=0)
adata = sc.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()
sc.pp.filter_genes(adata, min_cells=10)

# Doublet detection is performed per library. Very small libraries are retained and
# explicitly marked as not estimable rather than being forced through Scrublet.
adata.obs["doublet_score"] = np.nan
adata.obs["predicted_doublet"] = False
try:
    import scrublet as scr
    for gsm in adata.obs["gsm"].unique():
        ix = np.flatnonzero((adata.obs["gsm"] == gsm).to_numpy())
        if len(ix) < 100:
            continue
        scrub = scr.Scrublet(adata.layers["counts"][ix], expected_doublet_rate=0.06, random_state=SEED)
        scores, pred = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                             min_gene_variability_pctl=85,
                                             n_prin_comps=min(30, len(ix) - 1), verbose=False)
        adata.obs.iloc[ix, adata.obs.columns.get_loc("doublet_score")] = scores
        adata.obs.iloc[ix, adata.obs.columns.get_loc("predicted_doublet")] = pred
except Exception as exc:
    print(f"Scrublet warning: {exc}", flush=True)

n_doublets = int(adata.obs["predicted_doublet"].sum())
adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()
print(f"After QC/doublet filtering: {adata.n_obs} cells, {adata.n_vars} genes; doublets removed={n_doublets}", flush=True)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="batch_source")
analysis = adata[:, adata.var["highly_variable"]].copy()
sc.pp.combat(analysis, key="batch_source")
analysis.X = np.asarray(analysis.X).copy()
sc.pp.scale(analysis, max_value=10)
sc.tl.pca(analysis, n_comps=50, svd_solver="arpack", random_state=SEED)
sc.pp.neighbors(analysis, n_neighbors=20, n_pcs=35, random_state=SEED)
sc.tl.umap(analysis, min_dist=0.35, random_state=SEED)
sc.tl.leiden(analysis, resolution=0.65, key_added="leiden", random_state=SEED)

adata.obsm["X_pca"] = analysis.obsm["X_pca"]
adata.obsm["X_umap"] = analysis.obsm["X_umap"]
adata.obs["leiden"] = analysis.obs["leiden"].astype(str)

scores = module_zscores(adata, GENE_SETS)
broad_names = ["Epithelial", "T cell", "NK cell", "B cell", "Plasma cell", "Myeloid/DC",
               "Endothelial", "Fibroblast", "Pericyte", "Mast cell"]
cluster_scores = scores.join(adata.obs["leiden"]).groupby("leiden", observed=True)[broad_names].mean()
cluster_labels = cluster_scores.idxmax(axis=1).to_dict()
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_labels).astype("category")
for name in GENE_SETS:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    adata.obs[f"score_{safe_name}"] = scores[name].to_numpy()

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", pts=True, key_added="rank_cluster")
sc.get.rank_genes_groups_df(adata, group=None, key="rank_cluster").to_csv(
    OUT / "cluster_marker_genes.csv.gz", index=False, compression="gzip")
cluster_scores.assign(cell_type=pd.Series(cluster_labels)).to_csv(OUT / "cluster_annotation_scores.csv")

# High-confidence marker calls are retained alongside unsupervised clusters for DC comparisons.
expr = {g: direct_expression(adata, g) for g in
        ["PTPRC", "LST1", "TYROBP", "HLA-DRA", "CD74", "CLEC9A", "XCR1", "CD1C",
         "CLEC10A", "FCER1A", "CLEC4C", "IL3RA", "LAMP3", "CCR7", "FSCN1"]}
immune = (expr["PTPRC"] > 0) | (expr["LST1"] > 0) | (expr["TYROBP"] > 0)
apc = (expr["HLA-DRA"] > 0) | (expr["CD74"] > 0)
dc_label = np.full(adata.n_obs, "Other myeloid", dtype=object)
dc_label[immune & apc & ((expr["LAMP3"] > 0) | (expr["CCR7"] > 0) | (expr["FSCN1"] > 0))] = "Mature DC"
dc_label[immune & apc & ((expr["CLEC4C"] > 0) | (expr["IL3RA"] > 0))] = "pDC"
dc_label[immune & apc & ((expr["CD1C"] > 0) | (expr["CLEC10A"] > 0) | (expr["FCER1A"] > 0))] = "cDC2"
dc_label[immune & apc & ((expr["CLEC9A"] > 0) | (expr["XCR1"] > 0))] = "cDC1"
adata.obs["high_confidence_DC"] = dc_label

myeloid_mask = (adata.obs["cell_type"].astype(str) == "Myeloid/DC") | (immune & apc)
dc = adata[myeloid_mask].copy()
sc.pp.highly_variable_genes(dc, n_top_genes=min(2000, dc.n_vars - 1), flavor="seurat", batch_key="batch_source")
dc_work = dc[:, dc.var["highly_variable"]].copy()
sc.pp.combat(dc_work, key="batch_source")
dc_work.X = np.asarray(dc_work.X).copy()
sc.pp.scale(dc_work, max_value=10)
sc.tl.pca(dc_work, n_comps=min(40, dc_work.n_obs - 2, dc_work.n_vars - 1), svd_solver="arpack", random_state=SEED)
sc.pp.neighbors(dc_work, n_neighbors=min(15, dc_work.n_obs - 1), n_pcs=min(25, dc_work.obsm["X_pca"].shape[1]), random_state=SEED)
sc.tl.umap(dc_work, min_dist=0.30, random_state=SEED)
sc.tl.leiden(dc_work, resolution=0.55, key_added="dc_leiden", random_state=SEED)
dc.obsm["X_umap"] = dc_work.obsm["X_umap"]
dc.obs["dc_leiden"] = dc_work.obs["dc_leiden"].astype(str)

dc_scores = module_zscores(dc, {k: GENE_SETS[k] for k in ["Macrophage", "cDC1", "cDC2", "pDC", "Mature DC"]})
dc_cluster_scores = dc_scores.join(dc.obs["dc_leiden"]).groupby("dc_leiden", observed=True).mean()
dc_cluster_labels = dc_cluster_scores.idxmax(axis=1).to_dict()
dc.obs["dc_cluster_type"] = dc.obs["dc_leiden"].map(dc_cluster_labels).astype("category")
dc_cluster_scores.assign(dc_cluster_type=pd.Series(dc_cluster_labels)).to_csv(OUT / "dc_cluster_annotation_scores.csv")

# Pseudobulk raw counts: biological unit is patient, not cell.
use_groups = ["cDC1", "cDC2", "pDC", "Mature DC", "Other myeloid"]
dc.obs["patient_group"] = (dc.obs["patient"].astype(str) + "__" + dc.obs["specimen"].astype(str) +
                           "__" + dc.obs["high_confidence_DC"].astype(str))
pb_cols, pb_mats, pb_meta = [], [], []
for key, idx in dc.obs.groupby("patient_group", observed=True).indices.items():
    subobs = dc.obs.iloc[idx]
    group = str(subobs["high_confidence_DC"].iloc[0])
    if group not in use_groups:
        continue
    summed = np.asarray(dc.layers["counts"][idx].sum(axis=0)).ravel()
    pb_cols.append(key)
    pb_mats.append(summed)
    pb_meta.append({"pseudobulk_id": key, "patient": subobs["patient"].iloc[0],
                    "specimen": subobs["specimen"].iloc[0], "group": group,
                    "n_cells": len(idx), "n_samples": subobs["gsm"].nunique()})
pb = np.vstack(pb_mats).T.astype(np.int64)
pd.DataFrame(pb, index=dc.var_names, columns=pb_cols).to_csv(OUT / "patient_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")
pd.DataFrame(pb_meta).to_csv(OUT / "patient_pseudobulk_metadata.tsv", sep="\t", index=False)

# Tables and serialized data.
pd.DataFrame(qc_rows).to_csv(OUT / "sample_qc_summary.tsv", sep="\t", index=False)
adata.obs.to_csv(OUT / "cell_metadata.tsv.gz", sep="\t", compression="gzip")
dc.obs.to_csv(OUT / "dc_cell_metadata.tsv.gz", sep="\t", compression="gzip")
adata.write_h5ad(OUT / "GSE293189_standard_atlas.h5ad", compression="gzip")
dc.write_h5ad(OUT / "GSE293189_myeloid_dc.h5ad", compression="gzip")

# Figure 6: standard atlas.
fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.2))
axes[0, 0].hist(adata.obs["n_genes_by_counts"], bins=50, color="#4472C4", alpha=.85)
axes[0, 0].axvline(QC_MIN_GENES, color="#C00000", ls="--", lw=1)
axes[0, 0].set(xlabel="Detected genes per cell", ylabel="Cells", title="A  Post-QC library complexity")
um = adata.obsm["X_umap"]
cats = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
palette = plt.get_cmap("tab20")
for i, c in enumerate(cats):
    m = adata.obs["leiden"].to_numpy() == c
    axes[0, 1].scatter(um[m, 0], um[m, 1], s=1.5, alpha=.65, color=palette(i % 20), label=c, rasterized=True)
axes[0, 1].set(title="B  Unsupervised Leiden clusters", xlabel="UMAP 1", ylabel="UMAP 2")
axes[0, 1].legend(markerscale=4, fontsize=7, ncol=2, frameon=False)
cell_types = list(adata.obs["cell_type"].cat.categories)
for i, c in enumerate(cell_types):
    m = adata.obs["cell_type"].astype(str).to_numpy() == c
    axes[1, 0].scatter(um[m, 0], um[m, 1], s=1.7, alpha=.7, color=palette(i % 20), label=c, rasterized=True)
axes[1, 0].set(title="C  Major cell types", xlabel="UMAP 1", ylabel="UMAP 2")
axes[1, 0].legend(markerscale=4, fontsize=7, ncol=2, frameon=False)
for i, source in enumerate(sorted(adata.obs["batch_source"].unique())):
    m = adata.obs["batch_source"].to_numpy() == source
    axes[1, 1].scatter(um[m, 0], um[m, 1], s=1.5, alpha=.65, color=palette(i), label=source, rasterized=True)
axes[1, 1].set(title="D  Source batch after correction", xlabel="UMAP 1", ylabel="UMAP 2")
axes[1, 1].legend(markerscale=4, fontsize=8, frameon=False)
fig.suptitle("GSE293189 single-cell atlas: quality control, clustering and annotation", fontsize=15)
fig.tight_layout(rect=(0, 0, 1, .97))
savefig(fig, "Figure6_GSE293189_standard_atlas")

# Figure 7: dendritic/myeloid compartment.
dum = dc.obsm["X_umap"]
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
for i, c in enumerate(sorted(dc.obs["dc_leiden"].unique(), key=lambda x: int(x))):
    m = dc.obs["dc_leiden"].to_numpy() == c
    axes[0].scatter(dum[m, 0], dum[m, 1], s=4, alpha=.7, color=palette(i % 20), label=c, rasterized=True)
axes[0].set(title="A  Myeloid/DC subclusters", xlabel="UMAP 1", ylabel="UMAP 2")
axes[0].legend(markerscale=3, fontsize=7, frameon=False, ncol=2)
dc_types = ["Other myeloid", "cDC2", "cDC1", "pDC", "Mature DC"]
dc_colors = {"Other myeloid": "#BDBDBD", "cDC2": "#4DBBD5", "cDC1": "#E64B35", "pDC": "#00A087", "Mature DC": "#7E6148"}
for c in dc_types:
    m = dc.obs["high_confidence_DC"].astype(str).to_numpy() == c
    axes[1].scatter(dum[m, 0], dum[m, 1], s=4 if c == "Other myeloid" else 10,
                    alpha=.45 if c == "Other myeloid" else .85, color=dc_colors[c], label=c, rasterized=True)
axes[1].set(title="B  High-confidence DC identities", xlabel="UMAP 1", ylabel="UMAP 2")
axes[1].legend(markerscale=2.5, fontsize=7, frameon=False)
cross = dc.obs["score_Cross-presentation"].to_numpy()
p = axes[2].scatter(dum[:, 0], dum[:, 1], c=cross, cmap="viridis", s=5, alpha=.8, rasterized=True)
axes[2].set(title="C  Locked cross-presentation score", xlabel="UMAP 1", ylabel="UMAP 2")
fig.colorbar(p, ax=axes[2], label="Gene-wise z-score mean")
fig.suptitle("Dendritic-cell and myeloid compartment in GSE293189", fontsize=14)
fig.tight_layout(rect=(0, 0, 1, .94))
savefig(fig, "Figure7_GSE293189_DC_subclusters")

summary = {
    "raw_cells": int(sum(x["raw_cells"] for x in qc_rows)),
    "qc_cells_before_doublet": int(sum(x["qc_cells"] for x in qc_rows)),
    "doublets_removed": n_doublets,
    "final_cells": int(adata.n_obs), "genes": int(adata.n_vars),
    "patients": int(adata.obs["patient"].nunique()), "samples": int(adata.obs["gsm"].nunique()),
    "clusters": int(adata.obs["leiden"].nunique()),
    "myeloid_dc_cells": int(dc.n_obs),
    "tumor_cDC1_cells": int(((dc.obs["specimen"] == "tumor") & (dc.obs["high_confidence_DC"] == "cDC1")).sum()),
    "tumor_cDC1_patients": int(dc.obs.loc[(dc.obs["specimen"] == "tumor") & (dc.obs["high_confidence_DC"] == "cDC1"), "patient"].nunique()),
    "qc_thresholds": {"min_genes": QC_MIN_GENES, "max_genes": QC_MAX_GENES,
                      "min_counts": QC_MIN_COUNTS, "max_mito_percent": QC_MAX_MITO * 100},
}
(OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
