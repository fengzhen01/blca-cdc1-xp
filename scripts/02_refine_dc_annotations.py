from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

OUT = ROOT / "results" / "single_cell"
FIG = ROOT / "figures" / "single_cell"
ATLAS = OUT / "GSE293189_standard_atlas.h5ad"
DCFILE = OUT / "GSE293189_myeloid_dc.h5ad"


def savefig(fig, stem):
    for ext, dpi in [("png", 400), ("pdf", 300), ("svg", 300)]:
        fig.savefig(FIG / f"{stem}.{ext}", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def expr(a, gene):
    if gene not in a.var_names:
        return np.zeros(a.n_obs, dtype=np.float32)
    x = a[:, gene].X
    return np.asarray(x.toarray()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()


def label_dc(a):
    genes = ["PTPRC","LST1","TYROBP","HLA-DRA","CD74","CLEC9A","XCR1","CD1C",
             "CLEC10A","FCER1A","CLEC4C","IL3RA","LAMP3","CCR7","FSCN1"]
    x = {g: expr(a, g) for g in genes}
    immune = (x["PTPRC"] > 0) | (x["LST1"] > 0) | (x["TYROBP"] > 0)
    apc = (x["HLA-DRA"] > 0) | (x["CD74"] > 0)
    mature_n = (x["LAMP3"] > 0).astype(int) + (x["CCR7"] > 0).astype(int) + (x["FSCN1"] > 0).astype(int)
    lab = np.full(a.n_obs, "Other myeloid", dtype=object)
    lab[immune & apc & (mature_n >= 2)] = "Mature DC"
    lab[immune & apc & ((x["CLEC4C"] > 0) | (x["IL3RA"] > 0))] = "pDC"
    lab[immune & apc & ((x["CD1C"] > 0) | (x["CLEC10A"] > 0) | (x["FCER1A"] > 0))] = "cDC2"
    lab[immune & apc & ((x["CLEC9A"] > 0) | (x["XCR1"] > 0))] = "cDC1"
    return lab


adata = sc.read_h5ad(ATLAS)
dc = sc.read_h5ad(DCFILE)
adata.obs["high_confidence_DC"] = label_dc(adata)
dc.obs["high_confidence_DC"] = label_dc(dc)

# Recreate patient-level pseudobulk after refined marker definitions.
dc.obs["patient_group"] = (dc.obs["patient"].astype(str) + "__" + dc.obs["specimen"].astype(str) +
                           "__" + dc.obs["high_confidence_DC"].astype(str))
pb_cols, pb_mats, pb_meta = [], [], []
for key, idx in dc.obs.groupby("patient_group", observed=True).indices.items():
    z = dc.obs.iloc[idx]
    group = str(z["high_confidence_DC"].iloc[0])
    summed = np.asarray(dc.layers["counts"][idx].sum(axis=0)).ravel()
    pb_cols.append(key); pb_mats.append(summed)
    pb_meta.append({"pseudobulk_id": key, "patient": z["patient"].iloc[0],
                    "specimen": z["specimen"].iloc[0], "group": group,
                    "n_cells": len(idx), "n_samples": z["gsm"].nunique()})
pd.DataFrame(np.vstack(pb_mats).T.astype(np.int64), index=dc.var_names, columns=pb_cols).to_csv(
    OUT / "patient_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")
pd.DataFrame(pb_meta).to_csv(OUT / "patient_pseudobulk_metadata.tsv", sep="\t", index=False)
adata.obs.to_csv(OUT / "cell_metadata.tsv.gz", sep="\t", compression="gzip")
dc.obs.to_csv(OUT / "dc_cell_metadata.tsv.gz", sep="\t", compression="gzip")

# Updated DC figure.
dum = dc.obsm["X_umap"]
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
pal = plt.get_cmap("tab20")
for i, c in enumerate(sorted(dc.obs["dc_leiden"].astype(str).unique(), key=lambda z: int(z))):
    m = dc.obs["dc_leiden"].astype(str).to_numpy() == c
    axes[0].scatter(dum[m, 0], dum[m, 1], s=4, alpha=.7, color=pal(i % 20), label=c, rasterized=True)
axes[0].set(title="A  Myeloid/DC subclusters", xlabel="UMAP 1", ylabel="UMAP 2")
axes[0].legend(markerscale=3, fontsize=7, frameon=False, ncol=2)
types = ["Other myeloid", "cDC2", "cDC1", "pDC", "Mature DC"]
colors = {"Other myeloid":"#BDBDBD","cDC2":"#4DBBD5","cDC1":"#E64B35","pDC":"#00A087","Mature DC":"#7E6148"}
for c in types:
    m = dc.obs["high_confidence_DC"].astype(str).to_numpy() == c
    axes[1].scatter(dum[m,0], dum[m,1], s=4 if c=="Other myeloid" else 11,
                    alpha=.35 if c=="Other myeloid" else .9, color=colors[c], label=f"{c} (n={m.sum()})", rasterized=True)
axes[1].set(title="B  Marker-supported DC identities", xlabel="UMAP 1", ylabel="UMAP 2")
axes[1].legend(markerscale=2.4, fontsize=7, frameon=False)
cross = dc.obs["score_Cross-presentation"].to_numpy()
p = axes[2].scatter(dum[:,0], dum[:,1], c=cross, cmap="viridis", s=5, alpha=.8, rasterized=True)
axes[2].set(title="C  Locked cross-presentation score", xlabel="UMAP 1", ylabel="UMAP 2")
fig.colorbar(p, ax=axes[2], label="Gene-wise z-score mean")
fig.suptitle("Dendritic-cell and myeloid compartment in GSE293189", fontsize=14)
fig.tight_layout(rect=(0,0,1,.94))
savefig(fig, "Figure7_GSE293189_DC_subclusters")


def dotplot(ax, a, group_key, group_order, genes, title):
    rows = []
    for group in group_order:
        mask = a.obs[group_key].astype(str).to_numpy() == group
        for gene in genes:
            v = expr(a[mask], gene) if mask.any() else np.array([])
            rows.append((group, gene, float((v > 0).mean() * 100) if len(v) else 0,
                         float(v[v > 0].mean()) if (v > 0).any() else 0))
    d = pd.DataFrame(rows, columns=["group","gene","pct","mean"])
    # Color-standardize each gene to emphasize relative, not absolute, expression.
    d["zmean"] = d.groupby("gene")["mean"].transform(lambda x: (x-x.mean())/(x.std(ddof=0) or 1)).clip(-2.5,2.5)
    for yi, group in enumerate(group_order):
        z = d[d.group == group].set_index("gene").loc[genes]
        sca = ax.scatter(range(len(genes)), [yi]*len(genes), s=10 + z.pct.to_numpy()*2.0,
                         c=z.zmean.to_numpy(), cmap="RdBu_r", vmin=-2.5, vmax=2.5, edgecolor="grey", linewidth=.25)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(group_order))); ax.set_yticklabels(group_order, fontsize=8)
    ax.set_xlim(-.7, len(genes)-.3); ax.set_ylim(len(group_order)-.5,-.5); ax.set_title(title)
    ax.grid(axis="x", color="#EEEEEE", lw=.6)
    return sca

fig, axes = plt.subplots(1, 2, figsize=(16, 5.8), gridspec_kw={"width_ratios":[1.05,1.35]})
broad_order = [x for x in ["Epithelial","T cell","NK cell","B cell","Plasma cell","Myeloid/DC","Endothelial","Fibroblast","Pericyte","Mast cell"] if x in adata.obs["cell_type"].astype(str).unique()]
sca = dotplot(axes[0], adata, "cell_type", broad_order,
              ["EPCAM","KRT8","CD3D","NKG7","MS4A1","JCHAIN","LST1","C1QA","PECAM1","COL1A1","RGS5","TPSAB1"],
              "A  Canonical markers of major cell types")
dc_order = [x for x in ["cDC1","cDC2","pDC","Mature DC","Other myeloid"] if x in dc.obs["high_confidence_DC"].astype(str).unique()]
sca2 = dotplot(axes[1], dc, "high_confidence_DC", dc_order,
               ["XCR1","CLEC9A","WDFY4","CD1C","CLEC10A","FCER1A","CLEC4C","IL3RA","GZMB","LAMP3","CCR7","FSCN1","C1QA","SPP1","TAP1","TAPBP","SEC22B"],
               "B  Identity and functional markers in DC/myeloid cells")
cbar = fig.colorbar(sca2, ax=axes, fraction=.018, pad=.02); cbar.set_label("Relative mean expression (z)")
for size, label in [(30,"10%"),(110,"50%"),(210,"100%")]:
    axes[1].scatter([], [], s=size, c="white", edgecolor="grey", label=label)
axes[1].legend(title="Expressing cells", loc="upper left", bbox_to_anchor=(1.01,.70), frameon=False, fontsize=8)
fig.suptitle("Marker-based validation of single-cell annotations", fontsize=14)
fig.subplots_adjust(top=.88, bottom=.24, left=.09, right=.88, wspace=.30)
savefig(fig, "FigureS3_annotation_marker_dotplots")

summary_path = OUT / "analysis_summary.json"
summary = json.loads(summary_path.read_text(encoding="utf-8"))
summary.update({
    "tumor_cDC1_cells_refined": int(((dc.obs.specimen=="tumor") & (dc.obs.high_confidence_DC=="cDC1")).sum()),
    "tumor_cDC1_patients_refined": int(dc.obs.loc[(dc.obs.specimen=="tumor") & (dc.obs.high_confidence_DC=="cDC1"), "patient"].nunique()),
    "DC_definition_note": "cDC1/cDC2/pDC require immune+APC evidence and subtype markers; mature DC requires at least two of LAMP3/CCR7/FSCN1."
})
summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

# Keep serialized objects synchronized with the refined labels.
adata.write_h5ad(ATLAS, compression="gzip")
dc.write_h5ad(DCFILE, compression="gzip")
print(pd.crosstab(dc.obs["specimen"], dc.obs["high_confidence_DC"]))
