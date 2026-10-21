# Reproduction guide

## Software

The completed analysis used Python 3.12 / Scanpy 1.12.4 and R 4.4.1 with edgeR 4.4.2, limma 3.62.2 and survival 3.6.4. `environment.yml` provides a portable environment specification.

## Analysis order

### 1. Single-cell atlas

```bash
python scripts/01_standard_single_cell_atlas.py
```

Main decisions:

- retain cells with 150–5,000 detected genes;
- require at least 300 total counts;
- mitochondrial fraction ≤25%;
- run Scrublet per library when the library has sufficient cells;
- normalize to 10,000 counts, log-transform and select 3,000 highly variable genes;
- correct source batches using ComBat;
- build a 20-neighbor graph on 35 PCs;
- Leiden resolution 0.65 and UMAP minimum distance 0.35;
- random seed 20261020.

These permissive count thresholds are deliberate for mixed-source/Seq-Well libraries. A stricter 200-gene/500-count sensitivity run removed rare cDC1 cells disproportionately.

### 2. Refined DC annotations

```bash
python scripts/02_refine_dc_annotations.py
```

cDC1/cDC2/pDC calls require immune/APC evidence plus subtype markers. Mature DC requires at least two of `LAMP3`, `CCR7` and `FSCN1`. Marker dot plots are generated as an annotation audit.

### 3. Patient-paired pseudobulk

If the full single-cell workflow was run, the pseudobulk files are generated automatically. Otherwise, use the provided files under `data/derived/GSE293189/`.

```bash
Rscript scripts/03_pseudobulk_de_enrichment.R
```

The primary cDC1-versus-cDC2 analysis requires at least three cDC1 and ten cDC2 cells per patient. Eight paired patients pass this rule. A sensitivity analysis uses a three-cell comparator threshold. The model uses an edgeR quasi-likelihood GLM with patient blocking. GO Biological Process 2023 and KEGG 2021 Human gene sets are needed under `data/gene_sets/`.

### 4. IMvigor210

```bash
Rscript scripts/04_extract_imvigor210.R data/raw/IMvigor210CoreBiologies/data/cds.RData results/IMvigor210_validation
Rscript scripts/05_imvigor210_locked_scores.R results/IMvigor210_validation/IMvigor210_fixed_data.rds results/IMvigor210_validation
Rscript scripts/06_imvigor210_immune_context.R
```

The score uses gene-wise standardized log-CPM expression. Response is CR/PR versus SD/PD. Logistic regression reports response odds ratios; Cox regression reports overall-survival hazard ratios. Adjusted models test TMB, cDC1 identity and the T-cell-inflamed score separately because these scores are correlated.

## Expected checkpoints

- final GSE293189 cells: 24,715;
- detected genes: 25,594;
- Leiden clusters: 32;
- myeloid/DC cells: 10,670;
- tumor cDC1 cells: 82 from 14 patients;
- primary pseudobulk comparison: 8 paired patients;
- primary cDC1-vs-cDC2 DE genes: 779 at FDR<0.05.

Small numerical differences may occur with different UMAP/graph-library builds, but patient-level DE direction and downstream score conclusions should be stable.

## Validation

Check that each manuscript figure exists in PNG, PDF and SVG formats. Statistical outputs are CSV/TSV files under `results/`. Do not interpret the UMAP layout itself as evidence of biological distance.

