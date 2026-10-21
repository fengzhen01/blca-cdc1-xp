# Data access and provenance

## GSE293189: single-cell bladder cancer dataset

- Repository: NCBI Gene Expression Omnibus (GEO)
- Accession: **GSE293189**
- Role: standard single-cell atlas, myeloid/DC reclustering and patient-level pseudobulk analysis
- Expected local layout: `data/raw/GSE293189/`
- Required files: the per-library compressed DGE matrices and a sample metadata table containing GSM, patient, specimen and histology fields

The standard workflow constructs the union of features across libraries before matrix concatenation. This is required because the gene lists are not identical across all libraries; using the first library as the feature reference can remove rare DC markers such as `CLEC9A` and `XCR1`.

## IMvigor210

- Role: validation of locked immune scores against response, overall survival, TMB and tumor immune phenotype
- Expected local source: `data/raw/IMvigor210CoreBiologies/data/cds.RData`
- Redistribution: not included in this repository. Obtain the dataset from its original public distribution and comply with the original terms.

## Additional validation cohorts

The result tables and manuscript figures may reference GSE176307, GSE328930 and GSE135337 analyses. Raw inputs are not redistributed here. Retrieve them from GEO using their accession identifiers before rerunning the corresponding exploratory/validation scripts.

## Included derived data

`data/derived/GSE293189/` contains compact patient-level pseudobulk counts, pseudobulk metadata, sample-level QC summaries and analysis metadata. These files are sufficient to rerun the downstream edgeR and enrichment analyses without downloading the full single-cell matrix.

## Excluded files

- raw sequencing/expression matrices;
- full AnnData/H5AD objects (approximately 592 MB combined);
- third-party software packages and local dependency trees;
- article PDFs and copyrighted supplementary files;
- temporary range-download fragments and logs.

No patient names, direct identifiers or controlled clinical records are included. Public cohort sample identifiers are preserved only where required for reproducibility.

