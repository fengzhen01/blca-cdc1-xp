args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

e <- new.env(parent = emptyenv())
load(input, envir = e)
x <- e$cds
counts <- get("counts", envir = attr(x, "assayData"), inherits = FALSE)
pheno <- attr(attr(x, "phenoData"), "data")
feature <- attr(attr(x, "featureData"), "data")

cat("counts dim", paste(dim(counts), collapse="x"), "\n")
cat("pheno dim", paste(dim(pheno), collapse="x"), "\n")
cat("feature dim", paste(dim(feature), collapse="x"), "\n")
cat("pheno columns\n")
print(colnames(pheno))
cat("feature columns\n")
print(colnames(feature))
cat("rowname examples", paste(head(rownames(counts)), collapse=", "), "\n")
cat("feature head\n")
print(head(feature))

write.csv(pheno, file.path(outdir, "IMvigor210_pheno.csv"), row.names=TRUE, quote=TRUE)
write.csv(feature, file.path(outdir, "IMvigor210_feature_data.csv"), row.names=TRUE, quote=TRUE)

genes <- c(
  "CLEC9A","XCR1","BATF3","CADM1","IRF8","WDFY4",
  "TAP1","TAP2","TAPBP","B2M","HLA-A","HLA-B","HLA-C","SEC22B",
  "LAMP3","CCR7","FSCN1","CD80","CD86","CD40","IL12B","CCL19",
  "CD8A","CD8B","CXCL9","CXCL10","GZMA","GZMB","PRF1","IFNG","TBX21",
  "STAT1","ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","IDO1",
  "CCL5","CD27","CD274","CMKLR1","CXCR6","HLA-DQA1","HLA-DRB1","LAG3",
  "NKG7","PDCD1LG2","PSMB10","TIGIT","COL1A1","COL1A2","COL3A1","COL5A1",
  "COL6A1","ACTA2","TAGLN","FAP","TGFBI","SPP1","APOE","C1QA","C1QB",
  "C1QC","IL10","TGFB1"
)

# Row names are expected to be gene symbols in this package. If not, search all
# feature columns and retain exact matches.
row_ids <- rownames(counts)
selected_idx <- which(toupper(row_ids) %in% genes)
if (!length(selected_idx) && ncol(feature)) {
  for (j in seq_len(ncol(feature))) {
    vals <- toupper(as.character(feature[[j]]))
    selected_idx <- union(selected_idx, which(vals %in% genes))
  }
}
cat("selected rows", length(selected_idx), "\n")
print(cbind(row=rownames(counts)[selected_idx], feature[selected_idx,,drop=FALSE]))

selected <- counts[selected_idx,,drop=FALSE]
library_size <- colSums(counts)
write.csv(selected, file.path(outdir, "IMvigor210_fixed_gene_counts.csv"), row.names=TRUE, quote=FALSE)
saveRDS(list(counts=selected, library_size=library_size, pheno=pheno, feature=feature[selected_idx,,drop=FALSE]), file.path(outdir,"IMvigor210_fixed_data.rds"))
