options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(edgeR); library(limma)})

root <- "."
indir <- file.path(root, "data", "derived", "GSE293189")
outdir <- file.path(root, "results", "pseudobulk")
figdir <- file.path(root, "figures", "pseudobulk")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(figdir, recursive=TRUE, showWarnings=FALSE)

cts <- read.delim(gzfile(file.path(indir, "patient_pseudobulk_counts.tsv.gz")),
                  row.names=1, check.names=FALSE)
meta <- read.delim(file.path(indir, "patient_pseudobulk_metadata.tsv"), check.names=FALSE)
stopifnot(all(meta$pseudobulk_id %in% colnames(cts)))

read_gmt <- function(path) {
  x <- readLines(path, warn=FALSE)
  parts <- strsplit(x, "\t", fixed=FALSE)
  sets <- lapply(parts, function(z) unique(z[-c(1,2)]))
  names(sets) <- vapply(parts, `[`, "", 1)
  sets
}
go_sets <- read_gmt(file.path(root, "data", "gene_sets", "GO_Biological_Process_2023.gmt"))
kegg_sets <- read_gmt(file.path(root, "data", "gene_sets", "KEGG_2021_Human.gmt"))

run_pair <- function(comp, comparator, min_target=3, min_comp=10, tag="primary") {
  m <- meta[meta$specimen == "tumor" & meta$group %in% c(comp, comparator), ]
  keep_pat <- unique(m$patient)
  keep_pat <- keep_pat[vapply(keep_pat, function(p) {
    z <- m[m$patient == p, ]
    any(z$group == comp & z$n_cells >= min_target) &&
      any(z$group == comparator & z$n_cells >= min_comp)
  }, logical(1))]
  m <- m[m$patient %in% keep_pat, ]
  m <- m[order(m$patient, m$group), ]
  if (length(keep_pat) < 4) stop("Too few paired patients for ", comp, " vs ", comparator)
  x <- round(as.matrix(cts[, m$pseudobulk_id, drop=FALSE]))
  y <- DGEList(x)
  m$patient <- factor(m$patient)
  m$group <- relevel(factor(m$group), ref=comparator)
  design <- model.matrix(~patient + group, data=m)
  keep_gene <- filterByExpr(y, design=design, min.count=3, min.total.count=10)
  y <- y[keep_gene, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  coef_name <- grep("^group", colnames(design), value=TRUE)[1]
  qlf <- glmQLFTest(fit, coef=which(colnames(design) == coef_name))
  tab <- topTags(qlf, n=Inf, sort.by="none")$table
  tab$gene <- rownames(tab)
  tab <- tab[order(tab$PValue), c("gene","logFC","logCPM","F","PValue","FDR")]
  prefix <- paste0(tag, "_", comp, "_vs_", gsub(" ", "_", comparator))
  write.csv(tab, file.path(outdir, paste0(prefix, "_edgeR.csv")), row.names=FALSE)
  write.table(m, file.path(outdir, paste0(prefix, "_samples.tsv")), sep="\t", quote=FALSE, row.names=FALSE)

  stat <- sign(tab$logFC) * sqrt(pmax(tab$F, 0)); names(stat) <- tab$gene
  stat <- stat[is.finite(stat)]
  run_camera <- function(sets, collection) {
    idx <- lapply(sets, function(g) which(names(stat) %in% g))
    idx <- idx[vapply(idx, length, integer(1)) >= 10 & vapply(idx, length, integer(1)) <= 500]
    z <- cameraPR(stat, idx, use.ranks=TRUE, sort=TRUE)
    z$pathway <- rownames(z); rownames(z) <- NULL
    z$collection <- collection
    z[, c("collection","pathway","NGenes","Direction","PValue","FDR")]
  }
  enrich <- rbind(run_camera(go_sets, "GO Biological Process 2023"),
                  run_camera(kegg_sets, "KEGG 2021 Human"))
  enrich <- enrich[order(enrich$FDR, enrich$PValue), ]
  write.csv(enrich, file.path(outdir, paste0(prefix, "_ranked_enrichment.csv")), row.names=FALSE)
  list(tab=tab, enrich=enrich, meta=m, prefix=prefix,
       summary=data.frame(comparison=paste(comp,"vs",comparator), analysis=tag,
                          patients=length(keep_pat), genes=nrow(tab),
                          FDR_0.05=sum(tab$FDR < .05), FDR_0.10=sum(tab$FDR < .10),
                          pathways_FDR_0.05=sum(enrich$FDR < .05),
                          pathways_FDR_0.10=sum(enrich$FDR < .10)))
}

res1 <- run_pair("cDC1", "cDC2", 3, 10, "primary")
res2 <- run_pair("cDC1", "Mature DC", 3, 3, "exploratory")
res3 <- run_pair("cDC1", "cDC2", 3, 3, "sensitivity")
summary <- do.call(rbind, lapply(list(res1,res2,res3), `[[`, "summary"))
write.csv(summary, file.path(outdir, "pseudobulk_analysis_summary.csv"), row.names=FALSE)

plot_figure <- function(device, path, width, height, res=NULL) {
  if (device == "png") png(path, width=width*res, height=height*res, res=res)
  if (device == "pdf") pdf(path, width=width, height=height, useDingbats=FALSE)
  if (device == "svg") svg(path, width=width, height=height)
  op <- par(mfrow=c(1,2), mar=c(5,5,3,1), family="sans")
  tab <- res1$tab
  sig <- tab$FDR < .10 & abs(tab$logFC) >= .5
  col <- ifelse(sig & tab$logFC > 0, "#D73027", ifelse(sig, "#4575B4", "#BDBDBD"))
  plot(tab$logFC, -log10(pmax(tab$FDR, 1e-300)), pch=16, cex=.45, col=col,
       xlab=expression(log[2]~fold~change~"(cDC1 / cDC2)"), ylab=expression(-log[10]~FDR),
       main=paste0("A  Patient-paired pseudobulk DE (n=", res1$summary$patients, ")"))
  abline(v=c(-.5,.5), h=-log10(.10), lty=2, col="grey45")
  lab <- head(tab[order(tab$FDR, -abs(tab$logFC)), ], 12)
  text(lab$logFC, -log10(pmax(lab$FDR,1e-300)), labels=lab$gene, cex=.62, pos=3)

  e <- res1$enrich
  focus <- c("Antigen processing and presentation",
             "Response To Type II Interferon (GO:0034341)",
             "Inflammatory Response (GO:0006954)",
             "Phagocytosis (GO:0006909)",
             "Positive Regulation Of Phagocytosis (GO:0050766)",
             "Lysosome", "Cytokine-cytokine receptor interaction",
             "Cellular Response To Lipopolysaccharide (GO:0071222)")
  e <- e[e$pathway %in% focus, ]
  if (nrow(e) < 6) e <- head(res1$enrich, 10)
  e <- e[order(-log10(e$FDR)), ]
  labs <- gsub(" \\(GO:.*", "", e$pathway)
  labs <- substr(labs, 1, 52)
  val <- -log10(pmax(e$FDR, 1e-300))
  cc <- ifelse(e$Direction == "Up", "#D73027", "#4575B4")
  barplot(rev(val), names.arg=rev(labs), horiz=TRUE, las=1, col=rev(cc), border=NA,
          cex.names=.62, xlab=expression(-log[10]~FDR), main="B  Ranked pathway enrichment")
  legend("bottomright", legend=c("Higher in cDC1","Lower in cDC1"), fill=c("#D73027","#4575B4"),
         border=NA, bty="n", cex=.75)
  par(op); dev.off()
}
plot_figure("png", file.path(figdir,"Figure8_patient_pseudobulk_DE_and_enrichment.png"), 13.5, 6.2, 400)
plot_figure("pdf", file.path(figdir,"Figure8_patient_pseudobulk_DE_and_enrichment.pdf"), 13.5, 6.2)
plot_figure("svg", file.path(figdir,"Figure8_patient_pseudobulk_DE_and_enrichment.svg"), 13.5, 6.2)

cat("\nPseudobulk analyses completed\n")
print(summary)
cat("\nTop DE genes (primary cDC1 vs cDC2):\n")
print(head(res1$tab, 15))
cat("\nTop pathways:\n")
print(head(res1$enrich, 15))
