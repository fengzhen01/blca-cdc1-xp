suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

d <- readRDS(input)
counts <- d$counts
feature <- d$feature
pheno <- d$pheno
stopifnot(identical(colnames(counts), rownames(pheno)))

symbols <- as.character(feature$Symbol)
if (anyDuplicated(symbols)) stop("Duplicated symbols among fixed genes")
rownames(counts) <- symbols

# Locked before opening the IMvigor210 outcome associations.
sets <- list(
  cDC1_identity = c("CLEC9A","XCR1","BATF3","CADM1","IRF8","WDFY4"),
  cross_presentation_nonoverlap = c("TAP1","TAP2","TAPBP","SEC22B"),
  cross_presentation_original = c("WDFY4","TAP1","TAP2","TAPBP","B2M","HLA-A","HLA-B","HLA-C","SEC22B"),
  mature_migratory_DC = c("LAMP3","CCR7","FSCN1","CD80","CD86","CD40","IL12B","CCL19")
)

logcpm <- log2(sweep(counts, 2, d$library_size / 1e6, "/") + 1)
zexpr <- t(scale(t(logcpm)))
zexpr[!is.finite(zexpr)] <- 0

scores <- sapply(sets, function(gs) colMeans(zexpr[intersect(gs, rownames(zexpr)),,drop=FALSE]))
scores <- as.data.frame(scores, check.names=FALSE)
rownames(scores) <- colnames(counts)

dat <- cbind(sample_id=rownames(pheno), pheno, scores)
dat$response4 <- as.character(dat[["Best Confirmed Overall Response"]])
dat$response <- ifelse(dat$response4 %in% c("CR","PR"), 1,
                       ifelse(dat$response4 %in% c("SD","PD"), 0, NA))
dat$tmb <- suppressWarnings(as.numeric(as.character(dat[["FMOne mutation burden per MB"]])))
dat$ecog <- suppressWarnings(as.numeric(as.character(dat[["Baseline ECOG Score"]])))

auc_rank <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(x, ties.method="average")[y == 1]) - n1*(n1+1)/2) / (n1*n0)
}

auc_boot <- function(y, x, B=2000, seed=20261020) {
  ok <- is.finite(y) & is.finite(x); y <- y[ok]; x <- x[ok]
  pos <- which(y==1); neg <- which(y==0)
  set.seed(seed)
  vals <- replicate(B, {
    idx <- c(sample(pos, length(pos), replace=TRUE), sample(neg, length(neg), replace=TRUE))
    auc_rank(y[idx], x[idx])
  })
  unname(quantile(vals, c(.025,.975), na.rm=TRUE))
}

fit_glm <- function(formula, data, term) {
  fit <- glm(formula, data=data, family=binomial(), na.action=na.omit)
  sm <- summary(fit)$coefficients
  b <- coef(fit)[term]; se <- sm[term,"Std. Error"]
  data.frame(
    model=paste(deparse(formula), collapse=" "), term=term, n=nobs(fit), events=sum(model.response(model.frame(fit))==1),
    estimate=unname(b), effect=unname(exp(b)), ci_low=unname(exp(b-1.96*se)), ci_high=unname(exp(b+1.96*se)),
    p=unname(sm[term,"Pr(>|z|)"]), stringsAsFactors=FALSE
  )
}

fit_cox <- function(formula, data, term) {
  fit <- coxph(formula, data=data, na.action=na.omit, ties="efron")
  sm <- summary(fit)
  i <- match(term, rownames(sm$coefficients))
  data.frame(
    model=paste(deparse(formula), collapse=" "), term=term, n=fit$n, events=fit$nevent,
    estimate=sm$coefficients[i,"coef"], effect=sm$coefficients[i,"exp(coef)"],
    ci_low=sm$conf.int[i,"lower .95"], ci_high=sm$conf.int[i,"upper .95"],
    p=sm$coefficients[i,"Pr(>|z|)"], stringsAsFactors=FALSE
  )
}

score_names <- names(sets)
summary_rows <- list()
for (s in score_names) {
  idx <- is.finite(dat$response) & is.finite(dat[[s]])
  wt <- wilcox.test(dat[[s]][idx] ~ dat$response[idx], exact=FALSE, conf.int=FALSE)
  ci <- auc_boot(dat$response[idx], dat[[s]][idx])
  summary_rows[[s]] <- data.frame(
    score=s, n=sum(idx), responders=sum(dat$response[idx]==1), nonresponders=sum(dat$response[idx]==0),
    responder_mean=mean(dat[[s]][idx & dat$response==1]),
    nonresponder_mean=mean(dat[[s]][idx & dat$response==0]),
    auc=auc_rank(dat$response[idx], dat[[s]][idx]), auc_ci_low=ci[1], auc_ci_high=ci[2],
    wilcoxon_p=wt$p.value, stringsAsFactors=FALSE
  )
}
score_summary <- do.call(rbind, summary_rows)
score_summary$wilcoxon_fdr <- p.adjust(score_summary$wilcoxon_p, method="BH")

# Convert scores to per-SD effects within all available samples.
for (s in score_names) dat[[paste0(s,"_z")]] <- as.numeric(scale(dat[[s]]))
dat$log1p_tmb_z <- as.numeric(scale(log1p(dat$tmb)))

logit_results <- list()
for (s in score_names) {
  z <- paste0(s,"_z")
  logit_results[[paste0(s,"_unadj")]] <- fit_glm(as.formula(paste("response ~", z)), dat, z)
  cc <- dat[is.finite(dat$response) & is.finite(dat$log1p_tmb_z),,drop=FALSE]
  logit_results[[paste0(s,"_tmb_cc_unadj")]] <- fit_glm(as.formula(paste("response ~", z)), cc, z)
  logit_results[[paste0(s,"_tmb")]] <- fit_glm(as.formula(paste("response ~", z, "+ log1p_tmb_z")), dat, z)
}
# Primary incremental model: antigen-processing score beyond cDC1 abundance and TMB.
logit_results[["cross_incremental"]] <- fit_glm(
  response ~ cross_presentation_nonoverlap_z + cDC1_identity_z + log1p_tmb_z,
  dat, "cross_presentation_nonoverlap_z"
)
logit_results <- do.call(rbind, logit_results)

cox_results <- list()
for (s in score_names) {
  z <- paste0(s,"_z")
  cox_results[[paste0(s,"_unadj")]] <- fit_cox(as.formula(paste("Surv(os, censOS) ~", z)), dat, z)
  cc <- dat[is.finite(dat$os) & is.finite(dat$censOS) & is.finite(dat$log1p_tmb_z),,drop=FALSE]
  cox_results[[paste0(s,"_tmb_cc_unadj")]] <- fit_cox(as.formula(paste("Surv(os, censOS) ~", z)), cc, z)
  cox_results[[paste0(s,"_tmb")]] <- fit_cox(as.formula(paste("Surv(os, censOS) ~", z, "+ log1p_tmb_z")), dat, z)
}
cox_results[["cross_incremental"]] <- fit_cox(
  Surv(os, censOS) ~ cross_presentation_nonoverlap_z + cDC1_identity_z + log1p_tmb_z,
  dat, "cross_presentation_nonoverlap_z"
)
cox_results <- do.call(rbind, cox_results)

corr <- cor(dat[,score_names], use="pairwise.complete.obs", method="spearman")

write.csv(dat, file.path(outdir,"IMvigor210_locked_sample_scores.csv"), row.names=FALSE)
write.csv(score_summary, file.path(outdir,"IMvigor210_response_score_summary.csv"), row.names=FALSE)
write.csv(logit_results, file.path(outdir,"IMvigor210_logistic_results.csv"), row.names=FALSE)
write.csv(cox_results, file.path(outdir,"IMvigor210_OS_cox_results.csv"), row.names=FALSE)
write.csv(corr, file.path(outdir,"IMvigor210_score_spearman_correlations.csv"), row.names=TRUE)

cat("RESPONSE SCORE SUMMARY\n")
print(score_summary, row.names=FALSE)
cat("\nLOGISTIC RESULTS\n")
print(logit_results, row.names=FALSE)
cat("\nOS COX RESULTS\n")
print(cox_results, row.names=FALSE)
cat("\nSCORE SPEARMAN CORRELATIONS\n")
print(round(corr,3))
