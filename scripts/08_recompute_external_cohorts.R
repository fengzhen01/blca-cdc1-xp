args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[[1]] else "."
outdir <- file.path(root, "results", "locked_nonoverlap")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
genes <- c("TAP1","TAP2","TAPBP","SEC22B")

gene_z_score <- function(expr, transform=c("tpm","counts")) {
  transform <- match.arg(transform)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  if (transform == "tpm") logx <- log2(expr + 1) else {
    lib <- colSums(expr)
    logx <- log2(sweep(expr, 2, lib / 1e6, "/") + 1)
  }
  z <- t(scale(t(logx)))
  rowMeans_present <- rownames(z) %in% genes
  stopifnot(sum(rowMeans_present) >= 3)
  colMeans(z[rowMeans_present,,drop=FALSE], na.rm=TRUE)
}

auc_rank <- function(y,x) {
  ok <- is.finite(y)&is.finite(x); y<-y[ok];x<-x[ok]
  n1<-sum(y==1);n0<-sum(y==0)
  (sum(rank(x,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}

boot_auc <- function(y,x,B=5000,seed=20261020) {
  set.seed(seed); pos<-which(y==1);neg<-which(y==0)
  a <- replicate(B,{i<-c(sample(pos,length(pos),TRUE),sample(neg,length(neg),TRUE));auc_rank(y[i],x[i])})
  unname(quantile(a,c(.025,.975)))
}

summarize_response <- function(dat, cohort) {
  keep <- is.finite(dat$response)&is.finite(dat$cross_nonoverlap)
  y<-dat$response[keep];x<-dat$cross_nonoverlap[keep]
  ci<-boot_auc(y,x)
  wt<-wilcox.test(x~y,exact=FALSE)
  data.frame(cohort=cohort,n=length(y),responders=sum(y==1),nonresponders=sum(y==0),
             responder_median=median(x[y==1]),nonresponder_median=median(x[y==0]),
             auc=auc_rank(y,x),auc_low=ci[1],auc_high=ci[2],wilcoxon_p=wt$p.value)
}

# GSE176307
x176 <- read.delim(gzfile(file.path(root,"data","raw","GSE176307_salmon_tpm_gene.matrix.tsv.gz")),
                   check.names=FALSE,row.names=1,quote="",comment.char="")
score176_rs <- gene_z_score(x176,"tpm")
key <- read.csv(file.path(root,"data","raw","GSE176307_sample_key.txt"),
                check.names=FALSE,stringsAsFactors=FALSE,blank.lines.skip=TRUE)
key <- key[!is.na(key[[1]]) & nzchar(key[[1]]),,drop=FALSE]
map <- list()
for(i in seq_len(nrow(key))) {
  ids <- trimws(strsplit(as.character(key[i,2]),",",fixed=TRUE)[[1]])
  for(id in ids) map[[id]] <- as.character(key[i,1])
}
patient_ids <- vapply(names(score176_rs), function(z) if(!is.null(map[[z]])) map[[z]] else z, character(1))
score176 <- tapply(score176_rs, patient_ids, mean)
meta176 <- read.delim(file.path(root,"data","raw","GSE176307_sample_metadata.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
meta176$patient <- sub(".*(BACI[0-9]+).*","\\1",meta176$title)
d176 <- meta176[match(names(score176),meta176$patient),,drop=FALSE]
d176$sample_id <- names(score176)
d176$cross_nonoverlap <- as.numeric(score176)
d176$response <- ifelse(d176$io.response %in% c("CR","PR"),1,ifelse(d176$io.response %in% c("SD","PD"),0,NA))
d176$pfs_event <- ifelse(d176$progressed=="Yes",1,ifelse(d176$progressed=="No",0,NA))
d176$os_event <- ifelse(d176$alive=="No",1,ifelse(d176$alive=="Yes",0,NA))
write.csv(d176,file.path(outdir,"GSE176307_locked_nonoverlap_scores.csv"),row.names=FALSE)

# DUTRENEO
x328 <- read.delim(gzfile(file.path(root,"data","raw","GSE328930_counts.tsv.gz")),
                   check.names=FALSE,row.names=1,quote="",comment.char="")
score328 <- gene_z_score(x328,"counts")
meta328 <- read.delim(file.path(root,"data","raw","GSE328930_metadata.txt"),
                      check.names=FALSE,stringsAsFactors=FALSE)
d328 <- meta328[match(names(score328),meta328$ID),,drop=FALSE]
d328$sample_id <- names(score328)
d328$cross_nonoverlap <- as.numeric(score328)
d328$response <- ifelse(d328$Pathological_response %in% c("Pathological complete response","Pathological partial response"),1,
                        ifelse(d328$Pathological_response=="No response",0,NA))
d328$response[d328$Pathological_response=="No surgery"] <- NA
d328$evaluable <- d328$Group=="HOT-ICI" & is.finite(d328$response)
write.csv(d328,file.path(outdir,"GSE328930_locked_nonoverlap_scores.csv"),row.names=FALSE)

sum176 <- summarize_response(d176,"GSE176307")
sum328 <- summarize_response(d328[d328$evaluable,],"DUTRENEO HOT-ICI")
summary <- rbind(sum176,sum328)
write.csv(summary,file.path(outdir,"locked_nonoverlap_response_summary.csv"),row.names=FALSE)

# Regression/survival for GSE176307 using the identical score.
suppressPackageStartupMessages(library(survival))
d176$cross_z <- as.numeric(scale(d176$cross_nonoverlap))
d176$tmb_z <- as.numeric(scale(as.numeric(d176$tmb)))
g1 <- glm(response~cross_z,data=d176,family=binomial(),na.action=na.omit)
g2 <- glm(response~cross_z+tmb_z,data=d176,family=binomial(),na.action=na.omit)
c1 <- coxph(Surv(as.numeric(pfs),pfs_event)~cross_z,data=d176,na.action=na.omit)
c2 <- coxph(Surv(as.numeric(overall_survival),os_event)~cross_z,data=d176,na.action=na.omit)
extract_glm <- function(fit,label,term="cross_z") {
  sm<-summary(fit)$coef; b<-coef(fit)[term];se<-sm[term,"Std. Error"]
  data.frame(model=label,endpoint="response",term=term,n=nobs(fit),events=sum(model.response(model.frame(fit))==1),
             effect=exp(b),ci_low=exp(b-1.96*se),ci_high=exp(b+1.96*se),p=sm[term,"Pr(>|z|)"])
}
extract_cox <- function(fit,label) {
  sm<-summary(fit);data.frame(model=label,endpoint=label,term="cross_z",n=fit$n,events=fit$nevent,
  effect=sm$coef["cross_z","exp(coef)"],ci_low=sm$conf.int["cross_z","lower .95"],ci_high=sm$conf.int["cross_z","upper .95"],p=sm$coef["cross_z","Pr(>|z|)"])
}
reg <- rbind(extract_glm(g1,"Response, unadjusted"),extract_glm(g2,"Response, TMB-adjusted"),
             extract_cox(c1,"PFS"),extract_cox(c2,"OS"))
write.csv(reg,file.path(outdir,"GSE176307_locked_nonoverlap_models.csv"),row.names=FALSE)
print(summary)
print(reg)
