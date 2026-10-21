suppressPackageStartupMessages(library(survival))

root <- "."
input <- file.path(root,"results","IMvigor210_validation","IMvigor210_fixed_data.rds")
outdir <- file.path(root,"results","immunotherapy")
figdir <- file.path(root,"figures","immunotherapy")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE); dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
d <- readRDS(input); counts <- d$counts; pheno <- d$pheno; feature <- d$feature
rownames(counts) <- as.character(feature$Symbol)

sets <- list(
  cDC1_identity=c("CLEC9A","XCR1","BATF3","CADM1","IRF8","WDFY4"),
  cross_presentation=c("TAP1","TAP2","TAPBP","SEC22B"),
  T_cell_inflamed=c("CCL5","CD27","CD274","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1","HLA-DRB1","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT","GZMA","PRF1"),
  cytotoxic_T=c("CD8A","CD8B","GZMA","GZMB","PRF1","IFNG","TBX21","NKG7","CCL5"),
  IFN_response=c("STAT1","ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","CXCL9","CXCL10","IDO1"),
  fibroblast_TGFb=c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","ACTA2","TAGLN","FAP","TGFBI"),
  suppressive_myeloid=c("SPP1","APOE","C1QA","C1QB","C1QC","IL10","TGFB1")
)
logcpm <- log2(sweep(counts,2,d$library_size/1e6,"/")+1)
zexpr <- t(scale(t(logcpm))); zexpr[!is.finite(zexpr)] <- 0
scores <- sapply(sets,function(gs) colMeans(zexpr[intersect(gs,rownames(zexpr)),,drop=FALSE]))
scores <- as.data.frame(scores); rownames(scores)<-colnames(counts)
dat <- cbind(sample_id=rownames(pheno),pheno,scores)
dat$response4 <- as.character(dat[["Best Confirmed Overall Response"]])
dat$response <- ifelse(dat$response4 %in% c("CR","PR"),1,ifelse(dat$response4 %in% c("SD","PD"),0,NA))
dat$tmb <- suppressWarnings(as.numeric(as.character(dat[["FMOne mutation burden per MB"]])))
for(s in names(sets)) dat[[paste0(s,"_z")]] <- as.numeric(scale(dat[[s]]))
dat$log1p_tmb_z <- as.numeric(scale(log1p(dat$tmb)))

effect_glm <- function(f,data,term,label) {
  fit<-glm(f,data=data,family=binomial(),na.action=na.omit); sm<-summary(fit)$coef
  b<-coef(fit)[term]; se<-sm[term,"Std. Error"]
  data.frame(label=label,type="Response OR",n=nobs(fit),events=sum(model.response(model.frame(fit))==1),
             effect=exp(b),lo=exp(b-1.96*se),hi=exp(b+1.96*se),p=sm[term,"Pr(>|z|)"])
}
effect_cox <- function(f,data,term,label) {
  fit<-coxph(f,data=data,na.action=na.omit,ties="efron"); sm<-summary(fit); i<-match(term,rownames(sm$coef))
  data.frame(label=label,type="Overall survival HR",n=fit$n,events=fit$nevent,effect=sm$coef[i,"exp(coef)"],
             lo=sm$conf.int[i,"lower .95"],hi=sm$conf.int[i,"upper .95"],p=sm$coef[i,"Pr(>|z|)"])
}
rows <- list(
  effect_glm(response~cross_presentation_z,dat,"cross_presentation_z","Cross-presentation (unadjusted)"),
  effect_glm(response~cross_presentation_z+log1p_tmb_z,dat,"cross_presentation_z","Cross-presentation + TMB"),
  effect_glm(response~cross_presentation_z+log1p_tmb_z+cDC1_identity_z,dat,"cross_presentation_z","Cross-presentation + TMB + cDC1"),
  effect_glm(response~cross_presentation_z+log1p_tmb_z+T_cell_inflamed_z,dat,"cross_presentation_z","Cross-presentation + TMB + T-cell inflamed"),
  effect_cox(Surv(os,censOS)~cross_presentation_z,dat,"cross_presentation_z","Cross-presentation (unadjusted)"),
  effect_cox(Surv(os,censOS)~cross_presentation_z+log1p_tmb_z,dat,"cross_presentation_z","Cross-presentation + TMB"),
  effect_cox(Surv(os,censOS)~cross_presentation_z+log1p_tmb_z+cDC1_identity_z,dat,"cross_presentation_z","Cross-presentation + TMB + cDC1"),
  effect_cox(Surv(os,censOS)~cross_presentation_z+log1p_tmb_z+T_cell_inflamed_z,dat,"cross_presentation_z","Cross-presentation + TMB + T-cell inflamed")
)
effects <- do.call(rbind,rows); effects$FDR <- p.adjust(effects$p,"BH")
write.csv(effects,file.path(outdir,"IMvigor210_adjusted_effect_models.csv"),row.names=FALSE)

cors <- cor(dat[,names(sets)],method="spearman",use="pairwise.complete.obs")
write.csv(cors,file.path(outdir,"IMvigor210_immune_score_correlations.csv"))

group_tests <- function(var) {
  z<-dat[is.finite(dat$cross_presentation) & !is.na(dat[[var]]) & dat[[var]]!="",]
  if(length(unique(z[[var]]))<2) return(NULL)
  kw<-kruskal.test(z$cross_presentation~factor(z[[var]]))
  data.frame(variable=var,n=nrow(z),groups=length(unique(z[[var]])),kruskal_p=kw$p.value)
}
gt <- do.call(rbind,lapply(c("Immune phenotype","IC Level","TC Level","TCGA Subtype","Lund2"),group_tests))
gt$FDR<-p.adjust(gt$kruskal_p,"BH"); write.csv(gt,file.path(outdir,"IMvigor210_clinical_group_tests.csv"),row.names=FALSE)

# Interaction models are explicitly exploratory.
int_glm<-glm(response~cross_presentation_z*log1p_tmb_z,data=dat,family=binomial(),na.action=na.omit)
int_cox<-coxph(Surv(os,censOS)~cross_presentation_z*log1p_tmb_z,data=dat,na.action=na.omit)
ints<-rbind(data.frame(outcome="Response",n=nobs(int_glm),effect=exp(coef(int_glm)["cross_presentation_z:log1p_tmb_z"]),p=summary(int_glm)$coef["cross_presentation_z:log1p_tmb_z","Pr(>|z|)"]),
            data.frame(outcome="Overall survival",n=int_cox$n,effect=exp(coef(int_cox)["cross_presentation_z:log1p_tmb_z"]),p=summary(int_cox)$coef["cross_presentation_z:log1p_tmb_z","Pr(>|z|)"]))
write.csv(ints,file.path(outdir,"IMvigor210_cross_TMB_interactions_exploratory.csv"),row.names=FALSE)
write.csv(dat,file.path(outdir,"IMvigor210_extended_scores.csv"),row.names=FALSE)

plot_context <- function(dev,path,w=15,h=5.4,res=400) {
  if(dev=="png") png(path,w*res,h*res,res=res); if(dev=="pdf") pdf(path,w,h,useDingbats=FALSE); if(dev=="svg") svg(path,w,h)
  par(mfrow=c(1,3),mar=c(6,7,3,2),family="sans")
  image(seq_len(nrow(cors)),seq_len(ncol(cors)),t(cors[nrow(cors):1,]),col=colorRampPalette(c("#2166AC","white","#B2182B"))(101),zlim=c(-1,1),axes=FALSE,xlab="",ylab="",main="A  Spearman correlation of immune scores")
  axis(1,at=seq_len(ncol(cors)),labels=colnames(cors),las=2,cex.axis=.65); axis(2,at=seq_len(nrow(cors)),labels=rev(rownames(cors)),las=2,cex.axis=.65)
  for(i in seq_len(nrow(cors))) for(j in seq_len(ncol(cors))) text(j,nrow(cors)-i+1,sprintf("%.2f",cors[i,j]),cex=.55)
  z<-dat[is.finite(dat$response),]; gr <- factor(z$response,levels=c(0,1),labels=c("Non-response","Response")); boxplot(z$cross_presentation~gr,col=c("#A6CEE3","#FB9A99"),outline=FALSE,ylab="Cross-presentation score",main="B  Objective response")
  stripchart(z$cross_presentation~gr,vertical=TRUE,method="jitter",pch=16,col=adjustcolor("black",.35),add=TRUE)
  p<-wilcox.test(z$cross_presentation~gr,exact=FALSE)$p.value; mtext(sprintf("Wilcoxon P = %.3g",p),side=3,line=.2,cex=.8)
  z<-dat[is.finite(dat$cross_presentation)&!is.na(dat[["Immune phenotype"]])&dat[["Immune phenotype"]]!="",]; z[["Immune phenotype"]]<-factor(z[["Immune phenotype"]],levels=c("desert","excluded","inflamed"))
  boxplot(z$cross_presentation~z[["Immune phenotype"]],col=c("#FDD49E","#9ECAE1","#74C476"),outline=FALSE,ylab="Cross-presentation score",xlab="Immune phenotype",main="C  Tumor immune phenotype")
  stripchart(z$cross_presentation~z[["Immune phenotype"]],vertical=TRUE,method="jitter",pch=16,col=adjustcolor("black",.3),add=TRUE)
  p<-kruskal.test(z$cross_presentation~z[["Immune phenotype"]])$p.value; mtext(sprintf("Kruskal-Wallis P = %.3g",p),side=3,line=.2,cex=.8)
  dev.off()
}
plot_context("png",file.path(figdir,"Figure9_IMvigor210_immune_context.png")); plot_context("pdf",file.path(figdir,"Figure9_IMvigor210_immune_context.pdf")); plot_context("svg",file.path(figdir,"Figure9_IMvigor210_immune_context.svg"))

plot_forest <- function(dev,path,w=12,h=5.8,res=400) {
  if(dev=="png") png(path,w*res,h*res,res=res); if(dev=="pdf") pdf(path,w,h,useDingbats=FALSE); if(dev=="svg") svg(path,w,h)
  par(mar=c(5,19,3,3),family="sans"); e<-effects; y<-rev(seq_len(nrow(e))); col<-ifelse(e$type=="Response OR","#D95F02","#1B9E77")
  plot(e$effect,y,log="x",xlim=c(.45,3.5),pch=16,col=col,yaxt="n",ylab="",xlab="Effect per 1-SD increase in cross-presentation score",main="Adjusted associations in IMvigor210")
  segments(e$lo,y,e$hi,y,col=col,lwd=2); abline(v=1,lty=2,col="grey40"); axis(2,at=y,labels=paste0(e$label," [",ifelse(e$type=="Response OR","OR","HR"),"]"),las=2,cex.axis=.68)
  text(3.4,y,sprintf("P=%.3g",e$p),pos=2,cex=.72); legend("bottomleft",legend=c("Response OR","Overall survival HR"),col=c("#D95F02","#1B9E77"),pch=16,bty="n")
  dev.off()
}
plot_forest("png",file.path(figdir,"Figure10_IMvigor210_adjusted_forest.png")); plot_forest("pdf",file.path(figdir,"Figure10_IMvigor210_adjusted_forest.pdf")); plot_forest("svg",file.path(figdir,"Figure10_IMvigor210_adjusted_forest.svg"))

cat("\nAdjusted models\n"); print(effects,row.names=FALSE)
cat("\nClinical group tests\n"); print(gt,row.names=FALSE)
cat("\nInteractions (exploratory)\n"); print(ints,row.names=FALSE)
