suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else "."
figdir <- file.path(root,"figures","overview")
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)

options(stringsAsFactors=FALSE)
set.seed(20261020)

COL <- c(blue="#2F6B9A",orange="#D97732",teal="#258D7A",purple="#7A6AA6",
         red="#B94E5E",gold="#D5A021",gray="#66727D",light="#EAF0F4",dark="#24323D")
GROUP_COL <- c("Non-responder"=COL["gray"],"Responder"=COL["orange"],
               "Low"=COL["blue"],"High"=COL["orange"])

export_figure <- function(stem,width,height,draw) {
  png(file.path(figdir,paste0(stem,".png")),width=width,height=height,units="in",res=400,type="cairo",bg="white")
  draw(); dev.off()
  pdf(file.path(figdir,paste0(stem,".pdf")),width=width,height=height,useDingbats=FALSE,family="sans",paper="special")
  draw(); dev.off()
  svg(file.path(figdir,paste0(stem,".svg")),width=width,height=height,family="sans",bg="white")
  draw(); dev.off()
}

panel_label <- function(lbl,x=-0.13,y=1.08,cex=1.35) {
  usr<-par("usr")
  xx<-usr[1]+x*diff(usr[1:2]); yy<-usr[3]+y*diff(usr[3:4])
  if(par("xlog")) xx<-10^xx
  if(par("ylog")) yy<-10^yy
  text(xx,yy,lbl,font=2,cex=cex,xpd=NA)
}

p_fmt <- function(p) ifelse(p<0.001,"P < 0.001",paste0("P = ",formatC(p,format="f",digits=3)))

draw_box <- function(x1,y1,x2,y2,fill,label,sub=NULL,cex=1.0) {
  rect(x1,y1,x2,y2,col=fill,border=COL["dark"],lwd=1.2)
  text((x1+x2)/2,(y1+y2)/2+ifelse(is.null(sub),0,0.035),label,font=2,cex=cex)
  if(!is.null(sub)) text((x1+x2)/2,(y1+y2)/2-0.055,sub,cex=0.78,col=COL["gray"])
}

draw_arrow <- function(x1,y1,x2,y2) arrows(x1,y1,x2,y2,length=.08,lwd=1.5,col=COL["gray"])

score_boxplot <- function(dat,score_col,response_col="response",title="",ylab="Score",annotate_auc=NULL) {
  keep<-is.finite(dat[[response_col]])&is.finite(dat[[score_col]])
  y<-dat[[response_col]][keep];x<-dat[[score_col]][keep]
  lev<-c(0,1); labs<-c("Non-responder","Responder")
  boxplot(split(x,factor(y,levels=lev)),names=labs,col=c("#DCE2E7","#F4C7A7"),border=COL["dark"],
          outline=FALSE,ylab=ylab,main=title,las=1,cex.axis=.85,cex.main=1.0)
  for(i in 1:2) {
    v<-x[y==lev[i]]
    points(jitter(rep(i,length(v)),amount=.13),v,pch=21,bg=adjustcolor(c(COL["gray"],COL["orange"])[i],.48),
           col=adjustcolor("white",.4),cex=.55)
  }
  p<-wilcox.test(x~y,exact=FALSE)$p.value
  mtext(sprintf("n = %d vs %d; %s",sum(y==0),sum(y==1),p_fmt(p)),side=3,line=.05,cex=.74,col=COL["gray"])
  if(!is.null(annotate_auc)) legend("topleft",legend=sprintf("AUC %.3f",annotate_auc),bty="n",cex=.78,text.col=COL["blue"])
}

auc_rank <- function(y,x) {
  ok<-is.finite(y)&is.finite(x);y<-y[ok];x<-x[ok];n1<-sum(y==1);n0<-sum(y==0)
  (sum(rank(x,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}

roc_xy <- function(y,x) {
  ok<-is.finite(y)&is.finite(x);y<-y[ok];x<-x[ok]
  o<-order(x,decreasing=TRUE);y<-y[o]
  tp<-c(0,cumsum(y==1)/sum(y==1));fp<-c(0,cumsum(y==0)/sum(y==0))
  list(fpr=fp,tpr=tp)
}

plot_forest <- function(labels,effect,low,high,p=NULL,ref=1,xlim=NULL,xlab="Effect per 1-SD increase",cols=NULL,cex=.82) {
  n<-length(labels); y<-rev(seq_len(n))
  if(is.null(xlim)) xlim<-range(c(low,high),finite=TRUE)*c(.9,1.1)
  plot(NA,xlim=xlim,ylim=c(.4,n+.6),log="x",yaxt="n",xlab=xlab,ylab="",bty="n")
  abline(v=ref,lty=2,col="#9AA5AD")
  axis(2,at=y,labels=labels,las=1,tick=FALSE,cex.axis=cex)
  if(is.null(cols)) cols<-rep(COL["blue"],n)
  segments(low,y,high,y,lwd=2,col=cols)
  points(effect,y,pch=21,bg=cols,col="white",cex=1.25)
  if(!is.null(p)) {
    txt<-ifelse(p<.001,"P<0.001",sprintf("P=%.3f",p))
    text(rep(xlim[2],n),y,txt,pos=2,cex=.7,col=COL["gray"],xpd=NA)
  }
}

plot_forest_inside <- function(labels,effect,low,high,p=NULL,ref=1,xlim,xlab,cols=NULL,cex=.78) {
  n<-length(labels); y<-rev(seq_len(n))
  plot(NA,xlim=xlim,ylim=c(.4,n+.6),log="x",yaxt="n",xlab=xlab,ylab="",bty="n")
  abline(v=ref,lty=2,col="#9AA5AD")
  if(is.null(cols)) cols<-rep(COL["blue"],n)
  segments(low,y,high,y,lwd=2,col=cols)
  points(effect,y,pch=21,bg=cols,col="white",cex=1.25)
  text(rep(xlim[1]*1.07,n),y,labels,adj=0,cex=cex,col=COL["dark"])
  if(!is.null(p)) {
    txt<-ifelse(p<.001,"P<0.001",sprintf("P=%.3f",p))
    text(rep(xlim[2]/1.03,n),y,txt,adj=1,cex=.70,col=COL["gray"])
  }
}

km_panel <- function(time,event,score,title,xlab="Time",max_time=NULL) {
  ok<-is.finite(time)&is.finite(event)&is.finite(score)
  d<-data.frame(time=time[ok],event=event[ok],score=score[ok])
  d$group<-factor(ifelse(d$score>median(d$score),"High","Low"),levels=c("Low","High"))
  fit<-survfit(Surv(time,event)~group,data=d)
  if(is.null(max_time)) max_time<-max(d$time,na.rm=TRUE)
  plot(fit,col=c(COL["blue"],COL["orange"]),lwd=2,mark.time=TRUE,conf.int=FALSE,
       xlab=xlab,ylab="Survival probability",main=title,xlim=c(0,max_time),bty="l",las=1)
  lr<-survdiff(Surv(time,event)~group,data=d)
  p<-pchisq(lr$chisq,df=1,lower.tail=FALSE)
  legend("bottomleft",legend=c(sprintf("Low (n=%d)",sum(d$group=="Low")),sprintf("High (n=%d)",sum(d$group=="High"))),
         col=c(COL["blue"],COL["orange"]),lwd=2,bty="n",cex=.76)
  legend("topright",legend=p_fmt(p),bty="n",cex=.76)
  mtext("Median split for visualization only",side=3,line=.1,cex=.65,col=COL["gray"])
}

# Load data ------------------------------------------------------------------
locked_sum<-read.csv(file.path(root,"results","locked_nonoverlap","locked_nonoverlap_response_summary.csv"),check.names=FALSE)
d176<-read.csv(file.path(root,"results","locked_nonoverlap","GSE176307_locked_nonoverlap_scores.csv"),check.names=FALSE)
d328<-read.csv(file.path(root,"results","locked_nonoverlap","GSE328930_locked_nonoverlap_scores.csv"),check.names=FALSE)
d210<-read.csv(file.path(root,"results","immunotherapy","IMvigor210_locked_sample_scores.csv"),check.names=FALSE)
sum210<-read.csv(file.path(root,"results","immunotherapy","IMvigor210_response_score_summary.csv"),check.names=FALSE)
logit210<-read.csv(file.path(root,"results","immunotherapy","IMvigor210_logistic_results.csv"),check.names=FALSE)
cox210<-read.csv(file.path(root,"results","immunotherapy","IMvigor210_OS_cox_results.csv"),check.names=FALSE)
corr210<-as.matrix(read.csv(file.path(root,"results","immunotherapy","IMvigor210_score_spearman_correlations.csv"),row.names=1,check.names=FALSE))
scscores<-read.csv(file.path(root,"results","single_cell_targeted","GSE293189_patient_group_scores.csv"),check.names=FALSE)
sctests<-read.csv(file.path(root,"results","single_cell_targeted","GSE293189_patient_paired_tests.csv"),check.names=FALSE)
sccounts<-read.csv(file.path(root,"results","single_cell_targeted","GSE293189_strict_group_counts.csv"),check.names=FALSE)
scqc<-read.csv(file.path(root,"results","single_cell_targeted","GSE293189_qc_summary.csv"),check.names=FALSE)
reg176<-read.csv(file.path(root,"results","locked_nonoverlap","GSE176307_locked_nonoverlap_models.csv"),check.names=FALSE)
# Some source CSV files contain a UTF-8 BOM in the first header.
names(sccounts)[1] <- "comparison_group"
names(scscores)[1] <- "patient"
names(sctests)[1] <- "metric"

# Standardized presentation labels.
group_lab<-c(cDC1="cDC1",cDC2="cDC2",pDC="pDC",mature_migratory_DC="Mature DC",
             nonDC_myeloid="Other myeloid",epithelial="Epithelial")

# Figure 1: study design and sample flow -------------------------------------
draw_fig1<-function() {
  par(mfrow=c(1,2),mar=c(1.5,1.5,2,1.2),family="sans")
  plot.new();plot.window(c(0,1),c(0,1));panel_label("A",x=-.04,y=1.02)
  title("Study design",cex.main=1.1)
  draw_box(.08,.78,.92,.94,"#D9EAF4","Locked hypothesis","cDC1 function may outperform cDC1 abundance",.98)
  draw_arrow(.5,.78,.5,.69)
  draw_box(.05,.49,.46,.68,"#E1F0EA","Cellular localization","GSE293189 single-cell RNA-seq",.95)
  draw_box(.54,.49,.95,.68,"#FCE8DA","Clinical validation","Three real ICI cohorts",.95)
  draw_arrow(.26,.49,.42,.37);draw_arrow(.74,.49,.58,.37)
  draw_box(.25,.19,.75,.37,"#E8E2F2","Locked four-gene program","TAP1 / TAP2 / TAPBP / SEC22B",.95)
  draw_arrow(.5,.19,.5,.11)
  draw_box(.18,.00,.82,.105,"#EAF0F4","Response, survival and robustness analyses",NULL,.74)

  plot.new();plot.window(c(0,1),c(0,1));panel_label("B",x=-.04,y=1.02)
  title("Cohorts and analysis populations",cex.main=1.1)
  rows<-data.frame(
    cohort=c("GSE293189","GSE176307","DUTRENEO","IMvigor210"),
    role=c("Cell localization","ICI discovery","ICI direction check","Independent ICI validation"),
    flow=c("56,192 cells -> 26,468 QC -> 65 cDC1","90 GEO samples -> 88 response -> 87 TMB complete","66 tumors -> 28 HOT-ICI -> 26 evaluable","348 tumors -> 298 response -> 234 TMB complete"),
    stringsAsFactors=FALSE)
  yy<-c(.82,.61,.40,.19)
  fills<-c("#E1F0EA","#D9EAF4","#F8EED4","#FCE8DA")
  for(i in 1:4) {
    rect(.03,yy[i]-.085,.97,yy[i]+.085,col=fills[i],border=NA)
    text(.06,yy[i]+.03,rows$cohort[i],adj=0,font=2,cex=.92)
    text(.94,yy[i]+.03,rows$role[i],adj=1,cex=.72,col=COL["gray"])
    text(.06,yy[i]-.035,rows$flow[i],adj=0,cex=.76,col=COL["dark"])
  }
}
export_figure("Figure1_study_design_and_cohorts",11.5,5.3,draw_fig1)

# Figure 2: single-cell localization ----------------------------------------
draw_fig2<-function() {
  layout(matrix(1:4,2,2,byrow=TRUE));par(mar=c(4.4,4.6,2.7,1.0),family="sans")
  # A QC totals
  vals<-c(sum(scqc$raw_cells),sum(scqc$pass_cells))
  bp<-barplot(vals,names.arg=c("Raw cells","QC-passed"),col=c("#C8D5DE",COL["teal"]),border=NA,
              ylab="Number of cells",main="GSE293189 quality control",las=1,ylim=c(0,max(vals)*1.18))
  text(bp,vals,format(vals,big.mark=","),pos=3,cex=.82);panel_label("A")
  # B group counts
  ord<-c("cDC1","cDC2","pDC","mature_migratory_DC","nonDC_myeloid","epithelial")
  z<-sccounts[match(ord,sccounts$comparison_group),]
  bp<-barplot(log10(z$cells),names.arg=group_lab[z$comparison_group],las=2,
              col=c(COL["orange"],COL["blue"],COL["purple"],COL["teal"],"#8B6F47","#A6AEB5"),border=NA,
              ylab=expression(log[10](cell~count)),main="Strict marker-defined groups",ylim=c(0,4.25),cex.names=.73)
  text(bp,log10(z$cells),z$cells,pos=3,cex=.70);panel_label("B")
  # C patient-level score distributions
  ord<-c("cDC1","cDC2","pDC","mature_migratory_DC","nonDC_myeloid","epithelial")
  spl<-lapply(ord,function(g) scscores$raw_cross_presentation_nonoverlap_mean[scscores$comparison_group==g])
  boxplot(spl,names=group_lab[ord],las=2,col=c("#F4C7A7","#C9DDED","#DDD6ED","#CBE8E2","#E1D4C3","#E1E5E8"),
          border=COL["dark"],outline=FALSE,ylab="Patient-level mean score",main="Cross-presentation program",cex.names=.73)
  for(i in seq_along(spl)) points(jitter(rep(i,length(spl[[i]])),amount=.12),spl[[i]],pch=21,bg=adjustcolor(COL["dark"],.35),col="white",cex=.65)
  panel_label("C")
  # D paired effect summary
  tt<-sctests[sctests$metric=="raw_cross_presentation_nonoverlap_mean",]
  tt<-tt[match(c("cDC2","pDC","mature_migratory_DC","nonDC_myeloid","epithelial"),tt$comparator),]
  yy<-rev(seq_len(nrow(tt)))
  plot(tt$median_paired_difference,yy,xlim=c(-.02,max(tt$median_paired_difference)*1.6),ylim=c(.5,nrow(tt)+.5),
       yaxt="n",xlab="Median paired difference (cDC1 - comparator)",ylab="",pch=21,bg=COL["orange"],col="white",cex=1.35,
       main="Within-patient paired comparisons")
  abline(v=0,lty=2,col="#9AA5AD");axis(2,at=yy,labels=group_lab[tt$comparator],las=1,tick=FALSE,cex.axis=.8)
  segments(0,yy,tt$median_paired_difference,yy,col=adjustcolor(COL["orange"],.65),lwd=2)
  text(max(tt$median_paired_difference)*1.52,yy,sprintf("FDR %.3f",tt$fdr_all_tests),pos=2,cex=.70,col=COL["gray"])
  panel_label("D")
}
export_figure("Figure2_single_cell_localization",11.5,8.5,draw_fig2)

# Figure 3: cross-cohort response validation --------------------------------
imrow<-sum210[sum210$score=="cross_presentation_nonoverlap",]
aucdf<-rbind(
  data.frame(cohort=locked_sum$cohort,auc=locked_sum$auc,low=locked_sum$auc_low,high=locked_sum$auc_high,p=locked_sum$wilcoxon_p),
  data.frame(cohort="IMvigor210",auc=imrow$auc,low=imrow$auc_ci_low,high=imrow$auc_ci_high,p=imrow$wilcoxon_p)
)
write.csv(aucdf,file.path(figdir,"Figure3_cross_cohort_AUC_source_data.csv"),row.names=FALSE)

draw_fig3<-function() {
  layout(matrix(1:4,2,2,byrow=TRUE));par(mar=c(4.1,4.4,2.8,1),family="sans")
  score_boxplot(d176,"cross_nonoverlap",title="GSE176307",ylab="Locked cross-presentation score",annotate_auc=aucdf$auc[1]);panel_label("A")
  score_boxplot(d328[d328$evaluable,],"cross_nonoverlap",title="DUTRENEO HOT-ICI",ylab="Locked cross-presentation score",annotate_auc=aucdf$auc[2]);panel_label("B")
  score_boxplot(d210,"cross_presentation_nonoverlap",title="IMvigor210",ylab="Locked cross-presentation score",annotate_auc=aucdf$auc[3]);panel_label("C")
  yy<-rev(seq_len(nrow(aucdf)))
  plot(aucdf$auc,yy,xlim=c(.2,.86),ylim=c(.5,nrow(aucdf)+.5),yaxt="n",xlab="AUC for response",ylab="",pch=21,
       bg=c(COL["blue"],COL["gold"],COL["orange"]),col="white",cex=1.45,main="Cross-cohort effect consistency")
  abline(v=.5,lty=2,col="#9AA5AD");segments(aucdf$low,yy,aucdf$high,yy,lwd=2,col=c(COL["blue"],COL["gold"],COL["orange"]))
  axis(2,at=yy,labels=aucdf$cohort,las=1,tick=FALSE,cex.axis=.85)
  text(.85,yy,sprintf("AUC %.3f; %s",aucdf$auc,ifelse(aucdf$p<.001,"P<0.001",sprintf("P=%.3f",aucdf$p))),pos=2,cex=.72,col=COL["gray"])
  panel_label("D")
}
export_figure("Figure3_cross_cohort_response_validation",11.5,8.5,draw_fig3)

# Figure 4: IMvigor specificity and robustness ------------------------------
score_cols<-c("cDC1_identity","cross_presentation_nonoverlap","cross_presentation_original","mature_migratory_DC")
score_names<-c("cDC1 identity","Locked cross-presentation","Extended cross-presentation","Mature migratory DC")
score_colors<-c(COL["blue"],COL["orange"],COL["purple"],COL["teal"])
sm<-sum210[match(score_cols,sum210$score),]

draw_fig4<-function() {
  layout(matrix(1:4,2,2,byrow=TRUE));par(mar=c(4.3,4.5,2.8,1),family="sans")
  # A ROC
  plot(c(0,1),c(0,1),type="n",xlab="1 - specificity",ylab="Sensitivity",main="IMvigor210 ROC curves",asp=1,bty="l")
  abline(0,1,lty=2,col="#AAB2B8")
  for(i in seq_along(score_cols)) {r<-roc_xy(d210$response,d210[[score_cols[i]]]);lines(r$fpr,r$tpr,lwd=2,col=score_colors[i])}
  legend("bottomright",legend=sprintf("%s (AUC %.3f)",score_names,sum210$auc[match(score_cols,sum210$score)]),
         col=score_colors,lwd=2,bty="n",cex=.68);panel_label("A")
  # B AUC forest
  yy<-rev(seq_len(nrow(sm)))
  plot(sm$auc,yy,xlim=c(.35,.75),ylim=c(.5,nrow(sm)+.5),yaxt="n",xlab="AUC for response",ylab="",pch=21,
       bg=score_colors,col="white",cex=1.35,main="Program specificity")
  abline(v=.5,lty=2,col="#9AA5AD");segments(sm$auc_ci_low,yy,sm$auc_ci_high,yy,lwd=2,col=score_colors)
  axis(2,at=yy,labels=score_names,las=1,tick=FALSE,cex.axis=.74)
  text(.745,yy,sprintf("FDR %.3f",sm$wilcoxon_fdr),pos=2,cex=.68,col=COL["gray"]);panel_label("B")
  # C logistic model robustness
  par(mar=c(4.3,4.5,2.8,1.0))
  idx<-c(4,5,6,nrow(logit210))
  dd<-logit210[idx,]
  labs<-c("All evaluable\n(n=298)","Complete cases\n(n=234)","+ TMB\n(n=234)","+ TMB + cDC1\n(n=234)")
  plot_forest_inside(labs,dd$effect,dd$ci_low,dd$ci_high,dd$p,xlim=c(.25,2.25),xlab="Odds ratio for response per 1-SD")
  title("Robustness to case selection and covariates",cex.main=1.0);panel_label("C")
  # D correlation heatmap
  par(mar=c(4.3,4.5,2.8,1.0))
  labs2<-c("cDC1","Locked XP","Extended XP","Mature DC")
  image(seq_len(4),seq_len(4),t(corr210[4:1,]),zlim=c(-1,1),col=colorRampPalette(c("#3B72A1","white","#C95A4A"))(101),
        axes=FALSE,xlab="",ylab="",main="Spearman correlation between programs")
  axis(1,1:4,labs2,las=2,cex.axis=.75);axis(2,1:4,rev(labs2),las=1,cex.axis=.75,tick=FALSE)
  for(i in 1:4) for(j in 1:4) text(i,5-j,sprintf("%.2f",corr210[j,i]),cex=.8,col=ifelse(abs(corr210[j,i])>.65,"white",COL["dark"]))
  panel_label("D")
}
export_figure("Figure4_IMvigor_specificity_and_robustness",11.5,8.7,draw_fig4)

# Figure 5: survival associations -------------------------------------------
d176$pfs_event_num<-ifelse(d176$progressed=="Yes",1,ifelse(d176$progressed=="No",0,NA))
d176$os_event_num<-ifelse(d176$alive=="No",1,ifelse(d176$alive=="Yes",0,NA))

ddsurv<-rbind(
  data.frame(label="GSE176307 PFS",effect=reg176$effect[reg176$model=="PFS"],low=reg176$ci_low[reg176$model=="PFS"],high=reg176$ci_high[reg176$model=="PFS"],p=reg176$p[reg176$model=="PFS"]),
  data.frame(label="GSE176307 OS",effect=reg176$effect[reg176$model=="OS"],low=reg176$ci_low[reg176$model=="OS"],high=reg176$ci_high[reg176$model=="OS"],p=reg176$p[reg176$model=="OS"]),
  data.frame(label="IMvigor210 OS, unadjusted",effect=cox210$effect[4],low=cox210$ci_low[4],high=cox210$ci_high[4],p=cox210$p[4]),
  data.frame(label="IMvigor210 OS, TMB-adjusted",effect=cox210$effect[6],low=cox210$ci_low[6],high=cox210$ci_high[6],p=cox210$p[6]),
  data.frame(label="IMvigor210 OS, TMB + cDC1",effect=cox210$effect[nrow(cox210)],low=cox210$ci_low[nrow(cox210)],high=cox210$ci_high[nrow(cox210)],p=cox210$p[nrow(cox210)])
)
write.csv(ddsurv,file.path(figdir,"Figure5_survival_forest_source_data.csv"),row.names=FALSE)

draw_fig5<-function() {
  layout(matrix(1:4,2,2,byrow=TRUE));par(mar=c(4.4,4.5,2.8,1),family="sans")
  km_panel(as.numeric(d176$pfs),d176$pfs_event_num,d176$cross_nonoverlap,"GSE176307 progression-free survival","Days");panel_label("A")
  km_panel(as.numeric(d176$overall_survival),d176$os_event_num,d176$cross_nonoverlap,"GSE176307 overall survival","Days");panel_label("B")
  km_panel(d210$os,d210$censOS,d210$cross_presentation_nonoverlap,"IMvigor210 overall survival","Months",max_time=25);panel_label("C")
  par(mar=c(4.4,4.5,2.8,1.0))
  cols<-c(COL["blue"],COL["blue"],COL["orange"],COL["orange"],COL["gray"])
  short_labels<-c("GSE PFS","GSE OS","IMvigor OS","IMvigor + TMB","IMvigor + TMB + cDC1")
  plot_forest_inside(short_labels,ddsurv$effect,ddsurv$low,ddsurv$high,ddsurv$p,xlim=c(.25,1.25),xlab="Hazard ratio per 1-SD",cols=cols,cex=.72)
  title("Continuous-score Cox models",cex.main=1.0);panel_label("D")
}
export_figure("Figure5_survival_associations",11.5,8.7,draw_fig5)

# Supplementary Figure S1: patient-level single-cell heatmap ----------------
draw_s1<-function() {
  par(mar=c(6.3,7.1,2.5,1.2),family="sans")
  pats<-sort(unique(scscores$patient));grps<-c("cDC1","cDC2","pDC","mature_migratory_DC","nonDC_myeloid","epithelial")
  mat<-matrix(NA,nrow=length(grps),ncol=length(pats),dimnames=list(group_lab[grps],pats))
  for(i in seq_len(nrow(scscores))) mat[group_lab[scscores$comparison_group[i]],as.character(scscores$patient[i])]<-scscores$raw_cross_presentation_nonoverlap_mean[i]
  z<-t(scale(t(mat)));z[!is.finite(z)]<-NA
  image(seq_len(ncol(z)),seq_len(nrow(z)),t(z[nrow(z):1,,drop=FALSE]),zlim=c(-2,2),
        col=colorRampPalette(c("#315C8A","#D8E6F0","white","#F4D5C3","#B95345"))(101),axes=FALSE,xlab="Patient",ylab="",main="Patient-level cross-presentation score (row z-score)")
  axis(1,seq_len(ncol(z)),colnames(z),las=2,cex.axis=.72);axis(2,seq_len(nrow(z)),rev(rownames(z)),las=1,tick=FALSE,cex.axis=.82)
  for(i in seq_len(nrow(z))) for(j in seq_len(ncol(z))) if(is.na(z[i,j])) text(j,nrow(z)-i+1,"x",col="#B8B8B8",cex=.65)
  legend("topright",legend=c("Lower","Higher","Missing"),fill=c("#315C8A","#B95345",NA),border=NA,pch=c(NA,NA,4),bty="n",cex=.74,xpd=NA)
}
export_figure("FigureS1_single_cell_patient_heatmap",11.5,5.0,draw_s1)

# Supplementary Figure S2: score distributions in IMvigor210 ---------------
draw_s2<-function() {
  layout(matrix(1:4,2,2,byrow=TRUE));par(mar=c(4.0,4.4,2.7,1),family="sans")
  for(i in seq_along(score_cols)) {
    score_boxplot(d210,score_cols[i],title=score_names[i],ylab="Standardized program score",annotate_auc=sm$auc[i])
    panel_label(LETTERS[i])
  }
}
export_figure("FigureS2_IMvigor_program_distributions",11.5,8.2,draw_s2)

# Manifest ------------------------------------------------------------------
manifest<-data.frame(
  figure=c("Figure 1","Figure 2","Figure 3","Figure 4","Figure 5","Figure S1","Figure S2"),
  stem=c("Figure1_study_design_and_cohorts","Figure2_single_cell_localization","Figure3_cross_cohort_response_validation",
         "Figure4_IMvigor_specificity_and_robustness","Figure5_survival_associations","FigureS1_single_cell_patient_heatmap","FigureS2_IMvigor_program_distributions"),
  status=c("Main","Main","Main","Main","Main","Supplementary","Supplementary"),stringsAsFactors=FALSE)
write.csv(manifest,file.path(figdir,"figure_manifest.csv"),row.names=FALSE)
cat("Figures written to",normalizePath(figdir,winslash="/",mustWork=FALSE),"\n")
