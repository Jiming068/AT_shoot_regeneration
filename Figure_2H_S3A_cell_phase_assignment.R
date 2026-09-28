# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(dplyr)})
input_rds <- file.path(project_dir, "01_obj/01_WT_snRNA/03_snRNA_C0-2-3-subclustered_res0.4_SIM4-SIM16.rds")
gene_csv <- file.path(project_dir, "02_results/80_WT_C02781_cell_cycle_Kotaro_Torii/00_Kotaro_Torii_cell_cycle_gene_sets_long.csv")
outdir <- Sys.getenv("KT_C02418_OUTDIR",file.path(output_root, "137_WT_C02418_cell_cycle_Kotaro_Torii"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
cluster_order<-c("0","2","4","1","8");phase_order<-c("G1","S","G2/M")
phase_cols<-c("G1"="#E9B949","S"="#DE6E56","G2/M"="#6489B9");confidence_cols<-c("Higher-confidence"="#507C8C","Boundary/low-margin"="#D8D8D8")
save_plot<-function(p,n,w,h){ggsave(file.path(outdir,paste0(n,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white");ggsave(file.path(outdir,paste0(n,".png")),p,width=w,height=h,dpi=500,bg="white")}
theme_pub<-function(){theme_classic(base_size=10.5,base_family="Arial")+theme(axis.text=element_text(color="#222222",size=9.3),axis.text.x=element_text(face="bold"),axis.title=element_text(face="bold",size=10.3),legend.title=element_blank(),legend.text=element_text(size=9),legend.key.height=grid::unit(.35,"cm"),legend.key.width=grid::unit(.44,"cm"),plot.margin=margin(5,6,4,5))}

long<-read.csv(gene_csv,stringsAsFactors=FALSE,check.names=FALSE);long$Gene<-toupper(trimws(long$Gene));long$Phase<-trimws(long$Phase)
long<-long[long$Phase%in%phase_order & grepl("^AT[1-5CM]G[0-9]{5}$",long$Gene),,drop=FALSE]
gene_sets<-setNames(lapply(phase_order,function(ph)unique(long$Gene[long$Phase==ph])),phase_order)
dup_within<-long[duplicated(long[,c("Gene","Phase")]),,drop=FALSE];cross<-aggregate(Phase~Gene,long,function(x)paste(unique(x),collapse=";"));cross$n_phases<-lengths(strsplit(cross$Phase,";",fixed=TRUE));cross<-cross[cross$n_phases>1,,drop=FALSE]
write.csv(long,file.path(outdir,"00_Kotaro_Torii_cell_cycle_gene_sets_long.csv"),row.names=FALSE);write.csv(cross,file.path(outdir,"01_cross_phase_gene_conflicts.csv"),row.names=FALSE)

obj<-readRDS(input_rds);stopifnot("subcluster"%in%colnames(obj@meta.data));obj$analysis_cluster<-as.character(obj$subcluster);obj<-subset(obj,subset=analysis_cluster%in%cluster_order);obj$analysis_cluster<-factor(obj$analysis_cluster,levels=cluster_order)
DefaultAssay(obj)<-"RNA";obj<-NormalizeData(obj,verbose=FALSE)
present<-lapply(gene_sets,intersect,y=rownames(obj));qc<-data.frame(phase=phase_order,supplied=lengths(gene_sets),detected=lengths(present),detection_rate=lengths(present)/lengths(gene_sets),detected_genes=vapply(present,paste,collapse=";",FUN.VALUE=character(1)))
write.csv(qc,file.path(outdir,"02_gene_set_detection_QC.csv"),row.names=FALSE);stopifnot(all(lengths(present)>=10),nrow(cross)==0,nrow(dup_within)==0)
set.seed(20260814);safe<-c("G1"="G1","S"="S","G2/M"="G2M")
for(ph in phase_order){nm<-paste0("KT_CC_",safe[[ph]]);obj<-AddModuleScore(obj,features=list(present[[ph]]),assay="RNA",name=nm,seed=20260814);colnames(obj@meta.data)[colnames(obj@meta.data)==paste0(nm,"1")]<-nm}
score_cols<-paste0("KT_CC_",unname(safe[phase_order]));raw<-as.matrix(obj@meta.data[,score_cols,drop=FALSE]);zs<-scale(raw);colnames(zs)<-phase_order
obj@meta.data[,paste0("KT_z_",unname(safe[phase_order]))]<-zs;ord<-t(apply(zs,1,function(x)order(x,decreasing=TRUE)));obj$KT_dominant_phase<-factor(phase_order[ord[,1]],levels=phase_order)
obj$KT_phase_margin<-zs[cbind(seq_len(nrow(zs)),ord[,1])]-zs[cbind(seq_len(nrow(zs)),ord[,2])];obj$KT_phase_confidence<-factor(ifelse(obj$KT_phase_margin>=.25,"Higher-confidence","Boundary/low-margin"),levels=c("Higher-confidence","Boundary/low-margin"));obj$KT_cycling_score<-rowMeans(zs[,c("S","G2/M"),drop=FALSE])-zs[,"G1"]
md<-obj@meta.data;md$cell<-rownames(md);cell_cols<-c("cell","analysis_cluster",intersect(c("date","orig.ident"),colnames(md)),"KT_dominant_phase","KT_phase_margin","KT_phase_confidence","KT_cycling_score",score_cols,paste0("KT_z_",unname(safe[phase_order])))
write.csv(md[,cell_cols],file.path(outdir,"03_cell_level_Kotaro_Torii_cycle_scores.csv"),row.names=FALSE)

pc<-as.data.frame(table(subcluster=md$analysis_cluster,phase=md$KT_dominant_phase));pc$subcluster<-factor(pc$subcluster,levels=cluster_order);pc$phase<-factor(pc$phase,levels=phase_order);pc$proportion<-ave(pc$Freq,pc$subcluster,FUN=function(x)x/sum(x));pc$percentage<-100*pc$proportion;pc<-pc[order(pc$subcluster,pc$phase),];write.csv(pc,file.path(outdir,"04_phase_counts_and_proportions.csv"),row.names=FALSE)
cc<-as.data.frame(table(subcluster=md$analysis_cluster,confidence=md$KT_phase_confidence));cc$subcluster<-factor(cc$subcluster,levels=cluster_order);cc$proportion<-ave(cc$Freq,cc$subcluster,FUN=function(x)x/sum(x));write.csv(cc,file.path(outdir,"05_phase_assignment_confidence.csv"),row.names=FALSE)
ms<-aggregate(md[,score_cols],list(subcluster=md$analysis_cluster),mean);ms$subcluster<-factor(ms$subcluster,levels=cluster_order);ms<-ms[order(ms$subcluster),];write.csv(ms,file.path(outdir,"06_mean_module_scores_by_subcluster.csv"),row.names=FALSE)
summ<-do.call(rbind,lapply(cluster_order,function(cl){d<-md[as.character(md$analysis_cluster)==cl,,drop=FALSE];pp<-prop.table(table(factor(d$KT_dominant_phase,levels=phase_order)));data.frame(subcluster=cl,n_cells=nrow(d),G1_pct=100*pp["G1"],S_pct=100*pp["S"],G2M_pct=100*pp["G2/M"],cycling_pct=100*sum(pp[c("S","G2/M")]),mean_cycling_score=mean(d$KT_cycling_score),higher_confidence_pct=100*mean(d$KT_phase_confidence=="Higher-confidence"),median_phase_margin=median(d$KT_phase_margin),dominant_phase=names(which.max(pp)))}));rownames(summ)<-NULL;write.csv(summ,file.path(outdir,"07_subcluster_cycle_summary.csv"),row.names=FALSE)
if("date"%in%colnames(md)){dc<-as.data.frame(table(subcluster=md$analysis_cluster,date=md$date));dc$proportion<-ave(dc$Freq,dc$subcluster,FUN=function(x)x/sum(x));write.csv(dc,file.path(outdir,"08_date_composition_by_subcluster.csv"),row.names=FALSE);pd<-as.data.frame(table(subcluster=md$analysis_cluster,date=md$date,phase=md$KT_dominant_phase));pd$proportion<-ave(pd$Freq,interaction(pd$subcluster,pd$date),FUN=function(x)if(sum(x)==0)0 else x/sum(x));write.csv(pd,file.path(outdir,"09_phase_by_subcluster_and_date.csv"),row.names=FALSE)}

p<-ggplot(pc,aes(subcluster,proportion,fill=phase))+geom_col(width=.78,color="white",linewidth=.3)+scale_fill_manual(values=phase_cols,breaks=phase_order,drop=FALSE)+scale_y_continuous(labels=function(x)paste0(round(100*x),"%"),breaks=seq(0,1,.25),expand=c(0,0))+labs(x="Subcluster",y="Fraction of nuclei")+theme_pub()+theme(legend.position="top",legend.justification="left")
save_plot(p,"10_Kotaro_Torii_phase_composition",4.5,3.25);save_plot(p+geom_text(aes(label=ifelse(percentage>=4,sprintf("%.1f%%",percentage),"")),position=position_stack(vjust=.5),color="white",size=3.15,family="Arial",fontface="bold"),"11_Kotaro_Torii_phase_composition_with_percentages",4.5,3.25)
pcf<-ggplot(cc,aes(subcluster,proportion,fill=confidence))+geom_col(width=.78,color="white",linewidth=.3)+scale_fill_manual(values=confidence_cols,drop=FALSE)+scale_y_continuous(labels=function(x)paste0(round(100*x),"%"),breaks=seq(0,1,.25),expand=c(0,0))+labs(x="Subcluster",y="Fraction of nuclei")+theme_pub()+theme(legend.position="top",legend.justification="left");save_plot(pcf,"12_phase_assignment_confidence",4.7,3.25)
if("umap"%in%names(obj@reductions)){pu<-DimPlot(obj,reduction="umap",group.by="KT_dominant_phase",cols=phase_cols,pt.size=.25)+ggtitle("Dominant cell-cycle module")+theme_void(base_size=10.5,base_family="Arial")+theme(legend.position="right",legend.title=element_blank(),legend.text=element_text(size=9),plot.title=element_text(size=10.3,face="bold"));save_plot(pu,"13_Kotaro_Torii_phase_UMAP",4.7,3.8)}
saveRDS(obj,file.path(outdir,"14_C02418_Kotaro_Torii_cycle_scores.rds"),compress=TRUE)
writeLines(c("Kotaro/Torii G1, S and G2/M module scoring for subclusters 0 -> 2 -> 4 -> 1 -> 8.",paste0("Input: ",input_rds),"RNA was log-normalized and AddModuleScore was run with seed 20260814.","Raw module scores were standardized across the selected cells; the dominant phase is the highest Z-score.","A top-versus-second score margin >= 0.25 denotes a higher-confidence assignment.","Dominant module indicates relative transcriptomic resemblance and does not directly measure cell-cycle duration or division rate.","Output figures are PDF and PNG only."),file.path(outdir,"15_analysis_notes_and_guardrails.txt"));writeLines(capture.output(sessionInfo()),file.path(outdir,"16_sessionInfo.txt"));print(qc[,1:4]);print(summ)
