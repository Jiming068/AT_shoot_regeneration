# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({library(monocle3);library(Matrix);library(ggplot2);library(dplyr);library(tidyr);library(org.At.tair.db);library(AnnotationDbi)})

cds_file <- file.path(project_dir, "02_results/125_WT_res04_C02418_monocle3_trajectory/22_C02418_monocle3_cds_original_UMAP_with_pseudotime.rds")
dynamic_file <- file.path(project_dir, "02_results/130_WT_originalUMAP_C02418_pseudotime_gene_patterns/04_significant_dynamic_genes.csv")
outdir <- Sys.getenv("FOUR_MODULE_OUTDIR",file.path(output_root, "132_WT_originalUMAP_C02418_four_state_gene_modules"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE);set.seed(132)
states<-c("0","2","4","1","8");module_names<-c("Early module (cluster 0)","Cluster-2-specific module","Cluster-4-specific module","Late module (clusters 1/8)");n_bins<-60L

cds<-readRDS(cds_file);pt<-pseudotime(cds);keep<-is.finite(pt);cds<-cds[,keep];pt<-pt[keep];counts<-assays(cds)[["counts"]]
sig<-read.csv(dynamic_file,stringsAsFactors=FALSE);genes<-intersect(sig$gene_id,rownames(counts))
ord<-order(pt,names(pt));sorted_bins<-ceiling(seq_along(ord)*n_bins/length(ord));bin_id<-integer(length(ord));bin_id[ord]<-sorted_bins;names(bin_id)<-names(pt)
sf<-as.numeric(colData(cds)$Size_Factor);names(sf)<-colnames(cds);if(any(!is.finite(sf))||any(sf<=0)){lib<-Matrix::colSums(counts);sf<-lib/median(lib)}
norm<-counts[genes,,drop=FALSE]%*%Diagonal(x=1/sf[colnames(cds)]);rownames(norm)<-genes;colnames(norm)<-colnames(cds);norm@x<-log1p(norm@x)
be<-sapply(seq_len(n_bins),function(b)Matrix::rowMeans(norm[,names(bin_id)[bin_id==b],drop=FALSE]));rownames(be)<-genes;colnames(be)<-paste0("bin",seq_len(n_bins));bin_pt<-sapply(seq_len(n_bins),function(b)median(pt[names(bin_id)[bin_id==b]]));z<-t(scale(t(be)));z[!is.finite(z)]<-0;z<-pmax(pmin(z,2.5),-2.5)

# Unsupervised four-curve clustering, named only after ordering center peaks in pseudotime.
km<-kmeans(z,centers=4,nstart=100,iter.max=300);centers<-km$centers;center_peak_bin<-apply(centers,1,which.max);cluster_order<-as.integer(names(sort(center_peak_bin)));cluster_to_module<-setNames(module_names,cluster_order);module<-unname(cluster_to_module[as.character(km$cluster)])
ordered_centers<-centers[cluster_order,,drop=FALSE]

state_pt_raw<-sapply(states,function(s)median(pt[as.character(colData(cds)$subcluster)==s]));state_ref<-data.frame(state=c("0","2","4","1/8"),median_pseudotime=c(state_pt_raw["0"],state_pt_raw["2"],state_pt_raw["4"],mean(state_pt_raw[c("1","8")])));state_ref$bin<-vapply(state_ref$median_pseudotime,function(v)which.min(abs(bin_pt-v)),integer(1));write.csv(state_ref,file.path(outdir,"02_state_pseudotime_reference.csv"),row.names=FALSE)
center_summary<-data.frame(module=module_names,kmeans_cluster=cluster_order,peak_bin=apply(ordered_centers,1,which.max),peak_pseudotime=bin_pt[apply(ordered_centers,1,which.max)])
center_summary$nearest_state<-state_ref$state[vapply(center_summary$peak_pseudotime,function(v)which.min(abs(state_ref$median_pseudotime-v)),integer(1))]
write.csv(center_summary,file.path(outdir,"03_four_module_center_peak_validation.csv"),row.names=FALSE)

assignment<-data.frame(gene_id=rownames(z),module=factor(module,levels=module_names),kmeans_cluster=km$cluster,peak_bin=apply(z,1,which.max),peak_pseudotime=bin_pt[apply(z,1,which.max)],trajectory_autocorrelation=sig$trajectory_lag1_autocorrelation[match(rownames(z),sig$gene_id)],q_value=sig$q_value[match(rownames(z),sig$gene_id)])
assignment<-assignment[order(assignment$module,assignment$peak_bin,-assignment$trajectory_autocorrelation),];write.csv(assignment,file.path(outdir,"04_dynamic_gene_four_module_assignments.csv"),row.names=FALSE)
for(m in module_names)writeLines(assignment$gene_id[assignment$module==m],file.path(outdir,paste0("05_",gsub("[^A-Za-z0-9]+","_",m),"_genes.txt")))
module_counts<-as.data.frame(table(assignment$module));colnames(module_counts)<-c("module","n_genes");write.csv(module_counts,file.path(outdir,"06_gene_counts_by_module.csv"),row.names=FALSE)
write.csv(data.frame(pseudotime=bin_pt,t(ordered_centers),check.names=FALSE),file.path(outdir,"07_four_module_kmeans_centers.csv"),row.names=FALSE)

cols<-c("Early module (cluster 0)"="#3E5C91","Cluster-2-specific module"="#E5A04B","Cluster-4-specific module"="#D84A5B","Late module (clusters 1/8)"="#8E55A3")
clong<-bind_rows(lapply(seq_along(module_names),function(i)data.frame(module=factor(module_names[i],levels=module_names),pseudotime=bin_pt,mean_Zscore=ordered_centers[i,])))
pcurve<-ggplot(clong,aes(pseudotime,mean_Zscore,color=module))+geom_hline(yintercept=0,lty=2,color="#BDBDBD",linewidth=.35)+geom_vline(data=state_ref,aes(xintercept=median_pseudotime),color="#D0D0D0",linewidth=.4)+geom_line(linewidth=1.05)+scale_color_manual(values=cols)+scale_x_continuous(breaks=state_ref$median_pseudotime,labels=state_ref$state)+labs(x="Trajectory state",y="Mean gene-wise scaled expression",color=NULL)+theme_classic(base_family="Arial",base_size=10)+theme(axis.title=element_text(face="bold"),axis.text.x=element_text(face="bold"),legend.position="right",legend.text=element_text(size=8.2),plot.margin=margin(4,4,4,4))
save2<-function(p,n,w,h){ggsave(file.path(outdir,paste0(n,".pdf")),p,width=w,height=h,device=cairo_pdf);ggsave(file.path(outdir,paste0(n,".png")),p,width=w,height=h,dpi=500,bg="white")}
save2(pcurve,"08_four_gene_modules_over_pseudotime",6.5,3.8)

# All dynamic genes.
zl<-as.data.frame(as.table(z[assignment$gene_id,,drop=FALSE]),stringsAsFactors=FALSE);colnames(zl)<-c("gene_id","bin","Zscore");zl$module<-assignment$module[match(zl$gene_id,assignment$gene_id)];zl$gene<-factor(zl$gene_id,levels=rev(assignment$gene_id));zl$bin_num<-as.integer(sub("bin","",zl$bin))
pall<-ggplot(zl,aes(bin_num,gene,fill=Zscore))+geom_raster()+geom_vline(xintercept=state_ref$bin+.5,color="white",linewidth=.3)+facet_grid(module~.,scales="free_y",space="free_y",switch="y")+scale_fill_gradient2(low="#7896C8",mid="#F7F7F7",high="#E58A82",midpoint=0,limits=c(-2.5,2.5),name="Gene-wise\nZ-score")+scale_x_continuous(breaks=state_ref$bin,labels=state_ref$state,expand=c(0,0))+labs(x="Trajectory state",y=NULL)+theme_classic(base_family="Arial",base_size=9)+theme(axis.text.x=element_text(face="bold"),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(face="bold"),strip.text.y=element_text(face="bold",size=8),strip.background=element_rect(fill="#F2F2F2",colour=NA),panel.spacing.y=grid::unit(.7,"mm"),plot.margin=margin(4,5,4,4))
save2(pall,"09_all_dynamic_genes_four_module_heatmap",5.9,6.4)

# Top 20 genes per module with symbols.
top<-assignment%>%group_by(module)%>%slice_max(trajectory_autocorrelation,n=20,with_ties=FALSE)%>%ungroup();top<-top[order(top$module,top$peak_bin,-top$trajectory_autocorrelation),]
sm<-AnnotationDbi::select(org.At.tair.db,keys=unique(top$gene_id),keytype="TAIR",columns="SYMBOL");sm<-sm[!duplicated(sm$TAIR),];top$symbol<-sm$SYMBOL[match(top$gene_id,sm$TAIR)];top$display_gene<-ifelse(is.na(top$symbol)|top$symbol=="",top$gene_id,top$symbol);top$display_gene<-make.unique(top$display_gene,sep="_");write.csv(top,file.path(outdir,"10_top20_representative_genes_per_module.csv"),row.names=FALSE)
tl<-as.data.frame(as.table(z[top$gene_id,,drop=FALSE]),stringsAsFactors=FALSE);colnames(tl)<-c("gene_id","bin","Zscore");tl$module<-top$module[match(tl$gene_id,top$gene_id)];tl$gene<-factor(top$display_gene[match(tl$gene_id,top$gene_id)],levels=rev(top$display_gene));tl$bin_num<-as.integer(sub("bin","",tl$bin))
ptop<-ggplot(tl,aes(bin_num,gene,fill=Zscore))+geom_tile()+geom_vline(xintercept=state_ref$bin+.5,color="white",linewidth=.3)+facet_grid(module~.,scales="free_y",space="free_y",switch="y")+scale_fill_gradient2(low="#7896C8",mid="#F7F7F7",high="#E58A82",midpoint=0,limits=c(-2.5,2.5),name="Gene-wise\nZ-score")+scale_x_continuous(breaks=state_ref$bin,labels=state_ref$state,expand=c(0,0))+labs(x="Trajectory state",y=NULL)+theme_classic(base_family="Arial",base_size=9)+theme(axis.text.x=element_text(face="bold"),axis.text.y=element_text(face="italic",size=6.8),axis.ticks=element_blank(),axis.title.x=element_text(face="bold"),strip.text.y=element_text(face="bold",size=8),strip.background=element_rect(fill="#F2F2F2",colour=NA),panel.spacing.y=grid::unit(.7,"mm"),plot.margin=margin(4,5,4,4))
save2(ptop,"11_top20_genes_four_module_heatmap",5.9,9.5)

writeLines(c("Four temporal gene modules along the original-UMAP Monocle 3 trajectory 0 -> 2 -> 4 -> 1/8.","The significant dynamic-gene set and 60 equal-frequency pseudotime bins were inherited from analysis 130.","Gene-wise Z-scored curves were clustered by k-means (k=4, nstart=100) and modules were named according to the temporal order of center-curve peaks.","The center-peak validation table must be consulted before biological interpretation; module labels are not imposed on individual genes by cell identity.","Figures were exported as PDF and PNG only."),file.path(outdir,"12_method_information.txt"));writeLines(capture.output(sessionInfo()),file.path(outdir,"13_sessionInfo.txt"))
