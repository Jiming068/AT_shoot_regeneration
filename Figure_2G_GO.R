# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({library(clusterProfiler);library(org.At.tair.db);library(AnnotationDbi);library(ggplot2);library(dplyr);library(tidyr)})

assignment_file <- file.path(project_dir, "02_results/132_WT_originalUMAP_C02418_four_state_gene_modules/04_dynamic_gene_four_module_assignments.csv")
cds_file <- file.path(project_dir, "02_results/125_WT_res04_C02418_monocle3_trajectory/22_C02418_monocle3_cds_original_UMAP_with_pseudotime.rds")
outdir <- Sys.getenv("FOUR_GO_OUTDIR",file.path(output_root, "133_WT_originalUMAP_four_gene_modules_GO_BP"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
modules <- c("Early module (cluster 0)","Cluster-2-specific module","Cluster-4-specific module","Late module (clusters 1/8)")
short_labels <- c("0 early","2 specific","4 specific","1/8 late");names(short_labels)<-modules

assignment <- read.csv(assignment_file,stringsAsFactors=FALSE,check.names=FALSE);assignment$module<-factor(assignment$module,levels=modules)
cds <- readRDS(cds_file);valid <- keys(org.At.tair.db,keytype="TAIR");universe <- intersect(toupper(rownames(cds)),valid);writeLines(universe,file.path(outdir,"02_GO_universe_all_CDS_TAIR_genes.txt"))

all_res<-list();summaries<-list()
for(m in modules){
 genes<-intersect(unique(toupper(assignment$gene_id[assignment$module==m])),universe);tag<-gsub("[^A-Za-z0-9]+","_",m);writeLines(genes,file.path(outdir,paste0("03_",tag,"_TAIR_genes_for_GO.txt")))
 ego<-enrichGO(gene=genes,universe=universe,OrgDb=org.At.tair.db,keyType="TAIR",ont="BP",pAdjustMethod="BH",pvalueCutoff=.05,qvalueCutoff=1,minGSSize=10,maxGSSize=500,readable=FALSE)
 ed<-as.data.frame(ego);if(nrow(ed))ed<-ed[is.finite(ed$pvalue)&ed$pvalue<=.05&is.finite(ed$p.adjust)&ed$p.adjust<=.05,,drop=FALSE]
 if(nrow(ed)){ed<-ed[order(ed$p.adjust,-ed$Count),];ed$module<-m;ed$input_genes<-length(genes);all_res[[m]]<-ed;write.csv(ed,file.path(outdir,paste0("04_",tag,"_significant_GO_BP.csv")),row.names=FALSE)}
 summaries[[m]]<-data.frame(module=m,input_genes=length(genes),significant_GO_BP_terms=nrow(ed),minimum_FDR=if(nrow(ed))min(ed$p.adjust)else NA_real_,stringsAsFactors=FALSE)
}
res<-bind_rows(all_res);if(!nrow(res))stop("No GO BP terms passed FDR <= 0.05");res$module<-factor(res$module,levels=modules);res<-res[order(res$module,res$p.adjust,-res$Count),];write.csv(res,file.path(outdir,"05_all_four_modules_significant_GO_BP.csv"),row.names=FALSE)
summary_df<-bind_rows(summaries);write.csv(summary_df,file.path(outdir,"06_GO_BP_enrichment_summary.csv"),row.names=FALSE)

# Top 12 significant terms per module.
top12<-res%>%group_by(module)%>%slice_min(p.adjust,n=12,with_ties=FALSE)%>%ungroup();top12$GeneRatio_numeric<-vapply(strsplit(top12$GeneRatio,"/"),function(x)as.numeric(x[1])/as.numeric(x[2]),numeric(1));top12$term_label<-paste0(top12$Description," [",top12$ID,"]");top12<-top12%>%arrange(module,desc(p.adjust));top12$term_label<-factor(top12$term_label,levels=unique(top12$term_label));top12$module_short<-factor(short_labels[as.character(top12$module)],levels=short_labels[modules])
p1<-ggplot(top12,aes(GeneRatio_numeric,term_label,size=Count,color=-log10(p.adjust)))+geom_point(alpha=.9)+facet_grid(module_short~.,scales="free_y",space="free_y",switch="y")+scale_color_gradient(low="#F6C6A8",high="#B4452D",name=expression(-log[10](FDR)))+scale_size(range=c(1.8,4.8),name="Gene count")+labs(x="Gene ratio",y=NULL)+theme_classic(base_family="Arial",base_size=9)+theme(axis.text.y=element_text(size=7),axis.title.x=element_text(face="bold"),strip.text.y=element_text(face="bold",size=8.2),strip.background=element_rect(fill="#F2F2F2",colour=NA),panel.spacing.y=grid::unit(.8,"mm"),legend.position="right",plot.margin=margin(4,5,4,4))

# Top five terms from every module in a unified significance heatmap.
selected<-res%>%group_by(module)%>%slice_min(p.adjust,n=5,with_ties=FALSE)%>%ungroup();term_order<-unique(selected$ID);mat<-expand_grid(ID=term_order,module=factor(modules,levels=modules))%>%left_join(res%>%select(ID,Description,module,p.adjust),by=c("ID","module"));desc_map<-setNames(res$Description,res$ID);mat$term<-paste0(desc_map[mat$ID]," [",mat$ID,"]");mat$score<-ifelse(is.na(mat$p.adjust),0,pmin(-log10(mat$p.adjust),10));mat$term<-factor(mat$term,levels=rev(unique(paste0(desc_map[term_order]," [",term_order,"]"))));mat$module_short<-factor(short_labels[as.character(mat$module)],levels=short_labels[modules])
p2<-ggplot(mat,aes(module_short,term,fill=score))+geom_tile(color="white",linewidth=.5)+scale_fill_gradient(low="#F1F1F1",high="#D96B4B",limits=c(0,10),name=expression(-log[10](FDR)))+labs(x=NULL,y=NULL)+theme_classic(base_family="Arial",base_size=9)+theme(axis.text.x=element_text(angle=28,hjust=1,face="bold",size=8),axis.text.y=element_text(size=7),axis.ticks=element_blank(),panel.border=element_rect(fill=NA,colour="#555555",linewidth=.4),legend.position="right",plot.margin=margin(4,5,4,4))

save2<-function(p,n,w,h){ggsave(file.path(outdir,paste0(n,".pdf")),p,width=w,height=h,device=cairo_pdf);ggsave(file.path(outdir,paste0(n,".png")),p,width=w,height=h,dpi=500,bg="white")}
save2(p1,"07_top12_GO_BP_dotplot",7.3,10.0);save2(p2,"08_top5_GO_BP_comparison_heatmap",6.3,7.0)

# Plain-text top-term digest for rapid interpretation.
digest<-unlist(lapply(modules,function(m){x<-res[res$module==m,];c(paste0("## ",m),paste0(head(x$Description,15)," [",head(x$ID,15),"] | FDR=",format(head(x$p.adjust,15),digits=3,scientific=TRUE)),"")}));writeLines(digest,file.path(outdir,"09_top15_GO_terms_digest.txt"))
writeLines(c("GO Biological Process enrichment for four original-UMAP pseudotime gene modules.",paste0("Universe: ",length(universe)," valid TAIR IDs among all genes in the trajectory CDS."),"Method: clusterProfiler::enrichGO; keyType=TAIR; ont=BP; pAdjustMethod=BH; pvalueCutoff=0.05; qvalueCutoff=1; minGSSize=10; maxGSSize=500.","Reported terms satisfy both raw P <= 0.05 and BH-adjusted P <= 0.05.","Figures were exported as PDF and PNG only."),file.path(outdir,"10_method_information.txt"));writeLines(capture.output(sessionInfo()),file.path(outdir,"11_sessionInfo.txt"))
