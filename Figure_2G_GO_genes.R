# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({library(monocle3);library(Matrix);library(org.At.tair.db);library(AnnotationDbi);library(ggplot2);library(dplyr);library(tidyr)})

go_file <- file.path(project_dir, "02_results/133_WT_originalUMAP_four_gene_modules_GO_BP/05_all_four_modules_significant_GO_BP.csv")
assignment_file <- file.path(project_dir, "02_results/132_WT_originalUMAP_C02418_four_state_gene_modules/04_dynamic_gene_four_module_assignments.csv")
cds_file <- file.path(project_dir, "02_results/125_WT_res04_C02418_monocle3_trajectory/22_C02418_monocle3_cds_original_UMAP_with_pseudotime.rds")
outdir <- Sys.getenv("GO3GENE_OUTDIR",file.path(output_root, "136_WT_four_modules_15GO_three_genes_expression_heatmap"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

module_names <- c("Early module (cluster 0)","Cluster-2-specific module","Cluster-4-specific module","Late module (clusters 1/8)")
sel <- data.frame(
 ID=c("GO:0000280","GO:0007059","GO:0000911","GO:0000281","GO:0010073","GO:0010629","GO:0006643","GO:0006644","GO:0046486","GO:0002181","GO:0051604","GO:0009117","GO:0015979","GO:0019684","GO:0009768"),
 display=c("nuclear division","chromosome segregation","cytokinesis by cell plate formation","mitotic cytokinesis","meristem maintenance","negative regulation of gene expression","membrane lipid metabolic process","phospholipid metabolic process","glycerolipid metabolic process","cytoplasmic translation","protein maturation","nucleotide metabolism","photosynthesis","photosynthesis, light reaction","photosynthesis, light harvesting in photosystem I"),
 focal_module=c(rep(module_names[1],6),rep(module_names[2],3),rep(module_names[3],3),rep(module_names[4],3)),
 category=c(rep("Early | Cell division",4),rep("Early | Meristem/regulation",2),rep("Cluster 2 | Lipid metabolism",3),rep("Cluster 4 | Biosynthesis",3),rep("Late | Photosynthesis",3)),stringsAsFactors=FALSE)

go <- read.csv(go_file,stringsAsFactors=FALSE,check.names=FALSE)
sg <- left_join(sel,go %>% select(module,ID,Description,pvalue,p.adjust,geneID,Count),by=c("focal_module"="module","ID"="ID"))
if(any(is.na(sg$geneID))) stop("One or more selected GO terms are absent from their focal module")
assign <- read.csv(assignment_file,stringsAsFactors=FALSE,check.names=FALSE)

cds <- readRDS(cds_file); states <- c("0","2","4","1","8")
keep <- as.character(colData(cds)$subcluster) %in% states
cds <- cds[,keep]; counts <- assays(cds)[["counts"]]
sf <- as.numeric(colData(cds)$Size_Factor); names(sf) <- colnames(cds)
if(any(!is.finite(sf)) || any(sf<=0)){lib <- Matrix::colSums(counts);sf <- lib/median(lib)}
all_candidates <- unique(unlist(strsplit(sg$geneID,"/",fixed=TRUE)))
all_candidates <- intersect(all_candidates,rownames(counts))
norm <- counts[all_candidates,,drop=FALSE] %*% Diagonal(x=1/sf[colnames(cds)])
avg <- sapply(states,function(s) Matrix::rowMeans(norm[,as.character(colData(cds)$subcluster)==s,drop=FALSE]))
rownames(avg)<-all_candidates;colnames(avg)<-states
z <- t(scale(t(avg)));z[!is.finite(z)]<-0;z<-pmax(pmin(z,2.5),-2.5)

sm <- AnnotationDbi::select(org.At.tair.db,keys=all_candidates,keytype="TAIR",columns=c("SYMBOL","GENENAME"))
sm <- sm[!duplicated(sm$TAIR),];sym <- setNames(sm$SYMBOL,sm$TAIR);gname <- setNames(sm$GENENAME,sm$TAIR)
target_score <- function(g,m){
 if(m==module_names[1]) z[g,"0"] else if(m==module_names[2]) z[g,"2"] else if(m==module_names[3]) z[g,"4"] else mean(z[g,c("1","8")])
}
used <- character(); chosen <- list()
for(i in seq_len(nrow(sg))){
 cand <- intersect(unique(strsplit(sg$geneID[i],"/",fixed=TRUE)[[1]]),all_candidates)
 a <- assign[assign$gene_id %in% cand & assign$module==sg$focal_module[i],,drop=FALSE]
 if(nrow(a)<3) a <- assign[assign$gene_id %in% cand,,drop=FALSE]
 a <- a[!duplicated(a$gene_id),,drop=FALSE]
 a$state_specificity <- vapply(a$gene_id,target_score,numeric(1),m=sg$focal_module[i])
 a$has_symbol <- !is.na(sym[a$gene_id]) & nzchar(sym[a$gene_id])
 a$rank_score <- scale(a$state_specificity)[,1] + 0.35*scale(a$trajectory_autocorrelation)[,1] + ifelse(a$has_symbol,0.5,0)
 a$rank_score[!is.finite(a$rank_score)] <- a$state_specificity[!is.finite(a$rank_score)]
 a <- a[order(-a$rank_score,-a$state_specificity,-a$trajectory_autocorrelation),,drop=FALSE]
 fresh <- a$gene_id[!a$gene_id %in% used];pick <- head(unique(c(fresh,a$gene_id)),3)
 if(length(pick)<3) stop("Fewer than three candidates for ",sg$ID[i])
 used <- c(used,pick)
 chosen[[i]] <- data.frame(category=sg$category[i],GO_ID=sg$ID[i],GO_term=sg$display[i],focal_module=sg$focal_module[i],gene_id=pick,symbol=unname(sym[pick]),gene_name=unname(gname[pick]),state_specificity=vapply(pick,target_score,numeric(1),m=sg$focal_module[i]),trajectory_autocorrelation=assign$trajectory_autocorrelation[match(pick,assign$gene_id)],stringsAsFactors=FALSE)
}
genes <- bind_rows(chosen);genes$display_gene <- ifelse(is.na(genes$symbol)|genes$symbol=="",genes$gene_id,genes$symbol)
genes$display_gene <- make.unique(genes$display_gene,sep="_")
write.csv(sg %>% select(category,ID,display,focal_module,Description,pvalue,p.adjust,Count),file.path(outdir,"02_selected_15_GO_terms.csv"),row.names=FALSE)
write.csv(genes,file.path(outdir,"03_three_representative_genes_per_GO.csv"),row.names=FALSE)

zmat <- z[genes$gene_id,states,drop=FALSE]
row_labels <- paste0(genes$display_gene,"  |  ",genes$GO_term)
long <- as.data.frame(as.table(zmat),stringsAsFactors=FALSE);colnames(long)<-c("gene_id","state","Zscore")
long$state<-factor(long$state,levels=states);long$category<-factor(genes$category[match(long$gene_id,genes$gene_id)],levels=unique(sel$category));long$GO_term<-factor(genes$GO_term[match(long$gene_id,genes$gene_id)],levels=sel$display)
long$gene_label<-factor(row_labels[match(long$gene_id,genes$gene_id)],levels=rev(row_labels))
write.csv(cbind(genes,as.data.frame(zmat,check.names=FALSE)),file.path(outdir,"04_representative_gene_mean_expression_Zscore.csv"),row.names=FALSE)

p <- ggplot(long,aes(state,gene_label,fill=Zscore))+geom_tile(color="white",linewidth=.45)+facet_grid(GO_term~.,scales="free_y",space="free_y")+
 scale_fill_gradient2(low="#8FAFD3",mid="#F5F5F5",high="#DE7D73",midpoint=0,limits=c(-2.5,2.5),name="Gene-wise\nZ-score")+
 labs(x=NULL,y=NULL)+theme_classic(base_family="Arial",base_size=9)+theme(axis.text.x=element_text(face="bold",size=9),axis.text.y=element_text(face="italic",size=6.7),axis.ticks=element_blank(),strip.text.y=element_blank(),strip.background=element_blank(),panel.spacing.y=grid::unit(.65,"mm"),legend.position="right",plot.margin=margin(4,5,4,4))
ggsave(file.path(outdir,"05_three_genes_per_GO_expression_heatmap.pdf"),p,width=8.4,height=9.3,device=cairo_pdf)
ggsave(file.path(outdir,"05_three_genes_per_GO_expression_heatmap.png"),p,width=8.4,height=9.3,dpi=500,bg="white")

# A cleaner gene-symbol-only version; GO membership remains available in output table 03.
long$gene_label <- factor(genes$display_gene[match(long$gene_id,genes$gene_id)],levels=rev(genes$display_gene))
p2 <- p %+% long
ggsave(file.path(outdir,"06_three_genes_per_GO_expression_heatmap_gene_names_only.pdf"),p2,width=5.2,height=9.3,device=cairo_pdf)
ggsave(file.path(outdir,"06_three_genes_per_GO_expression_heatmap_gene_names_only.png"),p2,width=5.2,height=9.3,dpi=500,bg="white")
writeLines(c("Three representative genes were selected for each of 15 manually selected GO terms.","Candidates were required to belong to the GO enrichment gene set and were preferentially retained from the corresponding trajectory module.","Selection prioritized named genes, expression specificity for the focal state, and trajectory autocorrelation, while minimizing reuse across GO terms.","Expression is size-factor normalized, log1p transformed, averaged across cells in clusters 0, 2, 4, 1, and 8, then gene-wise Z-scored and capped at +/-2.5.","Figures were exported as PDF and PNG only."),file.path(outdir,"07_method_information.txt"))
writeLines(capture.output(sessionInfo()),file.path(outdir,"08_sessionInfo.txt"))
