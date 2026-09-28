# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

###### Ö²ÎïTFºÍ°Ð»ùÒò¹²±í´ï·ÖÎö£¨Ê¾ÀýÎªÄâÄÏ½æ£©single-cell regulatory network calculate in plant (SERCP)
###### Àî½à
###### 2022-10-17

# 1. ´«Èë²ÎÊý ------------------------------------------------------------------
print("** 1.1 ²ÎÊý´«Èë **")
args = commandArgs(T)
dir <- args[1]
outdir <- args[2]
rds <- args[3]
type <- args[4]
f2 <- args[5]
regulon <- args[6]
threshold.p <- as.vector(args[7])
threshold.n <- as.numeric(args[8])
#threshold.n <- -(threshold.n)
#dir <-file.path(external_data_dir, "26.TF_Plant")
#outdir <- "out_1021"
#rds <- file.path(external_data_dir, "test.rds")
#type <- "seurat_clusters"
#f2 <- file.path(external_data_dir, "Top20-DEGs.txt")
#regulon <- file.path(external_data_dir, "03.TF_target_all.txt")
#threshold.p <- 0.9
#threshold.n <- -0.9

print("¹¤×÷Â·¾¶: ")
print(dir)
print("½á¹ûÂ·¾¶: ")
print(outdir)
print("rdsÂ·¾¶: ")
print(rds)
print("°´ÕÕ·ÖÀà¼ÆËã¾ùÖµ: ")
print(type)
print("²îÒì»ùÒòÁÐ±í: ")
print(f2)
print("µ÷¿Ø¹ØÏµÁÐ±í: ")
print(regulon)
print("Õýµ÷¿ØµÄãÐÖµ: ")
print(threshold.p)
print("¸ºµ÷¿ØµÄãÐÖµ: ")
print(threshold.n)

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(openxlsx)
library(ComplexHeatmap)
library(grid)
library(circlize)
library("GetoptLong")
library(ape)
library(factoextra)
library(cowplot)
library(reshape2)
library("GetoptLong")
library(data.table)

plan()
plan("multiprocess", workers = as.numeric(4))
plan()
options(future.globals.maxSize= 429496729600)

dir.create(outdir,recursive = T)

# 2. µ¼ÈëÊý¾Ý ------------------------------------------------------------------
print("** 2. µ¼ÈëÊý¾Ý **")
obj <- readRDS(rds)
DefaultAssay(obj) <- "RNA"
degs <- read.csv(f2,header = T)
regu <- read.table(regulon,header = T)[,1:2]

## Ö»±£ÁôÔÚ²îÒì»ùÒòÖÐµÄTF
regu.degs <- regu[regu$TF %in% degs$gene & regu$target %in% rownames(obj),]
regu.degs.uniq <- regu.degs[!duplicated(regu.degs),]
#3. ¼ÆËã¾ùÖµ -------------------------------------------------------------------
print("** 3. ¼ÆËã¾ùÖµ **")
Exp <- as.data.frame(AverageExpression(obj,group.by = type,assays = "RNA"))
Exp <- Exp[rowSums(Exp) > 0,]
umap <- as.data.frame(obj@reductions$umap@cell.embeddings)
umap$cell <- rownames(umap)
meta.raw <- as.data.frame(obj@meta.data)
meta.raw$cell <- rownames(meta.raw)
meta <- meta.raw[,c("cell",type)]
rownames(meta) <- meta$cell

# 4. ¼ÆËãTFºÍ°Ð»ùÒòµÄÏà¹ØÐÔ ----------------------------------------------------
print("** 4. ¼ÆËãTFºÍ°Ð»ùÒòµÄÏà¹ØÐÔ **")
regu.degs.uniq <- regu.degs.uniq[regu.degs.uniq$target %in% rownames(Exp),]
outRes <- as.data.frame(matrix(nrow = nrow(regu.degs.uniq),ncol = 4,data = 0))
outRes[,1:2] <- regu.degs.uniq
for (i in 1:nrow(regu.degs.uniq)) {
  tf <- outRes[i,1]
  gene <- outRes[i,2]
  r <- cor.test(t(Exp[tf,]),t(Exp[gene,]))$estimate
  p <- cor.test(t(Exp[tf,]),t(Exp[gene,]))$p.value
  outRes[i,3] <- r
  outRes[i,4] <- p
}

# 5. ±£´æÏà¹ØÐÔµÄ½á¹û ----------------------------------------------------------
print("** 5. ±£´æÏà¹ØÐÔµÄ½á¹û **")
write.table(outRes,paste0(outdir,"/1.All_cor.rawdatra.txt"),row.names = F,sep = "\t",quote = F)
outRes[is.na(outRes)] <- 1
pc <- outRes[outRes$p.value < 0.05 & outRes$cor > threshold.p,]
nc <- outRes[outRes$p.value < 0.05 & outRes$cor < threshold.n,]
write.table(pc,paste0(outdir,"/2.Positive_cor.rawdatra.txt"),row.names = F,sep = "\t",quote = F)
write.table(nc,paste0(outdir,"/3.Negative_cor.rawdatra.txt"),row.names = F,sep = "\t",quote = F)
print(paste0(" p < 0.05 & cor > ",threshold.p," : ",nrow(outRes[outRes$p.value < 0.05 & outRes$cor > threshold.p,])))
print(paste0(" p < 0.05 & cor < ",threshold.n," : ",nrow(outRes[outRes$p.value < 0.05 & outRes$cor < threshold.n,])))

# 6. Í³¼ÆÃ¿¸öclusterÖÐµÄÕýµ÷¿Ø -------------------------------------------------
print("** 6. Í³¼ÆÃ¿¸öclusterÖÐµÄÕýµ÷¿Ø **")
TF.mean <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster))+3,nrow = 1))
TF.Regulon.mean <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster))+3,nrow = 1))

po.res <- list()
for (i in unique(degs$cluster)) {
  ## Ã¿¸öclusterµ¥¶À´´½¨Ä¿Â¼
  print(paste0("+++++++++++++++++++++ ## Run ",i," ## +++++++++++++++++++++"))
  outdir1 <- paste0(outdir,"/Positive.regulon/Cluster",i)
  dir.create(outdir1,recursive = T)

  cluster.degs <- degs[degs$cluster == i ,]
  cluster.regulon <- pc[pc$TF %in% cluster.degs$gene & pc$target %in% cluster.degs$gene,]

  ## ±£´æclusterÕýµ÷¿ØµÄ½á¹û
  print("-- ## ±£´æclusterÕýµ÷¿ØµÄ½á¹û ## --")
  write.table(cluster.regulon,paste0(outdir1,"/01.Cluster",i,"_all_Positive_TF.Targets.txt"),row.names = F,sep = "\t",quote = F)

  tf_genes <- unique(cluster.regulon$TF)
  if (length(tf_genes) > 0) {

    # ½«TFÍ³¼Æ½á¹û´æ·ÅÔÚ¸Ã±äÁ¿ÖÐ
    tfSummary <- as.data.frame(matrix(0,length(tf_genes),4))
    m = 0

    ## Í³¼Æµ¥¸öTFµÄµ÷¿Ø¹ØÏµ
    for (j in tf_genes) {
      m <- m + 1
      outdir2 <- paste0(outdir1,"/",j)
      dir.create(outdir2,recursive = T)
      resTF <- cluster.regulon[cluster.regulon$TF == j,]

      # ¿É×÷ÎªcytoscapeµÄÊäÈë
      write.table(resTF,paste0(outdir2,"/1.",j,"_positive.regulon.txt"),row.names = F,sep = "\t",quote = F)

      # Í³¼ÆÃ¿¸öTFµÄËùÓÐ°Ð»ùÒò
      tfSummary[m,1] <- i
      tfSummary[m,2] <- j
      tfSummary[m,3] <- length(resTF$target)
      tfSummary[m,4] <- paste(as.vector(resTF$target),sep = "",collapse = ",")

      # ¼ÆËãTFºÍµ÷¿Ø»ùÒòµÄregulon·ÖÖµ£¨Ã¿¸öÏ¸°ûµÄGSZA·ÖÖµ£©
      print("-- ## ¼ÆËãTFºÍµ÷¿Ø»ùÒòµÄregulon·ÖÖµ£¨Ã¿¸öÏ¸°ûµÄGSZA·ÖÖµ£© ## --")
      tf.genes <- c(j,resTF$target)
      cpm <- as.data.frame(obj@assays$RNA@data[tf.genes,])
      z <- apply(cpm,1,scale)

      # TFºÍËùÓÐ°Ð»ùÒòµÄµ÷¿Ø·ÖÖµ£¨TFºÍ°Ð»ùÒòÖÐÎ»ÊýµÄÆ½¾ùÖµ£©
      tfRegulon <- as.data.frame(matrix(0,nrow = nrow(z),ncol = ncol(z)))
      for (q in 1:ncol(tfRegulon)) {
        if (q == 1 & ncol(tfRegulon) > 2) {
          tfRegulon[,q] <- as.data.frame((z[,1] + apply(z[,-1], 1, median))/2)
        }else if(q == 1 & ncol(tfRegulon) == 2){
          tfRegulon[,q] <- as.data.frame((z[,1] + z[,-1])/2)
        }else {
          tfRegulon[,q] <- as.data.frame((z[,1] + z[,q])/2)
        }
      }
      tfRegulon.out <- cbind(z[,1],tfRegulon)
      colnames(tfRegulon.out)[1] <- j
      write.table(tfRegulon,paste0(outdir2,"/2.",j,"_All_regulon.txt"),row.names = T,sep = "\t",quote = F)

      # »æÖÆËùÓÐºÍµ¥¸öTFµÄregulonÈ¾É«Í¼
      print("-- ## »æÖÆËùÓÐºÍµ¥¸öTFµÄregulonÈ¾É«Í¼ ## --")
      tfRegulon.out$cell <- rownames(tfRegulon.out)
      data_plot <- merge(umap,tfRegulon.out,by="cell")
      write.table(data_plot,paste0(outdir2,"/3.",j,"_All_regulon_UMAP.plot.txt"),row.names = T,sep = "\t",quote = F)

      outdir3 <- paste0(outdir2,"/4.TF_regulon_UMAP.Plot")
      dir.create(outdir3,recursive = T)
      for (h in 5:ncol(data_plot)) {
        tmp <- data_plot[,c(1,2,3,h)]
        name <- colnames(tmp)[4]
        colnames(tmp)[4] <- "GSZA"
        p <- ggplot(data=tmp, aes(x=UMAP_1, y=UMAP_2,color=GSZA))+
          geom_point(size = 0.5) +
          theme_bw()+ theme(panel.grid=element_blank(),panel.border=element_blank(),
                            axis.line=element_line(size=0.5,colour="black")) +
          scale_color_viridis_c()
        #scale_color_viridis_c(option = "plasma")
        #scale_color_gradient(low = "#e0e0e0", high = "#b2182b")
        ggsave(filename = paste0(outdir3,"/",name,".pdf"),width = 15,height = 15,p)
        ggsave(filename = paste0(outdir3,"/",name,".png"),width = 15,height = 15,dpi = 300,p)
      }

      ## ¼ÆËãÃ¿¸öclusterµÄ¾ùÖµ
      print("-- ## ¼ÆËãÃ¿¸öclusterµÄ¾ùÖµ ## --")
      tfRegulon$cell <- rownames(tfRegulon)
      tfRegulon.cluster <- merge(meta,tfRegulon,by="cell")
      cluster.mean <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster))+3,nrow = ncol(tfRegulon.cluster)-2))
      cluster.mean$Cluster <- i
      cluster.mean$TF <- j
      cluster.mean$name <- colnames(tfRegulon.cluster)[c(-1,-2)]
      for (x in unique(degs$cluster)) {
        y <- which(colnames(cluster.mean) == paste0("C",x))
        cluster.mean[,y] <- colMeans(tfRegulon.cluster[tfRegulon.cluster[,2] == x, c(-1,-2)])
      }

      ## ËùÓÐÊý¾Ý´æ·Å
      TF.mean <- rbind(TF.mean,cluster.mean[1,])
      TF.Regulon.mean <- rbind(TF.Regulon.mean,cluster.mean[-1,])

      ## clusterÖÐËùÓÐTFºÍ°Ð»ùÒòµÄÈÈÍ¼
      print("-- ## »æÖÆclusterÖÐËùÓÐTFºÍ°Ð»ùÒòµÄÈÈÍ¼ ## --")
      rownames(cluster.mean) <- cluster.mean$name
      c1 <- colorRamp2(breaks = c(-1,0,1),c("#2166ac","#FFFFFF","#b2182b"))
      h <- Heatmap(as.matrix(cluster.mean[,1:(ncol(TF.mean)-3)]), name ="TF Regulon",
                   column_title = " ",row_title = " ",
                   col = c1,
                   cluster_rows = F,cluster_columns = F,
                   show_column_names = T,show_row_names = T,
                   row_names_side = "left",
                   row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
                   column_dend_side = "top",column_title_rot = 0,
                   row_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
                   column_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
                   use_raster =F)
      wh <- nrow(cluster.mean[,1:(ncol(cluster.mean)-3)])/ncol(cluster.mean[,1:(ncol(cluster.mean)-3)])
      pdf(paste0(outdir2,"/5.",j,"_Positive.regulon_TF.heatmap.pdf"),width = 15 ,height = 15*wh)
      print(h)
      dev.off()

      ## »æÖÆClusterµÄÕÛÏßÍ¼£¨TFºÍ°Ð»ùÒò£©
      outdir4 <- paste0(outdir2,"/6.TF_regulon_line.Plot")
      dir.create(outdir4,recursive = T)

      tfRegulon.p <- as.data.frame(matrix(0,nrow = nrow(z),ncol = ncol(z)))
      colnames(tfRegulon.p)[1] <- "Alltargets"
      colnames(tfRegulon.p)[-1] <- colnames(z)[-1]
      for (q in 1:ncol(tfRegulon.p)) {
        if (q == 1 & ncol(tfRegulon.p) > 2) {
          tfRegulon.p[,q] <- apply(z[,-1], 1, median)
        }else if(q == 1 & ncol(tfRegulon.p) == 2){
          tfRegulon.p[,q] <-  z[,-1]
        }else {
          tfRegulon.p[,q] <-  z[,q]
        }
      }

      tfRegulon.p.out <- cbind(z[,1],tfRegulon.p)
      colnames(tfRegulon.p.out)[1] <- j
      tfRegulon.p.out$cell <- rownames(tfRegulon.p.out)

      tfRegulon.cluster.tf <- merge(meta,tfRegulon.p.out,by="cell")
      cluster.mean.tf <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster))+1,nrow = ncol(tfRegulon.cluster.tf)-2))
      cluster.mean.tf$name <- colnames(tfRegulon.cluster.tf)[c(-1,-2)]
      for (x in unique(degs$cluster)) {
        y <- which(colnames(cluster.mean.tf) == paste0("C",x))
        cluster.mean.tf[,y] <- colMeans(tfRegulon.cluster.tf[tfRegulon.cluster.tf[,2] == x, c(-1,-2)])
      }

      ## ±£´æÕÛÏßÍ¼µÄÊý¾Ý
      print("-- ## ±£´æÕÛÏßÍ¼µÄÊý¾Ý ## --")
      write.table(cluster.mean.tf,paste0(paste0(outdir4,"/",j,"_All.TF_Regulon.lineplot.data.txt")),row.names = T,sep = "\t",quote = F)

      ## »æÖÆÕÛÏßÍ¼
      print("-- ## »æÖÆÕÛÏßÍ¼ ## --")
      df.plot <- melt(cluster.mean.tf)
      df.plot$name <- factor(df.plot$name,levels = cluster.mean.tf$name)
      l1 <-  ggplot(data = df.plot, mapping = aes(x = variable, y = value, colour = name,group = name)) +
        geom_line(size=1) +
        geom_point(shape=21,size=1.75,fill="white") +
        xlab('Cluster') + ylab("Z score") +
        scale_color_manual(values =c("#b2182b","#2166ac",rep("#e0e0e0",length(cluster.mean.tf$name)-2))) +
        theme_bw() +
        theme(panel.grid=element_blank(),panel.border=element_blank(),
              axis.line=element_line(size=0.5,colour="black"))
      ggsave(paste0(outdir4,"/01.",j,"_All.TF_Regulon.lineplot.pdf"),width = 15,height = 15,l1)
      ggsave(paste0(outdir4,"/01.",j,"_All.TF_Regulon.lineplot.png"),width = 15,height = 15,dpi = 600,l1)

      ## ·Ö±ðÊä³öTFºÍÃ¿¸ö°Ð»ùÒòµÄÕÛÏßÍ¼
      print("-- ## ·Ö±ðÊä³öTFºÍÃ¿¸ö°Ð»ùÒòµÄÕÛÏßÍ¼ ## --")
      for (v in cluster.mean.tf$name[-1]) {
        tmp.plot <- df.plot[df.plot$name %in% c(j,v),]
        l <- ggplot(data = tmp.plot, mapping = aes(x = variable, y = value, colour = name,group = name)) +
          geom_line(size=1) +
          geom_point(shape=21,size=1.75,fill="white") +
          xlab('Cluster') + ylab("Z score") +
          scale_color_manual(values =c("#b2182b","#2166ac")) +
          theme_bw() +
          theme(panel.grid=element_blank(),panel.border=element_blank(),
                axis.line=element_line(size=0.5,colour="black"))
        ggsave(paste0(outdir4,"/",j,"_",v,"_Regulon.lineplot.pdf"),width = 15,height = 15,l)
        ggsave(paste0(outdir4,"/",j,"_",v,"_Regulon.lineplot.png"),width = 15,height = 15,dpi = 600,l)
      }

    }

    ## Êä³öÃ¿¸öclusterÖÐtfSummaryµÄ½á¹û£¨txtºÍ±í¸ñ£©
    print("-- ## Êä³öÃ¿¸öclusterÖÐTFºÍ°Ð»ùÒòµÄ½á¹û ## --")
    write.table(tfSummary,paste0(outdir1,"/02.Cluster_",i,"_TF_target_summary.txt"),row.names = F,sep = "\t",quote = F)
    write.xlsx(tfSummary,paste0(outdir1,"/03.Cluster_",i,"_TF_target_summary.xlsx"))

  }else{
    tfSummary <- paste0("No TFs and target genes in the cluster",i)
    write.xlsx(tfSummary,paste0(outdir1,"/03.Cluster_",i,"_TF_target_summary.xlsx"))
  }

  ## ÕûÀíËùÓÐclusterÕýµ÷¿ØµÄ½á¹û
  print("-- ## ÕûÀíËùÓÐclusterÕýµ÷¿ØµÄ½á¹û ## --")
  po.res[[paste0("Cluster",i)]] <- as.data.frame(tfSummary)

}

## ±£´æÕýµ÷¿ØÖÐTF regulonÃ¿¸öclusterµÄ¾ùÖµ£¨ÓÃÓÚ»æÖÆÈÈÍ¼£©
print("-- ##  ±£´æÕýµ÷¿ØÖÐTF regulonÃ¿¸öclusterµÄ¾ùÖµ ## --")
write.xlsx(po.res,paste0(outdir,"/Positive.regulon/Positive.regulon_TF_target_summary.xlsx"))

## »æÖÆTF regulonµÄÈÈÍ¼
### TFºÍËùÓÐ°Ð»ùÒò
print("-- ##  »æÖÆTF regulonµÄÈÈÍ¼ ## --")
TF.mean <- TF.mean[-1,]
TF.mean$type <- paste0("C",TF.mean$Cluster,"_",TF.mean$TF)
rownames(TF.mean) <- TF.mean$type
#rownames(TF.mean) <- gsub("\\.1","",rownames(TF.mean))
write.table(TF.mean,paste0(outdir,"/Positive.regulon/Positive.regulon_TF.heatmap-data.txt"),row.names = T,sep = "\t",quote = F)
c1 <- colorRamp2(breaks = c(-1,0,1),
                 c("#2166ac","#FFFFFF","#b2182b"))
p1 <- Heatmap(as.matrix(TF.mean[,1:(ncol(TF.mean)-4)]), name ="TF Regulon",
              column_title = " ",row_title = " ",
              col = c1,split = TF.mean$Cluster,
              cluster_rows = F,cluster_columns = F,
              show_column_names = T,show_row_names = T,
              row_names_side = "left",
              row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
              column_dend_side = "top",column_title_rot = 0,
              #top_annotation= annotation2,
              row_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              column_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              use_raster =F)
wh <- nrow(TF.mean[,1:(ncol(TF.mean)-4)])/ncol(TF.mean[,1:(ncol(TF.mean)-4)])
pdf(paste0(outdir,"/Positive.regulon/Positive.regulon_TF.heatmap.pdf"),width = 15 ,height = 40*wh)
p1
dev.off()

### TFºÍµ¥¸ö°Ð»ùÒò
TF.Regulon.mean <- TF.Regulon.mean[-1,]
TF.Regulon.mean$type <- paste0("C",TF.Regulon.mean$Cluster,"_",TF.Regulon.mean$name)
rownames(TF.Regulon.mean) <- TF.Regulon.mean$type
write.table(TF.Regulon.mean,paste0(outdir,"/Positive.regulon/Positive.regulon_TF.targets.heatmap-data.txt"),row.names = T,sep = "\t",quote = F)
c2 <- colorRamp2(breaks = c(-1,0,1),
                 c("#2166ac","#FFFFFF","#b2182b"))
p2 <- Heatmap(as.matrix(TF.Regulon.mean[,1:(ncol(TF.Regulon.mean)-4)]),
              name ="TF Regulon",
              column_title = " ",row_title = " ",
              col = c1,split = TF.Regulon.mean$Cluster,
              cluster_rows = F,cluster_columns = F,
              show_column_names = T,show_row_names = T,
              row_names_side = "left",
              row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
              column_dend_side = "top",column_title_rot = 0,
              #top_annotation= annotation2,
              row_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              column_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              use_raster =F)

wh <- nrow(TF.Regulon.mean[,1:(ncol(TF.Regulon.mean)-4)])/ncol(TF.Regulon.mean[,1:(ncol(TF.Regulon.mean)-4)])/4
pdf(paste0(outdir,"/Positive.regulon/Positive.regulon_TF.targets.heatmap.pdf"),width = 15 ,height = 40*wh)
p2
dev.off()

# 7. Í³¼ÆÃ¿¸öclusterÖÐµÄ¸ºµ÷¿Ø -------------------------------------------------
print("** 7. Í³¼ÆÃ¿¸öclusterÖÐµÄ¸ºµ÷¿Ø **")
TF.mean.Neg <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster)) *2+3,nrow = 1))

TF.Regulon.mean.Neg <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster)) *2+3,nrow = 1))

neg.res <- list()
for (i in unique(degs$cluster)) {
  print(paste0("+++++++++++++++++++++ ## Run ",i," ## +++++++++++++++++++++"))
  ## Ã¿¸öclusterµ¥¶À´´½¨Ä¿Â¼
  outdir1 <- paste0(outdir,"/Negative.regulon/Cluster",i)
  dir.create(outdir1,recursive = T)

  cluster.degs <- degs[degs$cluster == i ,]
  cluster.regulon.tmp <- nc[nc$TF %in% cluster.degs$gene,]
  cluster.regulon <- cluster.regulon.tmp[!cluster.regulon.tmp$target %in% cluster.degs$gene,]

  ## ±£´æcluster¸ºµ÷¿ØµÄ½á¹û
  print("-- ##  »æÖÆTF regulonµÄÈÈÍ¼ ## --")
  write.table(cluster.regulon,paste0(outdir1,"/01.Cluster",i,"_all_Negative_TF.Targets.txt"),row.names = F,sep = "\t",quote = F)

  tf_genes <- unique(cluster.regulon$TF)
  if (length(tf_genes) > 0) {

    # ½«TFÍ³¼Æ½á¹û´æ·ÅÔÚ¸Ã±äÁ¿ÖÐ
    tfSummary <- as.data.frame(matrix(0,length(tf_genes),4))
    m = 0

    ## Í³¼Æµ¥¸öTFµÄµ÷¿Ø¹ØÏµ
    print("-- ## Í³¼Æµ¥¸öTFµÄµ÷¿Ø¹ØÏµ ## --")
    for (j in tf_genes) {
      m <- m + 1
      outdir2 <- paste0(outdir1,"/",j)
      dir.create(outdir2,recursive = T)
      resTF <- cluster.regulon[cluster.regulon$TF == j,]

      # ¿É×÷ÎªcytoscapeµÄÊäÈë
      write.table(resTF,paste0(outdir2,"/1.",j,"_negative.regulon.txt"),row.names = F,sep = "\t",quote = F)

      # Í³¼ÆÃ¿¸öTFµÄËùÓÐ°Ð»ùÒò
      print("-- ## Í³¼ÆÃ¿¸öTFµÄËùÓÐ°Ð»ùÒò ## --")
      tfSummary[m,1] <- i
      tfSummary[m,2] <- j
      tfSummary[m,3] <- length(resTF$target)
      tfSummary[m,4] <- paste(as.vector(resTF$target),sep = "",collapse = ",")

      # ¼ÆËãTFºÍµ÷¿Ø»ùÒòµÄregulon·ÖÖµ£¨Ã¿¸öÏ¸°ûµÄGSZA·ÖÖµ£©
      print("-- ## ¼ÆËãTFºÍµ÷¿Ø»ùÒòµÄregulon·ÖÖµ ## --")
      tf.genes <- c(j,resTF$target)
      data <- as.data.frame(obj@assays$RNA@data[tf.genes,])
      z <- apply(data,1,scale)

      # TFºÍËùÓÐ°Ð»ùÒòµÄµ÷¿Ø·ÖÖµ£¨TFºÍ°Ð»ùÒòÖÐÎ»ÊýµÄÆ½¾ùÖµ£©
      print("-- ## TFºÍËùÓÐ°Ð»ùÒòµÄµ÷¿Ø·ÖÖµ ## --")
      tfRegulon <- as.data.frame(matrix(0,nrow = nrow(z),ncol = ncol(z)))
      for (q in 1:ncol(tfRegulon)) {
        if (q == 1 & ncol(tfRegulon) > 2) {
          tfRegulon[,q] <- apply(z[,-1], 1, median)
        }else if(q == 1 & ncol(tfRegulon) == 2){
          tfRegulon[,q] <-  z[,-1]
        }else {
          tfRegulon[,q] <-  z[,q]
        }
      }

      tfRegulon.out <- cbind(z[,1],tfRegulon)
      colnames(tfRegulon.out)[1] <- j
      write.table(tfRegulon.out,paste0(outdir2,"/2.",j,"_All_regulon.txt"),row.names = T,sep = "\t",quote = F)

      # »æÖÆËùÓÐºÍµ¥¸öTFµÄregulonÈ¾É«Í¼
      print("-- ## »æÖÆËùÓÐºÍµ¥¸öTFµÄregulonÈ¾É«Í¼ ## --")
      tfRegulon.out$cell <- rownames(tfRegulon.out)
      data_plot <- merge(umap,tfRegulon.out,by="cell")
      write.table(data_plot,paste0(outdir2,"/3.",j,"_All_regulon_UMAP.plot.txt"),row.names = T,sep = "\t",quote = F)

      outdir3 <- paste0(outdir2,"/4.TF_regulon_Plot")
      dir.create(outdir3,recursive = T)
      tmp <- data_plot[,c(1,2,3,4)]
      name <- colnames(tmp)[4]
      colnames(tmp)[4] <- "TF"
      p1 <- ggplot(data=tmp, aes(x=UMAP_1, y=UMAP_2,color=TF))+
        geom_point(size = 0.5) +
        theme_bw()+ theme(panel.grid=element_blank(),panel.border=element_blank(),
                          axis.line=element_line(size=0.5,colour="black")) +
        scale_color_viridis_c()

      for (h in 5:ncol(data_plot)) {
        tmp <- data_plot[,c(1,2,3,h)]
        name <- colnames(tmp)[4]
        colnames(tmp)[4] <- "GSZA"
        p <- ggplot(data=tmp, aes(x=UMAP_1, y=UMAP_2,color=GSZA))+
          geom_point(size = 0.5) +
          theme_bw()+ theme(panel.grid=element_blank(),panel.border=element_blank(),
                            axis.line=element_line(size=0.5,colour="black")) +
          scale_color_viridis_c()
        #scale_color_viridis_c(option = "plasma")
        #scale_color_gradient(low = "#e0e0e0", high = "#b2182b")
        ggsave(filename = paste0(outdir3,"/",name,".pdf"),width = 15,height = 15,p1 + p)
        ggsave(filename = paste0(outdir3,"/",name,".png"),width = 15,height = 15,dpi = 300,p1 + p)
      }

      ## ¼ÆËãÃ¿¸öclusterµÄ¾ùÖµ
      print("-- ## ¼ÆËãÃ¿¸öclusterµÄ¾ùÖµ ## --")
      tfRegulon.cluster <- merge(meta,tfRegulon.out,by="cell")
      cluster.mean <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster)) *2+3,nrow = ncol(tfRegulon.cluster)-3))
      cluster.mean$Cluster <- i
      cluster.mean$TF <- j
      cluster.mean$name <- colnames(tfRegulon.cluster)[c(-1,-2,-3)]
      for (x in unique(degs$cluster)) {
        y <- which(colnames(cluster.mean) == paste0(x,".TF"))
        cluster.mean[,y] <- mean(tfRegulon.cluster[tfRegulon.cluster[,2] == x, 3])
        cluster.mean[,y+1] <- colMeans(tfRegulon.cluster[tfRegulon.cluster[,2] == x, c(-1,-2,-3)])
      }

      ## ËùÓÐÊý¾Ý´æ·Å
      print("-- ## ËùÓÐÊý¾Ý´æ·Å ## --")
      TF.mean.Neg <- rbind(TF.mean.Neg,cluster.mean[1,])
      TF.Regulon.mean.Neg <- rbind(TF.Regulon.mean.Neg,cluster.mean[-1,])

      ## clusterÖÐËùÓÐTFºÍ°Ð»ùÒòµÄÈÈÍ¼
      print("-- ## »æÖÆclusterÖÐËùÓÐTFºÍ°Ð»ùÒòµÄÈÈÍ¼ ## --")
      rownames(cluster.mean) <- cluster.mean$name
      c1 <- colorRamp2(breaks = c(-1,0,1),c("#2166ac","#FFFFFF","#b2182b"))
      h <- Heatmap(as.matrix(cluster.mean[,1:(ncol(TF.mean.Neg)-3)]), name ="TF Regulon",
                   column_title = " ",row_title = " ",
                   col = c1,
                   cluster_rows = F,cluster_columns = F,
                   show_column_names = T,show_row_names = T,
                   row_names_side = "left",
                   row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
                   column_dend_side = "top",column_title_rot = 0,
                   row_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
                   column_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
                   use_raster =F)
      wh <- nrow(cluster.mean[,1:(ncol(cluster.mean)-3)])/ncol(cluster.mean[,1:(ncol(cluster.mean)-3)])
      pdf(paste0(outdir2,"/5.",j,"_Negative.regulon_TF.heatmap.pdf"),width = 15 ,height = 15*wh)
      print(h)
      dev.off()

      ## »æÖÆcluster»ùÒòµÄÕÛÏßÍ¼
      print("-- ## »æÖÆcluster»ùÒòµÄÕÛÏßÍ¼ ## --")
      outdir4 <- paste0(outdir2,"/6.TF_regulon_line.Plot")
      dir.create(outdir4,recursive = T)

      cluster.mean.tf.neg <- as.data.frame(matrix(data = 0,ncol = length(unique(degs$cluster))+1,nrow = ncol(tfRegulon.cluster)-2))
      cluster.mean.tf.neg$name <- colnames(tfRegulon.cluster)[c(-1,-2)]
      cluster.mean.tf.neg$name <- gsub(paste0(j,"_"),"",cluster.mean.tf.neg$name)
      cluster.mean.tf.neg[1,1:length(unique(degs$cluster))] <- cluster.mean[1,seq(1,length(unique(degs$cluster)) *2,by=2)]
      cluster.mean.tf.neg[-1,1:length(unique(degs$cluster))] <- cluster.mean[,seq(2,length(unique(degs$cluster)) *2,by=2)]

      ## »æÖÆÕÛÏßÍ¼
      print("-- ## »æÖÆÕÛÏßÍ¼ ## --")
      df.plot.neg <- melt(cluster.mean.tf.neg)
      df.plot.neg$name <- factor(df.plot.neg$name,levels = cluster.mean.tf.neg$name)
      l1 <-  ggplot(data = df.plot.neg, mapping = aes(x = variable, y = value, colour = name,group = name)) +
        geom_line(size=1) +
        geom_point(shape=21,size=1.75,fill="white") +
        xlab('Cluster') + ylab("Z score") +
        scale_color_manual(values =c("#b2182b","#1a9850",rep("#e0e0e0",length(cluster.mean.tf.neg$name)-2))) +
        theme_bw() +
        theme(panel.grid=element_blank(),panel.border=element_blank(),
              axis.line=element_line(size=0.5,colour="black"))
      ggsave(paste0(outdir4,"/01.",j,"_All.TF_Regulon.lineplot.pdf"),width = 15,height = 15,l1)
      ggsave(paste0(outdir4,"/01.",j,"_All.TF_Regulon.lineplot.png"),width = 15,height = 15,dpi = 600,l1)

      ## ·Ö±ðÊä³öTFºÍÃ¿¸ö°Ð»ùÒòµÄÕÛÏßÍ¼
      print("-- ## ·Ö±ðÊä³öTFºÍÃ¿¸ö°Ð»ùÒòµÄÕÛÏßÍ¼ ## --")
      for (v in cluster.mean.tf.neg$name[-1]) {
        tmp.plot <- df.plot.neg[df.plot.neg$name %in% c(j,v),]
        l <- ggplot(data = tmp.plot, mapping = aes(x = variable, y = value, colour = name,group = name)) +
          geom_line(size=1) +
          geom_point(shape=21,size=1.75,fill="white") +
          xlab('Cluster') + ylab("Z score") +
          scale_color_manual(values =c("#b2182b","#1a9850")) +
          theme_bw() +
          theme(panel.grid=element_blank(),panel.border=element_blank(),
                axis.line=element_line(size=0.5,colour="black"))
        ggsave(paste0(outdir4,"/",j,"_",v,"_Regulon.lineplot.pdf"),width = 15,height = 15,l)
        ggsave(paste0(outdir4,"/",j,"_",v,"_Regulon.lineplot.png"),width = 15,height = 15,dpi = 600,l)
      }

    }

    ## Êä³öÃ¿¸öclusterÖÐtfSummaryµÄ½á¹û£¨txtºÍ±í¸ñ£©
    print("-- ## Êä³öÃ¿¸öclusterÖÐtfSummaryµÄ½á¹û ## --")
    write.table(tfSummary,paste0(outdir1,"/02.Cluster_",i,"_TF_target_summary.txt"),row.names = F,sep = "\t",quote = F)
    write.xlsx(tfSummary,paste0(outdir1,"/03.Cluster_",i,"_TF_target_summary.xlsx"))

  }else{
    tfSummary <- paste0("No TFs and target genes in the cluster",i)
    write.xlsx(tfSummary,paste0(outdir1,"/03.Cluster_",i,"_TF_target_summary.xlsx"))
  }

  ## ÕûÀíËùÓÐcluster¸ºµ÷¿ØµÄ½á¹û
  neg.res[[paste0("Cluster",i)]] <- as.data.frame(tfSummary)

}

## ±£´æÕýµ÷¿ØÖÐTF regulonÃ¿¸öclusterµÄ¾ùÖµ£¨ÓÃÓÚ»æÖÆÈÈÍ¼£©
write.xlsx(neg.res,paste0(outdir,"/Negative.regulon/Negative.regulon_TF_target_summary.xlsx"))

## »æÖÆTF regulonµÄÈÈÍ¼
### TFºÍËùÓÐ°Ð»ùÒò
TF.mean.Neg <- TF.mean.Neg[-1,]
TF.mean.Neg$type <- paste0("C",TF.mean.Neg$Cluster,"_",TF.mean.Neg$TF)
rownames(TF.mean.Neg) <- TF.mean.Neg$type
write.table(TF.mean.Neg,paste0(outdir,"/Negative.regulon/Negative.regulon_TF.heatmap-data.txt"),row.names = T,sep = "\t",quote = F)
c1 <- colorRamp2(breaks = c(-1,0,1),
                 c("#2166ac","#FFFFFF","#b2182b"))
p1 <- Heatmap(as.matrix(TF.mean.Neg[,1:(ncol(TF.mean.Neg)-4)]), name ="TF Regulon",
              column_title = " ",row_title = " ",
              col = c1,split = TF.mean.Neg$Cluster,
              cluster_rows = F,cluster_columns = F,
              show_column_names = T,show_row_names = T,
              row_names_side = "left",
              row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
              column_dend_side = "top",column_title_rot = 0,
              #top_annotation= annotation2,
              row_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              column_dend_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              use_raster =F)
wh <- nrow(TF.mean.Neg[,1:(ncol(TF.mean.Neg)-4)])/ncol(TF.mean.Neg[,1:(ncol(TF.mean.Neg)-4)])
pdf(paste0(outdir,"/Negative.regulon/Negative.regulon_TF.heatmap.pdf"),width = 15 ,height = 40*wh)
p1
dev.off()

### TFºÍµ¥¸ö°Ð»ùÒò
TF.Regulon.mean.Neg <- TF.Regulon.mean.Neg[-1,]
TF.Regulon.mean.Neg$type <- paste0(TF.Regulon.mean.Neg$Cluster,"_",TF.Regulon.mean.Neg$name)
rownames(TF.Regulon.mean.Neg) <- TF.Regulon.mean.Neg$type
write.table(TF.Regulon.mean.Neg,paste0(outdir,"/Negative.regulon/Negative.regulon_TF.targets.heatmap-data.txt"),row.names = T,sep = "\t",quote = F)
c2 <- colorRamp2(breaks = c(-1,0,1),
                 c("#2166ac","#FFFFFF","#b2182b"))
p2 <- Heatmap(as.matrix(TF.Regulon.mean.Neg[,1:(ncol(TF.Regulon.mean.Neg)-4)]),
              name ="TF Regulon",
              column_title = " ",row_title = " ",
              col = c1,split = TF.Regulon.mean.Neg$Cluster,
              cluster_rows = F,cluster_columns = F,
              show_column_names = T,show_row_names = T,
              row_names_side = "left",
              row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
              column_dend_side = "top",column_title_rot = 0,
              #top_annotation= annotation2,
              row_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              column_dend_gp = gpar(fontsize = 10, fontface="bold", col="black"),
              use_raster =F)

wh <- nrow(TF.Regulon.mean.Neg[,1:(ncol(TF.Regulon.mean.Neg)-4)])/ncol(TF.Regulon.mean.Neg[,1:(ncol(TF.Regulon.mean.Neg)-4)])/2
pdf(paste0(outdir,"/Negative.regulon/Negative.regulon_TF.targets.heatmap.pdf"),width = 15 ,height = 40*wh)
p2
dev.off()
