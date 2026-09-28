# Portable path configuration (override with environment variables).
project_dir <- Sys.getenv("SHOOT_REGEN_PROJECT_DIR", unset = ".")
output_root <- Sys.getenv("SHOOT_REGEN_OUTPUT_DIR", unset = file.path(project_dir, "02_results"))
external_data_dir <- Sys.getenv("SHOOT_REGEN_EXTERNAL_DATA_DIR", unset = file.path(project_dir, "external_data"))
work_dir <- Sys.getenv("SHOOT_REGEN_WORK_DIR", unset = ".")

suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(scales)})
input_rds<-file.path(project_dir, "01_obj/01_WT_snRNA/03_snRNA_C0-2-3-subclustered_res0.4_SIM4-SIM16.rds")
outdir<-Sys.getenv("PROPORTION_OUTDIR",file.path(output_root, "123_WT_res04_C02418_proportions_by_time"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
targets<-c("0","2","4","1","8");dates<-c("SIM4","SIM6","SIM8","SIM10","SIM14","SIM16");date_labs<-setNames(c("SIM4","SIM6","SIM8","SIM10","SIM14","SIM16"),dates)
# Colors are taken directly from subset_cell_ratio_change_barplot.R.
cols<-c("0"="#EBCC50","1"="#EBAB50","2"="#EB8950","4"="#EB505A","8"="#86431C")
obj<-readRDS(input_rds);md<-obj@meta.data;md$subcluster<-as.character(md$subcluster);md$date<-factor(as.character(md$date),levels=dates)
alltab<-as.data.frame(table(date=md$date,subcluster=md$subcluster));tot<-aggregate(Freq~date,alltab,sum);names(tot)[2]<-"total_cells"
d<-merge(alltab[alltab$subcluster%in%targets,],tot,by="date",all.x=TRUE);d$subcluster<-factor(d$subcluster,levels=targets);d$date<-factor(d$date,levels=dates);d$proportion_all_cells<-d$Freq/d$total_cells;d<-d[order(d$date,d$subcluster),]
write.csv(d,file.path(outdir,"02_C02418_counts_and_proportions_by_time.csv"),row.names=FALSE)
sumtab<-aggregate(cbind(selected_cells=Freq,proportion_all_cells=proportion_all_cells)~date,d,sum);write.csv(sumtab,file.path(outdir,"03_C02418_combined_fraction_by_time.csv"),row.names=FALSE)
theme_pub<-theme_classic(base_family="Arial",base_size=10)+theme(axis.text=element_text(size=9,colour="black"),axis.title=element_text(size=10,face="bold"),legend.title=element_text(size=9,face="bold"),legend.text=element_text(size=8.5),plot.margin=margin(5,6,5,5))
save3<-function(p,stem,w=4.7,h=3.45){ggsave(file.path(outdir,paste0(stem,".pdf")),p,width=w,height=h,device=cairo_pdf);ggsave(file.path(outdir,paste0(stem,".svg")),p,width=w,height=h,device=grDevices::svg);ggsave(file.path(outdir,paste0(stem,".png")),p,width=w,height=h,dpi=500,bg="white")}
pbar<-ggplot(d,aes(date,proportion_all_cells,fill=subcluster))+geom_col(position=position_stack(reverse=TRUE),width=.7,colour="white",linewidth=.3)+scale_fill_manual(values=cols,breaks=targets,drop=FALSE,name="Subcluster")+scale_x_discrete(labels=date_labs)+scale_y_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0,.03)))+labs(x=NULL,y="Proportion of all cells")+theme_pub+theme(panel.grid.major.y=element_line(colour="#ECECEC",linewidth=.28))
save3(pbar,"04_C02418_stacked_proportion_by_time")
plab<-pbar+geom_text(aes(label=ifelse(proportion_all_cells>=.025,percent(proportion_all_cells,accuracy=.1),"")),position=position_stack(reverse=TRUE,vjust=.5),size=2.55,family="Arial",colour="white",fontface="bold")
save3(plab,"05_C02418_stacked_proportion_by_time_labeled")
pline<-ggplot(d,aes(date,proportion_all_cells,group=subcluster,colour=subcluster))+geom_line(linewidth=.8)+geom_point(aes(fill=subcluster),shape=21,colour="#333333",stroke=.3,size=2.35)+scale_colour_manual(values=cols,breaks=targets,drop=FALSE,name="Subcluster")+scale_fill_manual(values=cols,breaks=targets,drop=FALSE)+scale_x_discrete(labels=date_labs)+scale_y_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0,.04)))+guides(fill="none")+labs(x=NULL,y="Proportion of all cells")+theme_pub+theme(panel.grid.major.y=element_line(colour="#ECECEC",linewidth=.28),legend.position="top")
save3(pline,"06_C02418_multiline_proportion_by_time",5.1,3.55)
pfacet<-ggplot(d,aes(date,proportion_all_cells,group=subcluster))+geom_line(colour="#555555",linewidth=.65)+geom_point(aes(fill=subcluster),shape=21,colour="#303030",stroke=.3,size=2.2)+facet_wrap(~subcluster,ncol=3,scales="free_y",drop=FALSE)+scale_fill_manual(values=cols,guide="none")+scale_x_discrete(labels=date_labs)+scale_y_continuous(labels=percent_format(accuracy=.1),expand=expansion(mult=c(.04,.1)))+labs(x=NULL,y="Proportion of all cells")+theme_pub+theme(strip.background=element_rect(fill="#F3F3F3",colour="#D0D0D0",linewidth=.35),strip.text=element_text(face="bold",size=9.5),axis.text.x=element_text(angle=45,hjust=1,size=7.5),axis.text.y=element_text(size=7.5),panel.grid.major.y=element_line(colour="#ECECEC",linewidth=.25),panel.spacing=grid::unit(.55,"lines"))
save3(pfacet,"07_C02418_faceted_proportion_by_time",5.6,4.5)
writeLines(c("WT res.0.4 subclusters 0,2,4,1,8 proportions across D4-D16.","The denominator is all nuclei retained at each time point; selected stacks therefore need not sum to 100%.","Colors follow <PROJECT_DIR>/03_script/subset_cell_ratio_change_barplot.R.",paste("Displayed order:",paste(targets,collapse=" -> "))),file.path(outdir,"08_analysis_notes.txt"))
cat("Completed:",outdir,"\n")
