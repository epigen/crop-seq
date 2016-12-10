library(data.table)
library(ggplot2)
library(reshape2)
theme_set(theme_bw())


data_dir="/scratch/lab_bock/shared/projects/crop-seq/results_pipeline/Drop-seq_HEK293T-3T3/"
out_dir="/scratch/lab_bock/shared/projects/crop-seq/mixing_hum_mouse"


dir.create(out_dir)

refFlat=fread("/data/groups/lab_bock/shared/resources/genomes/hg19_mm10_transgenes/hg19_mm10_transgenes.refFlat",drop=c(7:11))
setnames(refFlat,names(refFlat),c("GENE","transcript","chrom","strand","start","end"))
refFlat=refFlat[!grepl("ERCC",chrom),]

umiPerGene=fread(file.path(data_dir,"cell_umi_barcodes.500genes.tsv"))
dge=fread(file.path(data_dir,"digital_expression.500genes.tsv"))



#plot number of duplicate UMIs per cell
duplUMIs=umiPerGene[Num_Obs>1,.N,by="Cell Barcode"]
duplUMIs[,N_topped:=ifelse(N>1000,1000,N),]
duplScore=nrow(umiPerGene[Num_Obs>1])/nrow(umiPerGene)
pdf(paste0(out_dir,"/duplUMIs_per_cell_hist.pdf"),height=6,width=6)
ggplot(duplUMIs,aes(x=log10(N)))+geom_histogram(bins=40,col="black",fill="grey",size=0.2)+xlab("Number of duplicate UMIs (log10)")+ylab("Number of cells")+ggtitle(paste0("Duplication rate: ",round(duplScore,3)))
dev.off()


#plot percentage of duplicate UMIs per cell
percent_dup=umiPerGene[,.N,by=c("Num_Obs","Cell Barcode")][,list(percent_dupl=sum(N[Num_Obs>1]/sum(N)*100)),by="Cell Barcode"]
pdf(paste0(out_dir,"/duplUMIs_percent_cell_hist.pdf"),height=6,width=6)
ggplot(percent_dup,aes(x=percent_dupl))+geom_histogram(bins=50,col="black",fill="grey",size=0.2)+xlab("Percent duplicate UMIs")+ylab("Number of cells")+ggtitle(paste0("Duplication rate: ",round(duplScore,3)))
dev.off()


#actually analyze the expression data
dge=merge(refFlat,dge,by="GENE")
dge_long=melt(dge,id.vars=c("GENE","transcript","chrom","strand","start","end"),value.name="count",variable.name="cell_barcode")
dge_long[,species:=sub("_.*","",chrom),]

dge_long_uniq=unique(dge_long[,-grep("transcript|chrom|strand|start|end",names(dge_long)),with=FALSE]) 

dge_long_mean=dge_long_uniq[,list(g1=sum(count>=1),g5=sum(count>=5),g10=sum(count>=10),mean=mean(count[count>0])),by=c("cell_barcode","species")]

dge_long_mean[,diff_mh_g1:=g1[species=="MOUSE"]-g1[species=="HUMAN"],by="cell_barcode"]
dge_long_mean[,ratio_mh_g1:=g1[species=="MOUSE"]/g1[species=="HUMAN"],by="cell_barcode"]
dge_long_mean[,percent_g1:=g1/sum(g1)*100,by="cell_barcode"]
dge_long_mean[,total_g1:=sum(g1),by="cell_barcode"]
dge_long_mean[,species_dom:=ifelse(diff_mh_g1>0,"MOUSE",ifelse(diff_mh_g1<0,"HUMAN","tie")),]

dge_wide_mean=reshape(dge_long_mean,idvar="cell_barcode",timevar="species",direction="wide")
dge_wide_mean[,doublet:=ifelse(percent_g1.MOUSE>20&percent_g1.MOUSE<80,TRUE,FALSE),]
double_frac=(nrow(dge_wide_mean[percent_g1.MOUSE>20&percent_g1.MOUSE<80])*2)/nrow(dge_wide_mean)


#ratio plots
dge_long_mean[,cell_barcode:=factor(cell_barcode,levels=unique(cell_barcode[order(ratio_mh_g1)])),]
pdf(paste0(out_dir,"/species_ratio_g1.pdf"),height=7,width=5)
ggplot(dge_long_mean[percent_g1!=0], aes(x=cell_barcode,y=percent_g1,fill=species,alpha=log10(total_g1)))+geom_bar(position="stack",stat="identity",col="transparent")+ggtitle(paste0("doublet-estimate: ",round(double_frac,3)*100,"%"))+theme(axis.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+coord_flip()+xlab("Cells")+ylab("% of total detected genes")
dev.off()

#mouse vs. human
pdf(paste0(out_dir,"/species_cor_g1.pdf"),height=6,width=5)
ggplot(dge_wide_mean, aes(x=log10(g1.MOUSE),y=log10(g1.HUMAN),fill=species_dom.MOUSE,col=doublet))+geom_abline(slope=1)+geom_point(shape=21,alpha=0.4,size=2)+theme(legend.position="bottom")+ggtitle(paste0("doublet-estimate: ",round(double_frac,3)*100,"%"))+xlim(c(0,4))+ylim(c(0,4))+scale_color_manual(values=c("TRUE"="black","FALSE"="grey"))
dev.off()

pdf(paste0(out_dir,"/species_cor_g1_noLog.pdf"),height=6,width=5)
ggplot(dge_wide_mean, aes(x=g1.MOUSE,y=g1.HUMAN,fill=species_dom.MOUSE,col=doublet))+geom_abline(slope=1)+geom_point(shape=21,alpha=0.4,size=2)+theme(legend.position="bottom")+ggtitle(paste0("doublet-estimate: ",round(double_frac,3)*100,"%"))+scale_color_manual(values=c("TRUE"="black","FALSE"="grey"))
dev.off()

