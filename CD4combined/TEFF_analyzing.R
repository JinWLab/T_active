# ------------- CD4combined_scRNAseq_TEFF_analyzing -------------
# This code shows how to explore difference of different effector T cells subsets.

# TEFF specific genes
TEFF <- subset(cd4.combined,idents = c('cytohi TEFF','CTLA4hi TEFF','Conv TEFF','Prolif TEFF'))
saveRDS(TEFF, file = 'TEFF.rds')

# FigS3
#CTLA4hi Th1
Dat <- FindMarkers(object = TEFF , ident.1 = 'CTLA4hi Th1', ident.2 = 'Conv TEFF',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 1, ifelse(Dat$avg_log2FC>= 1 ,'CTLA4hi Th1','Conv TEFF'),'NoSignifi'),levels=c('CTLA4hi Th1','Conv TEFF','NoSignifi'))
repel.data <- Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>1,]
pdf('FigS3-A1.pdf',width = 5,height = 3)
ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>1,],
    aes(label = allele),
    size = 3,
    segment.color = "black", show.legend = FALSE,
    max.overlaps = Inf)+#gene names of the point of interest
  theme_bw()+#
  theme(
    legend.title = element_blank()#no legend title
  )+
  ylab('-log10 (p-adj)')+#y-axis name
  xlab('log2 (FoldChange)')+#x-axis name
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+  #添加竖线p_val_adj<0.05
  ggtitle(label = 'CTLA4hi Th1 vs Conv TEFF' )
dev.off()
write.csv(Dat,'CTLA4hi Th1_Dat.csv')

#Th1
Dat <- FindMarkers(object = TEFF , ident.1 = 'Th1', ident.2 = 'Conv TEFF',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 1, ifelse(Dat$avg_log2FC>= 1 ,'Th1','Conv TEFF'),'NoSignifi'),levels=c('Th1','Conv TEFF','NoSignifi'))
repel.data <- Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>1,]
pdf('FigS3-A2.pdf',width = 5,height = 3)
ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>1,],
    aes(label = allele),
    size = 3,
    segment.color = "black", show.legend = FALSE,
    max.overlaps = Inf)+#
  theme_bw()+#Modify background
  theme(
    legend.title = element_blank()#no legend title
  )+
  ylab('-log10 (p-adj)')+#y-axis name
  xlab('log2 (FoldChange)')+#x-axis name
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+  #添加竖线p_val_adj<0.05
  ggtitle(label = 'Th1 vs Conv TEFF' )
dev.off()
write.csv(Dat,'Th1_Dat.csv')

#Prolif
Dat <- FindMarkers(object = TEFF , ident.1 = 'Prolif TEFF', ident.2 = 'Conv TEFF',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 1, ifelse(Dat$avg_log2FC>= 1 ,'Prolif TEFF','Conv TEFF'),'NoSignifi'),levels=c('Prolif TEFF','Conv TEFF','NoSignifi'))
repel.data <- Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>1,]
pdf('FigS3-C.pdf',width = 5,height = 3)
ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>1,],
    aes(label = allele),
    size = 3,
    segment.color = "black", show.legend = FALSE,
    max.overlaps = Inf)+#gene names of the point of interest
  theme_bw()+#Modify background
  theme(
    legend.title = element_blank()#no legend title
  )+
  ylab('-log10 (p-adj)')+#y-axis name
  xlab('log2 (FoldChange)')+#x-axis name
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+  #添加竖线p_val_adj<0.05
  ggtitle(label = 'Prolif TEFF vs Conv TEFF' )
dev.off()
write.csv(Dat,'Prolif_Dat.csv')