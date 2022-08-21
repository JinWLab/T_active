# ------------- CD4combined_scRNAseq_HSPhiT -------------
# This code shows the analyzing of HSPhi T cells, including expression differences 
# of HSP specific genes, and re-reduction and reclassification of cells by using 
# or not using HSP specific genes.  

# reclustering
DefaultAssay(cd4.combined) <- "integrated"
cd4.combined <- ScaleData(cd4.combined, verbose = FALSE)
cd4.combined <- RunPCA(cd4.combined, npcs = 30, verbose = FALSE)
cd4.combined <- RunUMAP(cd4.combined, reduction = "pca", dims = 1:20, return.model = TRUE)
DefaultAssay(cd4.combined) <- 'RNA'
DimPlot(cd4.combined, reduction = "umap", label = F, cols = set3)

# Take hsphi and other cells separately
cd4.combined$celltype <- cd4.combined@active.ident
other <- subset(cd4.combined, idents =  c('TN','TCM','TEM','ISAGhi','TEMRA','Treg','cytohi TEFF','CTLA4hi TEFF','Conv TEFF','Prolif TEFF'))
hsp <- subset(cd4.combined, idents = 'HSPhi')
DimPlot(other, reduction = "umap", label = F, cols = set3)
DimPlot(hsp, reduction = "umap", label = F, cols = set3)

# Map hsp to other
Dat <- FindMarkers(object = cd4.combined , ident.1 = 'HSPhi',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 2, ifelse(Dat$avg_log2FC>= 2 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Dat[Dat$threshold == 'Up',] %>% write.csv(file = 'gene_hsphi.csv')

hsp.gene <- Dat[Dat$threshold == 'Up',]$allele
other <- FindVariableFeatures(other, selection.method = "vst", nfeatures = 2000, assay = 'RNA')
used.gene <- setdiff(VariableFeatures(other),intersect(hsp.gene,VariableFeatures(other)))  # Delete hsphi genes to rerun umap

anchors <- FindTransferAnchors(reference = other,query = hsp,reference.reduction = "pca", dims = 1:20,features = used.gene)
hsp <- MapQuery(anchorset = anchors,query = hsp,reference = other,refdata = list(celltype = "celltype",predicted_ADT = "ADT"),reference.reduction = "pca", reduction.model = "umap")
DimPlot(hsp, reduction = "ref.umap", group.by = 'predicted.celltype',label = F, repel = TRUE)# + NoLegend()

saveRDS(other, file = 'other.rds')
saveRDS(hsp, file = 'hsp.rds')

#Re-merge to get hsp.combined
hsp@reductions[["umap"]] <- hsp@reductions[["ref.umap"]]
hsp.combined <- merge(other, hsp, merge.data = T, merge.dr = 'umap')
hsp.combined@meta.data[which(hsp.combined@active.ident == 'HSPhi'),'hsp'] <- 'HSPhi'
hsp.combined@meta.data[which(hsp.combined@active.ident != 'HSPhi'),'hsp'] <- 'Other'
# Idents(hsp.combined) <- factor(Idents(hsp.combined), levels = c("TN","TCM","TEM","TEMRA","Treg","cytohi TEFF","Conv TEFF","Prolif TEFF",'HSPhi'))


# Dot plot of HSP genes
genelist <- c('HSPA1A','HSPA1B','HSPA6','HSPB1','HSPA2','DNAJB1','DNAJB4')
genelist <- c('HSPA1A','HSPA1B','HSPA6','HSPB1','HSPA2','DNAJB1','DNAJB4','HSP90AA1','HSP90AB1','HSPH1','HSPD1','HSPA8','HSPBP1','DNAJA1','DNAJA2','DNAJC19','DNAJC8')
pdf('Fig2-A.pdf',width = 6, height = 5)
DotPlot(cd4.combined, features = genelist,dot.scale = 8, cols =c('gray','red'),col.min = 0) + RotatedAxis() +
  theme(axis.title=element_blank(),legend.text = element_text(size = 12),legend.title=element_text(size=12),
        axis.text = element_text(size=18),legend.position = 'right')
dev.off()


# Expression level of HSPA6 and HSPA2 on UMAP projection of CD4+ T cells
pdf('Fig2-B.pdf',width = 6,height = 3)
FeaturePlot(cd4.combined,c('HSPA6','HSPA2'),cols = c('lightgray','red'),ncol = 2, order = T)
dev.off()


# Cell composition of HSPhi T based on cell sources.
df <- as.data.frame(table(hsp$sample))
names(df) <- c('category','count')
df$fraction<-df$count/sum(df$count)
df$ymax<-cumsum(df$fraction)
df$ymin<-c(0,head(df$ymax,n=-1))
df$labelPosition<-(df$ymax + df$ymin)/2
df$label<-paste0(df$category,"\n value: ",df$count)

pdf('Fig2-C.pdf',width = 3.5,height = 2.5)
ggplot(df,aes(ymax=ymax,ymin=ymin,
              xmax=4,xmin=3))+
  geom_rect(aes(fill=category))+
  # geom_label(x=3.5,aes(y=labelPosition,label=label),size=4)+
  scale_fill_manual(values = c('#7e76c9','#a7cee2','#e21f26','#f6999a')) +
  coord_polar(theta = "y")+
  xlim(2,4)+
  theme_void()+
  theme(legend.position = "right",text = element_text(size=16)) +
  ggtitle('HSPhi') 
dev.off()


# Volcano plot of DEGs between HSPhi T and the other T cells (p < 0.01). Red points represent HSPhi T cell specific high expressed genes, while blue points represent the other T cell specific highly expressed genes
Dat <- FindMarkers(object = cd4.combined , ident.1 = 'HSPhi',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)  #将行名作为第一列‘genes’

Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 1.5, ifelse(Dat$avg_log2FC>= 1.5 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Dat[Dat$threshold == 'Up',] %>% write.csv(file = 'gene_hsphi.csv')

pdf('Fig2-D.pdf',width = 5,height = 3)
g<-ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#point color
  theme_bw()+#Modify background
  ylab('-log10 (p-adj)')+#y-axis name
  xlab('log2 (FoldChange)')+#x-axis name
  geom_vline(xintercept=c(-1.5,1.5),lty=3,col="black",lwd=0.5) +#Add horizontal line|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+  #Add vertical line p_val_adj<0.05
  ggtitle(label = 'HSPhi' )  ;
g <- LabelPoints(plot = g,points = Dat[Dat$threshold != 'NoSignifi','allele'], repel = T) ; g
dev.off()


# GO enrichment analysis of HSPhi T specific expressed genes
hsp_go.data <- data.frame('GO_terms' = HSP_GO[1:9,1],
                          'p_value' = HSP_GO[1:9,2])
hsp_go.data$term <- factor(hsp_go.data$term,levels = rev(hsp_go.data$term))

pdf('Fig2-E.pdf',width = 11, height = 5)
ggplot(data=hsp_go.data, aes(x=term, y=..log10p)) +
  geom_bar(stat="identity", width=0.2, fill = 'black') + coord_flip() + 
  scale_fill_manual(values = 'black') + theme_test()  +
  xlab(label = element_blank()) + ylab("Log10(p-value)") +
  theme(axis.text=element_text(color="black", size = 15),panel.border = element_blank()) +　
  ylim(0,-15)
dev.off()


# UMAP projection of CD4 T cells using gene sets that removed HSP related genes, in which HSPhi T cell cluster disappears
Idents(hsp.combined) <- factor(Idents(hsp.combined), levels = levels(cd4.combined@active.ident))

pdf('Fig2-F.pdf',width = 5.5,height = 5.5)
DimPlot(hsp.combined,reduction = 'umap',cols = set3, label = F)  +
  theme(legend.position='none',legend.text = element_blank()) 　+ labs(title  = element_blank()) + scale_y_reverse()
dev.off()
table(hsp@meta.data[["predicted.celltype"]])



# GO enrichment analysis of HSPhi T specific expressed genes
# hsp.combined@active.ident <- factor(hsp.combined@active.ident, levels = c("TN","TCM","TEM","TEMRA","Treg","cytohi TEFF","Conv TEFF","Prolif TEFF",'HSPhi'))
pdf('Fig2-E.pdf',width = 4,height = 4)
DimPlot(hsp.combined,group.by = 'hsp',reduction = 'umap',cols = c('red','lightgray'),pt.size = 0.1)  +
  theme(legend.position=c(0.7,0.9)) 　+ labs(title  = element_blank()) + scale_y_reverse()
dev.off()
#Merged result (hsp becomes other type)


# recluster no_hsp_gene pie chart (proportion of hsp cells redistributed in other subpopulations)

df<-data.frame(celltype=c('TN','TCM','Treg','cytohi TEFF','CTLA4hi TEFF','Conv TEFF'),
               count=c(164,6,7,1,14,36))
df$celltype <- factor(df$celltype,levels = c('TN','TCM','Treg','cytohi TEFF','CTLA4hi TEFF','Conv TEFF'))

library(ggpubr)
pdf('Fig2-H.pdf',width = 4,height = 4)
ggpie(data = df,x = "count",label = rep("",6), fill = 'celltype',color = 'white', palette = c(set3[1],set3[2],set3[3],set3[4],set3[6],set3[7])) + theme(legend.position = "right",legend.text = element_text(size = 15),legend.title = element_text(size = 15)) 
dev.off()

c(164,6,7,1,14,36)/sum(c(164,6,7,1,14,36))


# Re-dimensionality reduction with HSPhi T-specific genes
Resting_1 <- subset(cd4.combined, subset = orig.ident == 'Resting_1')
Resting_2 <- subset(cd4.combined, subset = orig.ident == 'Resting_2')
Activated_DB1 <- subset(cd4.combined, subset = orig.ident == 'Activated_DB1')
Activated_DB2 <- subset(cd4.combined, subset = orig.ident == 'Activated_DB2')

try.anchors <- FindIntegrationAnchors(object.list = c(Resting_1,Resting_2,Activated_DB1,Activated_DB2)
                                      , dims = 1:20,anchor.features = HSP.gene)
try.combined <- IntegrateData(anchorset = try.anchors, dims = 1:20,features.to.integrate = HSP.gene)
DefaultAssay(try.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
try.combined <- ScaleData(try.combined, verbose = FALSE,features = HSP.gene)#
try.combined <- RunPCA(try.combined, npcs = 30, verbose = FALSE,features = HSP.gene) #
# t-SNE and Clustering
try.combined <- RunUMAP(try.combined, reduction = "pca",dims = c(1:20))  # 

rownames(cd4.combined@meta.data[which(cd4.combined$celltype == 'HSPhi'),]) %in%
  rownames(try.combined@meta.data)
rownames(try.combined@meta.data)%in%rownames(cd4.combined@meta.data[which(cd4.combined$celltype == 'HSPhi'),]) 
try.combined$hsp <- 'other'
try.combined@meta.data[rownames(try.combined@meta.data)%in%rownames(cd4.combined@meta.data[which(cd4.combined$celltype == 'HSPhi'),]) ,'hsp'] <- 'HSPhi'

try.combined$hsp <- factor(try.combined$hsp, levels = c('HSPhi','other'))
pdf('Fig2_I.pdf',width = 4, height = 4)
DimPlot(try.combined,group.by = 'hsp',cols = c('red','lightgray')) +
  theme(legend.position=c(0.7,0.95)) + labs(title  = element_blank())
dev.off()
