# ------------- CD4combined_scRNAseq_inertT -------------
# This code shows the expression difference of inert T cells to activated or resting
# T cells. The cells were first classfied into three categories and some specific 
# gene expression were analyzed. The expressions of CXCR4 were also focused.

# classfying cells 
cd4.combined@meta.data[["inert"]] <- factor(cd4.combined@meta.data[["inert"]],levels = c('resting','inert','activated'))
cd4.combined$inert <- 'resting'
levels(cd4.combined$celltype)
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'TN'),'inert'] <- 'inert'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'TCM'),'inert'] <- 'inert'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'TEM'),'inert'] <- 'inert'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'ISAGhi'),'inert'] <- 'inert'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'TEMRA'),'inert'] <- 'inert'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'Treg'),'inert'] <- 'inert'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'cytohi TEFF'),'inert'] <- 'activated'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'CTLA4hi TEFF'),'inert'] <- 'activated'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'Conv TEFF'),'inert'] <- 'activated'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'Prolif TEFF'),'inert'] <- 'activated'
cd4.combined@meta.data[which(cd4.combined$activated == 'stim' & cd4.combined$celltype == 'HSPhi'),'inert'] <- 'activated'
table(cd4.combined@meta.data[["inert"]])

cd4.combined$inert <- factor(cd4.combined$inert, levels = c("resting",'inert','activated'))
CD4.combined <- subset(cd4.combined, subset = inert != 'activated')
CD4.combined$inert <- factor(CD4.combined$inert, levels = c("resting","inert"))

pdf('Fig3-A.pdf',width = 4,height = 4)
DimPlot(cd4.combined,group.by = 'inert',cols = rev(c(set3[4],set3[5],'lightgray'))) + theme(legend.position=c(0.7,0.8)) 
dev.off()

## Ternary diagram of specific genes expression in three categories 
library('ggtern')

# Calculate specifically expressed genes
Idents(cd4.combined) <- 'inert'
Dat <- FindMarkers(object = cd4.combined , ident.1 = 'inert', ident.2 = 'resting',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)  
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 1, ifelse(Dat$avg_log2FC>= 1 ,'inert','resting'),'NoSignifi'),levels=c('inert','resting','NoSignifi'))
Dat %>% write.csv(file = 'gene_i_r.csv')

Dat <- FindMarkers(object = cd4.combined , ident.1 = 'activated', ident.2 = 'inert',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)  
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 1, ifelse(Dat$avg_log2FC>= 1 ,'activated','inert'),'NoSignifi'),levels=c('activated','inert','NoSignifi'))
Dat　%>% write.csv(file = 'gene_a_i.csv')

Dat <- FindMarkers(object = cd4.combined , ident.1 = 'activated', ident.2 = 'resting',assay = 'RNA',logfc.threshold = 0)
Dat=cbind(allele=row.names(Dat), Dat)  
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 1, ifelse(Dat$avg_log2FC>= 1 ,'activated','resting'),'NoSignifi'),levels=c('activated','resting','NoSignifi'))
Dat　%>% write.csv(file = 'gene_a_r.csv')

sig.gene <- read.csv(file = 'CD4_Th/gene_a_i.csv',header = T)
sig.gene <- sig.gene[sig.gene$threshold != 'NpSignifi',]
gene.ac1 <- sig.gene[sig.gene$threshold == 'activated',]$allele
gene.in1 <- sig.gene[sig.gene$threshold == 'inert',]$allele

sig.gene <- read.csv(file = 'CD4_Th/gene_i_r.csv',header = T)
sig.gene <- sig.gene[sig.gene$threshold != 'NpSignifi',]
gene.re1 <- sig.gene[sig.gene$threshold == 'resting',]$allele
gene.in2 <- sig.gene[sig.gene$threshold == 'inert',]$allele

sig.gene <- read.csv(file = 'CD4_Th/gene_a_r.csv',header = T)
sig.gene <- sig.gene[sig.gene$threshold != 'NpSignifi',]
gene.ac2 <- sig.gene[sig.gene$threshold == 'activated',]$allele
gene.re2 <- sig.gene[sig.gene$threshold == 'resting',]$allele

gene.in <- unique(c(gene.in1,gene.in2))
gene.re <- unique(c(gene.re1,gene.re2))
gene.ac <- unique(c(gene.ac1,gene.ac2)) 
rm(gene.in1,gene.in2,gene.re1,gene.re2,gene.ac1,gene.ac2)

# combined all gene
sig.gene <- unique(c(gene.ac,gene.re,gene.in))
cluster.averages <- AverageExpression(cd4.combined,assays = 'RNA',group.by = 'inert')
cluster.averages <- as.data.frame(cluster.averages$RNA)
cluster.averages$gene <- rownames(cluster.averages)
sig.Dat <-  cluster.averages[cluster.averages$gene %in% sig.gene,]

# classfying genes, based on the max value of gene expression of three conditions
for (i in c(1:nrow(sig.Dat))) {
  maxvalue <- max(sig.Dat[i,1:3]) 
  sig.Dat[i,'group'] <- colnames(sig.Dat)[which(sig.Dat[i,] == maxvalue)]
}

# Ternary diagram of specific genes
pdf('Fig3-B.pdf')
ggtern(data=sig.Dat,aes(resting, inert, activated))+     
  geom_point(size=2.5,aes(color=group),alpha=0.8)+
  scale_colour_manual(values = c(set3[4],set3[5],'lightgray'))+ 
  theme_rgbw(base_size = 17)+   
  theme(plot.title = element_text(size=15,hjust = 0.5),legend.position = c(0.9,0.8),legend.text = element_text(size=13))
dev.off()

## Venn diagram showing the correlation of high expressing genes in five subsets
Idents(cd4.combined) <- cd4.combined$celltype 
for (i in list('TN','TCM','TEM','TEMRA','Treg')) {
  mydata = subset(cd4.combined, idents = i)
  Dat <- FindMarkers(object = mydata , ident.1 = 'inert', ident.2 = 'resting',group.by = 'inert',assay = 'RNA',logfc.threshold = 1,fc.name = 'avg_log2FC',only.pos = T)
  write.csv(Dat,paste0("vol_",i,"_res.csv"))
}
rm(mydata,Dat,i)

venn.diagram(
  x = list(na.omit(inert_1$TN), na.omit(inert_1$TCM), na.omit(inert_1$TEM), na.omit(inert_1$TEMRA), na.omit(inert_1$Treg)),
  category.names = c('TN', 'TCM', 'TEM', 'TEMRA', 'Treg'),
  fill = c(paired[1],paired[7],paired[8],paired[3],paired[9]),
  alpha = 0.8,
  cex = 1.5,
  filename = 'Fig3-C.png',
  cat.col = c(paired[1],paired[7],paired[8],paired[3],paired[9])
)

## Heatmap of differential gene of inert T and activated T cell in various subsets
#  Calculate the average expression of genes
Idents(CD4.combined) <- 'celltype'
cluster.averages <- AverageExpression(CD4.combined,assays = 'RNA',add.ident = 'inert')

# normalizing expression
head(cluster.averages[['RNA']])
cluster.averages$RNA<-log(cluster.averages$RNA)
is.na(cluster.averages$RNA) <- sapply(cluster.averages$RNA, is.infinite) 

# gene selecting based on biological knowledge
markers.match<- match(c("CTLA4","MIR155HG","IL2RA","IL2",               # 激活
                        "HSPD1","HSP90AB1","HSPE1","DNAJA1","HSP90AA1",               # HSP 热休克蛋白及相关蛋白
                        "CCL4L2","XCL1","XCL2","CCL3L3","CCL4","CCL3","CCL20","IFNG","LTA","MIF","ARID5A",    # 细胞因子
                        "ISG15","STAT1","IRF4","SOCS1","ZBED2",     # IFN干扰素响应
                        "FABP5","CREM","NME1","DDX21","TUBA1B","CCND2","CYTOR","DDIT4",   # 增殖分裂、细胞周期、抑制凋亡相关
                        "HILPDA","MTHFD2","PGK1","SLC3A2","GBP2",      # 增殖分裂、细胞周期、抑制凋亡相关
                        "TNF","TNFRSF4","TNFRSF9","TNFRSF18",  # TNF家族
                        "ENO1","EIF5A","NCL","BATF",#转录翻译
                        "LDHA","ALDOA","PKM","TPI1","GAPDH",# 糖酵解
                        "PIM3","SERBP1","NAMPT","BIRC3","LGALS1","GZMB","PRDX1",  # 凋亡、衰老
                        "RPS26","SRM","TXN","NFKBIA","AC017002.1","TYMP","FTL"),rownames(cluster.averages$RNA))  #知道指定marker的行数位置

# genes filtering 
markers.to.plot <- rownames(cluster.averages$RNA)[markers.match]  
new.average.expression <- cluster.averages$RNA[markers.match,]  
new.average.expression[is.na(new.average.expression)] <- 0  
new.average.expression <- as.data.frame(t(new.average.expression))
new.average.expression <- new.average.expression[c("Treg_resting","TEMRA_resting","TEM_resting","TCM_resting","TN_resting"   
                                                   ,"Treg_inert","TEMRA_inert","TEM_inert","TCM_inert","TN_inert"),]  # df 按行名重新排序
new.average.expression <- as.data.frame(t(new.average.expression))

bk = unique(c(seq(-1.5,1.5, length=100)))    # Set the color range

# Filter the results of the clusters and rearrange the dataframe to remove the 
# effect of the cluster tree drawing
p1 <- pheatmap::pheatmap(new.average.expression,breaks = bk,scale = 'row',
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), border=FALSE,
                         cluster_row =T,cluster_col =F)  ; p1
# Reorder by selected genes and sample numbers
gn <- rownames(new.average.expression)[p1$tree_row[["order"]]] ; gn
new.average.expression <- new.average.expression[gn,]
labels_col = c("Treg","TEMRA","TEM" ,"TCM" ,"TN","Treg","TEMRA","TEM" ,"TCM" ,"TN")
annotation_col <- data.frame(type = rep(c("activated", "inert"), each = 5))
rownames(annotation_col) <- colnames(new.average.expression)[1:10]

# drew the heatmap 
p1 <- pheatmap::pheatmap(new.average.expression,breaks = bk,scale = 'row',
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), border=FALSE,
                         cluster_row =F,cluster_col =F,  labels_col = labels_col,  
                         annotation_col = annotation_col, annotation_colors = list(type=c(inert = set3[5],activated = set3[4])),
                         angle_col = "45",gaps_col = c(5,10))

# add flag genes 
source("add_flag.R")
library(grid)
gene_name<-c("CTLA4","MIR155HG","IL2RA",              
             "HSPD1","HSPE1",               # HSP heat shock proteins and related proteins
             "CCL4","CCL3","IFNG",    #Cytokine
             "ISG15","IRF4",     # IFN-interferon response
             "FABP5","CREM","NME1","TUBA1B",  # Proliferation and division, cell cycle, inhibition of apoptosis
             "TNF",  # TNF family
             "ENO1","EIF5A",#Transcription and translation
             "LDHA","ALDOA",
             "PIM3","LGALS1")
add.flag(p1,kept.labels = gene_name,repel.degree = 0.2) %>% ggsave(filename = 'Fig3-D.png',width = 6,height = 6)

# vlnplot of specific genes in three categories
pdf('Fig3-E.pdf',width = 7.5,height = 5)
VlnPlot(cd4.combined, features = c("MIR155HG","IL2","FABP5","IER3",'CD52','CXCR4','CORO1A','TXNIP'),group.by = 'inert',rev(c(set3[4],set3[5],'lightgray')),pt.size = 0,ncol = 4) +
  theme(legend.position='none', axis.title = element_blank())
dev.off()

# CXCR4 expression of stimulated or resting T cells 
# CXCR4 Featureplot
pdf('Fig3-F1.pdf',width = 6,height =3)
FeaturePlot(cd4.combined, features = c("CXCR4"), split.by = "activated",cols = c('#add8e6','#f08080'))
dev.off()

# Vlnplot
test.data <- subset(cd4.combined, subset = activated=='stim' | celltype %in%  c('TN',  'TCM',  'TEM', 'TEMRA', 'Treg')) 
pdf('Fig3-F2.pdf',width = 8,height =3)
VlnPlot(test.data, features = 'CXCR4',c('#add8e6','#f08080'),pt.size = 0,split.by = 'activated') +
  theme(legend.position='right', axis.title = element_blank())
dev.off()
