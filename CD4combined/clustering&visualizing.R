# ------------- CD4combined_scRNAseq_clustering_visualizing -------------
# This code shows how cells were clustered and annotated.

cd4.combined@meta.data$activated <- factor(cd4.combined@meta.data$activated, levels = c('rest','stim'))

# UMAP projection of CD4+ T cells pre-and post-anti-CD3/CD28 stimulation, colored by whether it was stimulated
pdf('Fig1-B.pdf',width = 4, height = 4)
DimPlot(cd4.combined, reduction = "umap", group.by = "activated",label = F,cols = c('#add8e6','#f08080'),pt.size = 0.1) +
  theme(legend.position=c(0.6,0.8))　+ labs(title  = element_blank())
dev.off()

# IL2RA KLF2 expression 
cd4.combined@meta.data$sample <- factor(cd4.combined@meta.data$sample, levels = c('Resting_1','Resting_2','Stimulated_1','Stimulated_2'))
pdf('Fig1-C.pdf',width = 3.5, height = 5)
VlnPlot(cd4.combined, features = c('IL2RA','KLF2'),group.by = 'sample',rev(c('#f6999a','#e21f26','#a7cee2','#7e76c9')),pt.size = 0,ncol = 1) +
  theme(legend.position='none', axis.title.x = element_blank())
dev.off()

# rename cluster idents
a <- c('TN','TEM','TN','TCM','TCM','TCM','CTLA4hi TEFF','Treg','Conv TEFF','TEM','Conv TEFF','cytohi TEFF','TEMRA','Prolif TEFF','TEMRA','Treg','TEMRA','TN','TEMRA','HSPhi','ISAGhi','Prolif TEFF')
#cd4.combined <- RenameIdents(cd4.combined, 'Th0' = 'IFNhi')
names(a) <- levels(cd4.combined)
cd4.combined <- RenameIdents(cd4.combined, a)
table(cd4.combined@active.ident)
cd4.combined@active.ident <- factor(cd4.combined@active.ident, levels = c('TN','TCM','TEM','ISAGhi','TEMRA','Treg','cytohi TEFF','CTLA4hi TEFF','Conv TEFF','Prolif TEFF','HSPhi'))
DimPlot(cd4.combined, reduction = "umap", label = T,cols = set3)
rm(a)
cd4.combined@meta.data[["celltype"]] <- cd4.combined@active.ident

# significant gene expression
markers.to.plot <-c('HSPA6','HSPA1A','HSPB1'  # HSPhi
                    ,'TUBB',"TUBA1B","FABP5","MIR155HG"  # proliferation T
                    ,'IL2','IFNG'# Th1
                    ,"IFNG","IL2","CCL20"  # cytokine T
                    ,"IL2RA","CTLA4","FOXP3"  #Treg
                    ,"GZMA","GZMK",'CST7'  #TEMRA
                    ,"IL7R" # TM
                    ,'TIMP1','S100A11',"ANXA1","ISG15" # TEM/TCM
                    ,"SELL","CCR7") #naive T
pdf('Fig1-D.pdf',width = 11, height = 5)
DotPlot(cd4.combined, features = rev(unique(markers.to.plot)),dot.scale = 8, cols =c('gray','red'),col.min = 0) + RotatedAxis() +
  theme(axis.title=element_blank(),legend.text = element_text(size = 12),legend.title=element_text(size=12),
        axis.text = element_text(size=18),legend.position = 'right')
dev.off()

# UMAP projection of CD4+ T cells pre-and post-anti-CD3/CD28 stimulation, colored by cell subsets
pdf('Fig1-E.pdf',width = 4, height = 4)
DimPlot(cd4.combined, reduction = "umap", label = F,cols = set3) + theme(legend.position='none')
dev.off()
