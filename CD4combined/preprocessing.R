# ------------- CD4combined_scRNAseq_data_preprocessing -------------
# This code shows how data were integrated and filtered.

# LOADING R LIBRARIES

library(ggplot2)
library(RColorBrewer) 
library(ggrepel)
library(rgl)
library(scales)
library(Seurat)
library(cowplot)
library(reticulate)
library(magrittr) #%>%
library(dplyr)  #group_by()
set3 = brewer.pal(n = 12,name = 'Set3')
set3 = c(set3,"#A6CEE3")
paired = brewer.pal(n = 12,name = 'Paired')
accent = brewer.pal(n = 8,name = 'Accent')
set.seed(12345)


## LOADING 10X DATA
Resting_1.data <-  Read10X(data.dir = "./cd4r_filtered_feature_bc_matrix/")
Activated_DB1.data <- Read10X(data.dir = "./cd4a_filtered_feature_bc_matrix/")

# CREATING SEURAT OBJECTS
Resting_1 <- CreateSeuratObject(counts = Resting_1.data, project = "Resting_1", min.cells = 5)
Activated_DB1 <- CreateSeuratObject(counts = Activated_DB1.data, project = "Activated_DB1", min.cells = 5)
Resting_2 <- readRDS("./GSM4450387_unstimulated_full_seurat.rds")
Activated_DB2 <- readRDS("./GSM4450386_stimulated_full_seurat.rds")

# UPDATEING SEURAT OBJECT
Resting_2<-UpdateSeuratObject(Resting_2)
Activated_DB2<-UpdateSeuratObject(Activated_DB2)

# sample
Resting_1$sample <- "Resting_1"
Activated_DB1$sample <- "Stimulated_1"
Resting_2$sample <- "Resting_2"
Activated_DB2$sample <- "Activated_DB2"

# activated
Resting_1$activated <- "rest"
Activated_DB1$activated <- "stim"
Resting_2$activated <- "rest"
Activated_DB2$activated <- "stim"

## Integrating data
DefaultAssay(cd4.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
cd4.combined <- ScaleData(cd4.combined, verbose = FALSE)
cd4.combined <- FindVariableFeatures(cd4.combined, selection.method = "vst", nfeatures = 2000, assay = 'RNA')
cd4.combined <- RunPCA(cd4.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
cd4.combined <- FindNeighbors(cd4.combined, reduction = "pca", dims = 1:20)
cd4.combined <- FindClusters(cd4.combined, resolution = 1.4)
DefaultAssay(cd4.combined) <- 'RNA'

## FILTERING OUT LOW QUALITY CELLS CELLS
# Red Blood Cell
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(cd4.combined@assays$RNA)) 
HB.genes <- rownames(cd4.combined@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
cd4.combined[["percent.HB"]]<-PercentageFeatureSet(cd4.combined, features=HB.genes) 
# Mitochondrial content
cd4.combined[["percent.mt"]] <- PercentageFeatureSet(cd4.combined, pattern = "^MT-")
# CD8 T cell
CD8.genes <- c("CD8A","CD8B")  
CD8_m <- match(CD8.genes, rownames(cd4.combined@assays$RNA)) 
CD8.genes <- rownames(cd4.combined@assays$RNA)[CD8_m]  
CD8.genes <- CD8.genes[!is.na(CD8.genes)] 
cd4.combined[["percent.CD8"]]<-PercentageFeatureSet(cd4.combined, features=CD8.genes) 

cd4.combined <- subset(cd4.combined,subset = nFeature_RNA > 200 & nCount_RNA < 60000 
                       & percent.CD8 < 0.01 & percent.HB<0.1 & percent.mt <7.5)

saveRDS(object = cd4.combined,file = "CD4combined.rds")
