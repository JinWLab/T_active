# ------------- CD4combined_scRNAseq_lineage -------------
# This code shows how cell lineage were inferred.

# loading packages
library(monocle3)
library(ggplot2)
library(patchwork)
library(dplyr)  

# creating CDS obeject and pre-processing data 
sdata <- cd4.combined
data <- GetAssayData(sdata, assay = 'RNA', slot = 'counts') 
cell_metadata <- sdata@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))  
rownames(gene_annotation) <- rownames(data) 
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

# loading integrated UMAP from Seurat object
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sdata, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

# Learn the trajectory 
cds <- learn_graph(cds)
# Order the cells in pseudotime
cds <- order_cells(cds)   # Manually selecting root nodes
p <- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5)
ggsave('mono3_ordercells.png',p,width = 6,height = 5)
