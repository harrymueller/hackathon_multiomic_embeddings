# integration of unpaired scRNA-seq and ATAC-seq datasets using cca dimensionality reduction #

# load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2) 
library(sctransform) 
library(dplyr)
library(harmony)
library(Signac)
library(cowplot)

# load data
## scRNA-seq dataset
rna <- readRDS("/Users/a1798837/repos/R_projects/hackathon/Annotated_Processed_unpaired_RNA.RDS")
rna <- UpdateSeuratObject(rna)


## scATAC-seq dataset 
atac <- readRDS("/Users/a1798837/repos/R_projects/hackathon/Annotated_processed_unpaired_ATAC.RDS")

# make metadata column names consistent 
## add modality column to distinguish between rna and atac for visualisation 
rna$modality <- "rna"
atac$modality <- "atac"
## add a clusters column
rna$clusters <- rna$seurat_clusters
atac$clusters <- atac$seurat_clusters
## add a column named cell_type 
rna$cell_type <- rna$Celltype
atac$cell_type <- atac$Celltype

head(rna@meta.data) # check these additions were successfully added into the rna meta data
head(atac@meta.data) # check these additions were successfully added into the atac meta data

# plot the original unintegrated data 
p1 <- DimPlot(rna, group.by = "Celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac, group.by = "cell_type", label = FALSE, repel = TRUE) + NoLegend() + ggtitle("ATAC")
p1 + p2 # plot umaps

# ----------------------------------------------------------------------------------------------------------------------- #
# integration of unpaired scRNA- and ATAC-seq datasets using rpca dim reduction #
# quantify gene activity
DefaultAssay(atac) <- "ACTIVITY"  # switch to the gene activity assay
gene.activities <- GetAssayData(atac, assay = "ACTIVITY", slot = "data")

# add gene activities as a new assay
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))
rna <- RunPCA(rna, features = VariableFeatures(rna), npcs = 30)

# identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "rpca")

# annotate scATAC-seq cells via label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$cell_type,
                                     weight.reduction = atac[["lsi"]], dims = 2:30)

atac <- AddMetaData(atac, metadata = celltype.predictions) # add cell type predicted labels to atac meta data

atac$annotation_correct <- atac$predicted.id == atac$cell_type # label agreement check

# plot the predicted annotation and org annotation to compared
p1 <- DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(atac, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2

# ----------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------- #
# assessment of integration - biological conservation metrics #

# load relevant libraries
library(Matrix)
library(cluster)
library(ape)  # for Moran's I and Geary's C
library(philentropy)  # for pairwise distances
library(FNN)  # for nearest neighbors

## marker-based metrics ##
# Define marker sets
t_cell_markers <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IL7R", "CCR7") # t cells
b_cell_markers <- c("MS4A1", "CD19", "CD79A") # b cells
nk_cell_markers <- c("NKG7", "GNLY", "KLRD1", "NCAM1", "FCGR3A") # NK cells 
monocyte_markers <- c("CD14", "LYZ", "FCN1", "FCGR3A") # monocytes
dendritic_cell_markers <- c("FCER1A", "CLEC9A", "BATF3", "IRF8", "CD1C") # dendritic cells
megakaryocyte_markers <- c("PPBP", "PF4", "ITGA2B", "GP1BA") # Megakaryocytes / Platelets

# score marker genes
atac <- AddModuleScore(object = atac, features = list(b_cell_markers), name = "BCell_Score", ctrl = 50, nbins = 25, 
                       seed = 0 ) # b cell markers
FeaturePlot(atac, features = "BCell_Score1", cols = c("lightgrey", "blue")) # plot gene scores
FeaturePlot(atac, features = b_cell_markers, ncol = 3) # plot gene markers 

# moran's I
## build KNN graph from UMAP
coords_atac <- Embeddings(atac, reduction = "umap")
knn_graph <- FNN::get.knn(coords_atac, k = 15)
w <- matrix(0, nrow = nrow(coords_atac), ncol = nrow(coords_atac))
for (i in 1:nrow(coords_atac)) {
  w[i, knn_graph$nn.index[i, ]] <- 1
}
diag(w) <- 0

# Moran’s I
moran <- Moran.I(atac$BCell_Score1, weight = w)
print(moran)


# custom marker
## create function to assess neighbours 
marker_spatial_scores <- function(object, gene_set, reduction = "umap", 
                                  n_neighbors = 15, n_distant = 15, n_random = 15, 
                                  agg = "mean", seed = 42) {
  coords <- Embeddings(object, reduction = reduction)
  
  # Average gene expression per cell
  expr_mat <- GetAssayData(object, slot = "data")[gene_set, , drop = FALSE]
  if (agg == "mean") {
    expr <- Matrix::colMeans(expr_mat)
  } else if (agg == "sum") {
    expr <- Matrix::colSums(expr_mat)
  } else {
    stop("agg must be 'mean' or 'sum'")
  }
  
  # Local neighbors
  nn <- FNN::get.knn(coords, k = n_neighbors)$nn.index
  local_means <- apply(nn, 1, function(idx) mean(expr[idx]))
  local_corr <- cor(expr, local_means, method = "spearman")
  
  # Distant neighbors
  dist_mat <- as.matrix(philentropy::distance(coords, method = "euclidean"))
  diag(dist_mat) <- -Inf
  far_idx <- t(apply(dist_mat, 1, function(x) order(x, decreasing = TRUE)[1:n_distant]))
  distant_means <- apply(far_idx, 1, function(idx) mean(expr[idx]))
  distant_corr <- cor(expr, distant_means, method = "spearman")
  
  # Random neighbors
  set.seed(seed)
  rand_idx <- replicate(length(expr), sample(seq_along(expr), n_random))
  rand_means <- apply(rand_idx, 2, function(idx) mean(expr[idx]))
  random_corr <- cor(expr, rand_means, method = "spearman")
  
  return(list(
    local_corr = local_corr, 
    distant_corr = distant_corr, 
    random_corr = random_corr
  ))
}

# run the metric
scores <- marker_spatial_scores(atac, gene_set = b_cell_markers)
print(scores)

# ----------------------------------------------------------------------------------------------------------------------- #
## annotation-based metrics ##

# load libraries 
library(Seurat)
library(cluster)         # For silhouette
library(FNN)             # For kNN
library(lisi)            # For LISI

# create the function to assess silhouette score and kNN purity
evaluate_embedding_metrics <- function(seurat_obj, 
                                       label_key,
                                       embedding_key = "integrated", 
                                       knn_k = 30) {
  # 1. Extract embedding and labels
  embedding <- Embeddings(seurat_obj, reduction = embedding_key)
  labels <- seurat_obj[[label_key]][, 1]
  
  # Check for valid input
  if (is.null(embedding)) stop("Embedding not found in Seurat object.")
  if (is.null(labels)) stop("Label column not found in metadata.")
  
  # 2. Silhouette score
  dist_mat <- dist(embedding)
  sil <- silhouette(as.integer(factor(labels)), dist_mat)
  silhouette_score <- mean(sil[, "sil_width"])
  
  # 3. kNN Purity
  knn <- get.knn(embedding, k = knn_k + 1)$nn.index[, -1]  # remove self
  knn_purity <- mean(sapply(1:nrow(knn), function(i) {
    mean(labels[knn[i, ]] == labels[i])
  }))
  
  return(list(
    silhouette = silhouette_score,
    knn_purity = knn_purity
  ))
}

# run the function on the dataset
results_atac <- evaluate_embedding_metrics(
  seurat_obj = atac,
  label_key = "predicted.id",
  embedding_key = "umap",
  knn_k = 15
)
print(results_atac)















