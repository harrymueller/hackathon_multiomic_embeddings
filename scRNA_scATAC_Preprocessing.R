
library(Signac)
library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(SeuratData)
library(dplyr)
library(patchwork)
library(rliger)


RNA <- readRDS("D:/Hackathon/Annotated_Processed_unpaired_RNA.RDS")
ATAC <- readRDS("D:/Hackathon/Annotated_processed_unpaired_ATAC.RDS")


load("D:/Hackathon/Annotated_object_Monaco_TimP.RData")


RNA <- readRDS("D:/Hackathon/pbmc_rna_unpaired.RDS")

RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = "^MT-")

VlnPlot(RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(RNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(RNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2

RNA <- subset(RNA, subset = nFeature_RNA < 4500 & percent.mt < 22 & nCount_RNA < 14000)

RNA <- FindVariableFeatures(RNA, nfeatures = 2000)
RNA <- NormalizeData(RNA)
RNA <- ScaleData(RNA)
RNA <- RunPCA(RNA, npcs = 20)
ElbowPlot(RNA)

RNA <- RunUMAP(RNA, dims = 1:20)
RNA <- FindNeighbors(RNA, dims = 1:20)
RNA <- FindClusters(RNA, resolution = 0.3)

DimPlot(RNA,group.by  = "Celltype", reduction = "umap")

FeaturePlot(RNA, features = "BANK1")

common.cells <- intersect(colnames(RNA), colnames(sobj))

RNA <- AddMetaData(
  object   = RNA,
  metadata = sobj@meta.data[common.cells, "Celltype", drop = FALSE]
)


pbmc.markers <- FindAllMarkers(RNA, only.pos = TRUE)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(RNA, features = top10$gene) + NoLegend()


saveRDS(RNA,file='D:/Hackathon/Annotated_Processed_unpaired_RNA.RDS')





ATAC <- readRDS("D:/Hackathon/pbmc_atac_unpaired.RDS")

granges(ATAC)
peaks.keep <- seqnames(granges(ATAC)) %in% standardChromosomes(granges(ATAC))
ATAC <- ATAC[as.vector(peaks.keep), ]


library(AnnotationHub)
ah <- AnnotationHub()

query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]

annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(ATAC) <- annotations

frags <- ATAC@assays$ATAC@fragments
frags[[1]]@path <- "D:/Hackathon/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
ATAC@assays$ATAC@fragments <- frags


ATAC <- NucleosomeSignal(object = ATAC)

ATAC <- TSSEnrichment(object = ATAC)



metrics <- read.csv("D:/Hackathon/pbmc_granulocyte_sorted_3k_per_barcode_metrics.csv",
                    stringsAsFactors = FALSE)
rownames(metrics) <- metrics$barcode
meta_to_add <- metrics[, c("atac_fragments", "atac_peak_region_fragments")]

matched <- intersect(colnames(ATAC), rownames(meta_to_add))
length(matched)            
length(matched) == ncol(ATAC)


ATAC <- AddMetaData(
  ATAC,
  metadata = meta_to_add
)

ATAC$pct_reads_in_peaks <- 
  ATAC$atac_peak_region_fragments / ATAC$atac_fragments * 100



ATAC$blacklist_ratio <- FractionCountsInRegion(
  object  = ATAC,
  assay   = "ATAC",               
  regions = blacklist_hg38_unified
)

VlnPlot(
  ATAC,
  features = c("pct_reads_in_peaks", "blacklist_ratio",
               "TSS.enrichment", "nucleosome_signal","nCount_ATAC"),
  ncol     = 5
) + NoLegend()


DensityScatter(ATAC, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

ATAC <- subset(
  x = ATAC,
  subset = nCount_ATAC > 9000 &
    nCount_ATAC < 55000 &
    pct_reads_in_peaks > 60 &
    blacklist_ratio < 0.0025 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 4
)

ATAC <- RunTFIDF(ATAC)
ATAC <- FindTopFeatures(ATAC, min.cutoff = 'q0')
ATAC <- RunSVD(ATAC)

DepthCor(ATAC)

ATAC <- RunUMAP(object = ATAC, reduction = 'lsi', dims = 2:30)
ATAC <- FindNeighbors(object = ATAC, reduction = 'lsi', dims = 2:30)
ATAC <- FindClusters(object = ATAC, verbose = FALSE, algorithm = 3,resolution = 0.5)


common.cells <- intersect(colnames(ATAC), colnames(sobj))

ATAC <- AddMetaData(
  object   = ATAC,
  metadata = sobj@meta.data[common.cells, "Celltype", drop = FALSE]
)

DimPlot(object = ATAC, group.by = "Celltype",label = TRUE) + NoLegend()

gene.activities <- GeneActivity(ATAC)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
ATAC[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)


DefaultAssay(ATAC) <- "ACTIVITY"
ATAC <- NormalizeData(ATAC)
ATAC <- ScaleData(ATAC)




transfer.anchors <- FindTransferAnchors(
  reference = RNA,
  query = ATAC,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA$Celltype,
  weight.reduction = ATAC[['lsi']],
  dims = 2:30
)


ATAC <- AddMetaData(object = ATAC, metadata = predicted.labels)

plot1 <- DimPlot(
  object = RNA,
  group.by = 'Celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = ATAC,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

saveRDS(ATAC,file='D:/Hackathon/Annotated_processed_unpaired_ATAC.RDS')



p1 <- DimPlot(ATAC, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(ATAC, group.by = "Celltype", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2



genes.use <- VariableFeatures(RNA)
refdata <- GetAssayData(RNA, assay = "RNA", layer = "data")[genes.use, ]


imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = ATAC[["lsi"]],
                           dims = 2:30)

ATAC[["RNA"]] <- imputation
RNA$dataset <- "RNA"
ATAC$dataset <- "ATAC"

coembed <- merge(x = RNA, y = ATAC)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
genes.use <- VariableFeatures(RNA)
coembed <- RunCCA(RNA, ATAC, features = genes.use)

DimPlot(coembed, group.by = "dataset", reduction='umap')


integrated_cca_mat <- Embeddings(coembed, reduction = "cca")

integrated_pca_mat  <- Embeddings(coembed, "pca")

integrated_pca_umap_mat <- Embeddings(coembed, "umap")


lsi_mat  <- Embeddings(ATAC, "lsi")

umap_mat <- Embeddings(ATAC, "umap")

dim(lsi_mat)   
dim(umap_mat)
dim(cca_mat)

write.csv(lsi_mat,
          file = "D:/Hackathon/ATAC_lsi_embedding.csv",
          quote = FALSE,
          row.names = TRUE,    # 细胞条码
          col.names = NA)

write.csv(umap_mat,
          file = "D:/Hackathon/ATAC_umap_embedding.csv",
          quote = FALSE,
          row.names = TRUE,
          col.names = NA)

write.csv(integrated_cca_mat,
          file = "D:/Hackathon/integrated_cca_embedding.csv",
          quote = FALSE,
          row.names = TRUE,    # 细胞条码
          col.names = NA)

write.csv(integrated_pca_mat,
          file = "D:/Hackathon/integrated_pca_embedding.csv",
          quote = FALSE,
          row.names = TRUE,    # 细胞条码
          col.names = NA)

write.csv(integrated_pca_umap_mat,
          file = "D:/Hackathon/integrated_pca_umap_embedding.csv",
          quote = FALSE,
          row.names = TRUE,    # 细胞条码
          col.names = NA)

genes.use <- VariableFeatures(coembed)

rna_matrix <- coembed[["RNA"]]$data[genes.use, ]
rna_df <- as.data.frame(t(as.matrix(rna_matrix)))  # 转置为 cell × gene


activity_genes <- intersect(genes.use, rownames(coembed[["ACTIVITY"]]$data))

activity_matrix <- coembed[["ACTIVITY"]]$data[activity_genes, ]
activity_df <- as.data.frame(t(as.matrix(activity_matrix)))  


metadata_df <- coembed[[]]  
write.csv(metadata_df, "D:/Hackathon/integrated_metadata.csv")

common_genes <- intersect(colnames(rna_df), colnames(activity_df))
rna_common <- rna_df[, common_genes, drop = FALSE]
activity_common <- activity_df[, common_genes, drop = FALSE]
rna_common$source <- "RNA"
activity_common$source <- "ACTIVITY"
combined_df <- rbind(rna_common, activity_common)
write.csv(combined_df, "D:/Hackathon/integrated_cell_gene_matrix.csv")







library(rliger)  


bmmcLiger <- createLiger(list(atac = ATAC, rna = RNA),
                         modal = c("atac", "rna"))

bmmcLiger <- bmmcLiger %>%
  normalize() %>%
  selectGenes(useDatasets = "rna") %>%
  scaleNotCenter()

bmmcLiger <- runIntegration(bmmcLiger, k = 20)

bmmcLiger <- quantileNorm(bmmcLiger)

bmmcLiger <- runCluster(bmmcLiger, nNeighbors = 30, resolution = 0.2)

bmmcLiger <- runUMAP(bmmcLiger, nNeighbors = 30, minDist = 0.3)

options(ligerDotSize = 0.5)
plotByDatasetAndCluster(bmmcLiger)


nmf_embedding <- bmmcLiger@H.norm  
nmf_df <- as.data.frame(nmf_embedding)
write.csv(nmf_df, "D:/Hackathon/integrated_NMF_embedding.csv")


meta <- as.data.frame(bmmcLiger@cellMeta)
umap <- as.data.frame(bmmcLiger@dimReds$UMAP)
colnames(umap) <- c("UMAP_1", "UMAP_2")

meta_umap <- cbind(meta, umap)
write.csv(meta_umap, "D:/Hackathon/liger_metadata.csv")

rna_mat <- as.matrix(bmmcLiger@datasets$rna@normData)  # gene × cell
rna_mat <- t(rna_mat)  # 转置为 cell × gene
rna_df <- as.data.frame(rna_mat)

atac_mat <- as.matrix(bmmcLiger@datasets$atac@normData)
atac_mat <- t(atac_mat)
atac_df <- as.data.frame(atac_mat)

common_genes <- intersect(colnames(rna_df), colnames(atac_df))
rna_df_sub <- rna_df[, common_genes]
atac_df_sub <- atac_df[, common_genes]

expr_combined <- rbind(rna_df_sub, atac_df_sub)
write.csv(expr_combined, "D:/Hackathon/liger_cell_gene_expression_combined.csv")



