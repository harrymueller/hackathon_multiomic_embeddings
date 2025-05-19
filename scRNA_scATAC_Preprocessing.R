library(Signac) # Load Signac for single-cell chromatin analysis
library(Seurat) # Load Seurat for single-cell data handling
library(ggplot2) # Load ggplot2 for plotting
library(SeuratDisk) # Load SeuratDisk for HDF5 conversions
library(SeuratData) # Load SeuratData for example datasets
library(dplyr) # Load dplyr for data manipulation
library(patchwork) # Load patchwork for combining ggplots
library(rliger) # Load rliger for joint matrix factorization

# Load pre-annotated Seurat object containing RNA data and metadata
load("D:/Hackathon/Annotated_object_Monaco_TimP.RData")

# Calculate percentage of mitochondrial gene expression per cell
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = "^MT-")

# Visualize QC metrics: number of features, counts, and mitochondrial percentage
VlnPlot(RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots for QC: counts vs mitochondrial percentage and counts vs feature count
plot1 <- FeatureScatter(RNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(RNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Combine the QC scatter plots
plot1 + plot2

# Filter cells based on QC thresholds
RNA <- subset(RNA, subset = nFeature_RNA < 4500 & percent.mt < 22 & nCount_RNA < 14000)

# Identify highly variable features for downstream analysis
RNA <- FindVariableFeatures(RNA, nfeatures = 2000)
# Normalize the data
RNA <- NormalizeData(RNA)
# Scale the data to remove unwanted sources of variation
RNA <- ScaleData(RNA)
# Perform principal component analysis
RNA <- RunPCA(RNA, npcs = 20)
# Plot elbow plot to decide number of PCs
ElbowPlot(RNA)

# Run UMAP for visualization using top PCs
RNA <- RunUMAP(RNA, dims = 1:15)
# Construct nearest neighbors graph
RNA <- FindNeighbors(RNA, dims = 1:15)
# Cluster the cells at specified resolution
RNA <- FindClusters(RNA, resolution = 0.3)

# Plot UMAP colored by existing Celltype annotation
DimPlot(RNA, group.by = "Celltype", reduction = "umap")

# Plot expression of BANK1 on UMAP
FeaturePlot(RNA, features = "BANK1")

# Identify cells common to both RNA and ATAC datasets
common.cells <- intersect(colnames(RNA), colnames(sobj))

# Transfer Celltype annotation from 'sobj' to RNA object
RNA <- AddMetaData(
  object = RNA,
  metadata = sobj@meta.data[common.cells, "Celltype", drop = FALSE]
)

# Find marker genes for each cluster
pbmc.markers <- FindAllMarkers(RNA, only.pos = TRUE)

# Select top 10 markers per cluster based on log fold change
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
# Plot heatmap of top markers
DoHeatmap(RNA, features = top10$gene) + NoLegend()

# Save the processed RNA object
saveRDS(RNA, file = 'D:/Hackathon/Annotated_Processed_unpaired_RNA.RDS')


# Load unpaired ATAC Seurat object
ATAC <- readRDS("D:/Hackathon/pbmc_atac_unpaired.RDS")

# Keep only peaks on standard chromosomes
peaks.keep <- seqnames(granges(ATAC)) %in% standardChromosomes(granges(ATAC))
ATAC <- ATAC[as.vector(peaks.keep), ]

# Load genomic annotations from AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()
# Query for EnsDb.Hsapiens.v98
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]
# Convert EnsDb to GRanges
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
# Set UCSC style seqlevels and genome
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
# Add annotation to ATAC assay
Annotation(ATAC) <- annotations

# Correct fragment file path for ATAC object
frags <- ATAC@assays$ATAC@fragments
frags[[1]]@path <- "D:/Hackathon/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
ATAC@assays$ATAC@fragments <- frags

# Compute nucleosome signal per cell
ATAC <- NucleosomeSignal(object = ATAC)
# Compute TSS enrichment per cell
ATAC <- TSSEnrichment(object = ATAC)

# Read per-barcode metrics for FRiP calculation
metrics <- read.csv(
  "D:/Hackathon/pbmc_granulocyte_sorted_3k_per_barcode_metrics.csv",
  stringsAsFactors = FALSE
)
# Set rownames to barcode
rownames(metrics) <- metrics$barcode
# Subset desired metrics
meta_to_add <- metrics[, c("atac_fragments", "atac_peak_region_fragments")]

# Add metrics to ATAC metadata
ATAC <- AddMetaData(ATAC, metadata = meta_to_add)
# Calculate fraction of reads in peaks (FRiP)
ATAC$pct_reads_in_peaks <- ATAC$atac_peak_region_fragments / ATAC$atac_fragments * 100

# Compute blacklist ratio for quality control
ATAC$blacklist_ratio <- FractionCountsInRegion(
  object  = ATAC,
  assay   = "ATAC",
  regions = blacklist_hg38_unified
)

# Visualize QC metrics for ATAC
VlnPlot(
  ATAC,
  features = c("pct_reads_in_peaks", "blacklist_ratio",
               "TSS.enrichment", "nucleosome_signal", "nCount_ATAC"),
  ncol     = 5
) + NoLegend()

# Density scatter of counts vs TSS enrichment
DensityScatter(
  ATAC,
  x = 'nCount_ATAC',
  y = 'TSS.enrichment',
  log_x = TRUE,
  quantiles = TRUE
)

# Subset ATAC cells based on QC thresholds
ATAC <- subset(
  x = ATAC,
  subset = nCount_ATAC > 9000 &
           nCount_ATAC < 55000 &
           pct_reads_in_peaks > 60 &
           blacklist_ratio < 0.0025 &
           nucleosome_signal < 1.5 &
           TSS.enrichment > 4
)

# Normalize ATAC data with TF-IDF
ATAC <- RunTFIDF(ATAC)
# Identify top features
ATAC <- FindTopFeatures(ATAC, min.cutoff = 'q0')
# Perform SVD (LSI)
ATAC <- RunSVD(ATAC)

# Compute correlation between depth and LSI components
DepthCor(ATAC)

# Run UMAP, neighbors, and clustering on LSI
ATAC <- RunUMAP(object = ATAC, reduction = 'lsi', dims = 2:30)
ATAC <- FindNeighbors(object = ATAC, reduction = 'lsi', dims = 2:30)
ATAC <- FindClusters(object = ATAC, verbose = FALSE, algorithm = 3, resolution = 0.5)

# Transfer cell type annotation from sobj to ATAC
common.cells <- intersect(colnames(ATAC), colnames(sobj))
ATAC <- AddMetaData(
  object   = ATAC,
  metadata = sobj@meta.data[common.cells, "Celltype", drop = FALSE]
)
# Plot ATAC UMAP by cell type
DimPlot(object = ATAC, group.by = "Celltype", label = TRUE) + NoLegend()

# Compute gene activity matrix from ATAC data
gene.activities <- GeneActivity(ATAC)
# Add gene activity as new assay
ATAC[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)

# Normalize and scale gene activity assay
DefaultAssay(ATAC) <- "ACTIVITY"
ATAC <- NormalizeData(ATAC)
ATAC <- ScaleData(ATAC)

# Find anchors for label transfer from RNA to ATAC
transfer.anchors <- FindTransferAnchors(
  reference = RNA,
  query = ATAC,
  reduction = 'cca'
)
# Transfer cell type labels
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA$Celltype,
  weight.reduction = ATAC[['lsi']],
  dims = 2:30
)
# Add predicted labels to ATAC metadata
ATAC <- AddMetaData(object = ATAC, metadata = predicted.labels)

# Plot comparison of RNA and ATAC cell type assignments
plot1 <- DimPlot(
  object = RNA,
  group.by = 'Celltype',
  label = TRUE,
  repel = TRUE
) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(
  object = ATAC,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE
) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2

# Save the processed ATAC object
saveRDS(ATAC, file = 'D:/Hackathon/Annotated_processed_unpaired_ATAC.RDS')
