# Analyse different dimensionality reduction methods to see which one preserves biology best across unpaired modality 
## 2025 OZ Single Cell Hackathon 
### Authors: Harry Mueller, Maria Nucera, Tayla Albertini, Chen Zhu

## - *STEP 1*

  Preprocessing scRNA-seq and scATAC-seq data
  [Code](https://github.com/harrymueller/hackathon_multiomic_embeddings/blob/main/scRNA_scATAC_Preprocessing.R)

  Converting to anndata

  Adding cell type annotation
  
## - *STEP 2*
  
  Integration with Seurat

  Testing different dimensionality reduction methods

## - *STEP 3*
    
  #### Computing metrics to evaluate biological conservation

  ##### Indipendent from labels (no annotation required just basic knowledge of gene expression patterns in the system)

  Identify set of markers genes related to a cell type

  The idea is that those markers should be localised together in the embedding as they are spefic to a cell type

  The expression of those genes should not be sparse and not localised in 2 ore more neighborhood distant between each other

  - Moran’s I
    [doc](https://scanpy.readthedocs.io/en/stable/generated/scanpy.metrics.morans_i.html)
    
    What it measures:
    The global spatial autocorrelation of a variable (e.g., a gene’s expression) "how do similar values cluster together across the embedding?"
  - Geary’s C
    [doc](https://scanpy.readthedocs.io/en/stable/generated/scanpy.metrics.gearys_c.html)

    What it measures:
   The local spatial heterogeneity " how different each cell’s value is from its neighbors"

    | Metric        | Measures…           | Good for…                   | Ideal Value |
    | ------------- | ------------------- | --------------------------- | ----------- |
    | **Moran’s I** | Global similarity   | Coherent gene expression    | Near **+1** |
    | **Geary’s C** | Local dissimilarity | Patchiness or abrupt shifts | Near **0**  |

  - Custom function to compute spearman correlation between set of genes in close cells, far cells and random cells


   [Code](https://github.com/harrymueller/hackathon_multiomic_embeddings/blob/main/Biological_conservation_metrics%20(2).ipynb)


##### Using the cell types labels/annotation  (no clustering required) :

| **Metric**                               | **Measures…**                                      | **Good for…**                                                  | **Ideal Value**               |
| ---------------------------------------- | -------------------------------------------------- | -------------------------------------------------------------- | ----------------------------- |
| **Silhouette Score**                     | Separation of known cell types in embedding        | Detecting whether cell types remain distinct after integration | Near **+1**                   |
| **kNN Purity**                           | Fraction of a cell’s neighbors with the same label | Checking local consistency of cell type or modality            | Near **1.0**                  |
| **LISI (Local Inverse Simpson’s Index)** | Diversity of labels in local neighborhoods         | Assessing batch mixing (*iLISI*) or conservation (*cLISI*)     | High (*iLISI*), Low (*cLISI*) |


##### Using the cell types labels/annotation and clustering:

| **Metric**                              | **Measures…**                                     | **Good for…**                                    | **Ideal Value**   |
| --------------------------------------- | ------------------------------------------------- | ------------------------------------------------ | ----------------- |
| **Adjusted Rand Index (ARI)**           | Similarity between clustering and true labels     | Evaluating how well clustering recovers biology  | Near **+1**       |
| **NMI (Normalized Mutual Information)** | Overall agreement between labels and clusters     | Evaluating information shared between partitions | Near **1.0**      |



  - How to test for overintegration in RNA + ATAC integration
  Before integration, you remove one cell type (e.g. T cells) from the RNA data.
  So now your dataset has:

  All cell types in ATAC

  All cell types except T cells in RNA

  Then you run your integration method.
  After integration, you look at the cells that were originally T cells in the ATAC data.

   You check:
   Are the nearest neighbors of those ATAC T cells mostly RNA cells?

  If yes →  Overintegration
  The method falsely matched ATAC T cells with unrelated RNA cells
  If no (they’re mostly surrounded by other ATAC cells or stay isolated) →  Good integration
