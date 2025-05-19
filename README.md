# Analyse different dimensionality reduction methods to see which one preserves biology best across unpaired modality (e.g. CCA, WNN, UMAP, multiVI)
## 2025 OZ Single Cell Hackathon 
### Authors: Harry Mueller, Maria Nucera, Tayla Albertini, Chen Zhu

- *STEP 1*

  Preprocessing scRNA-seq and scATAC-seq data

  Converting to anndata

- *STEP 2*
  
  Integration with Harmony

- *STEP 3*

  Testing different dimensionality reduction methods

- *STEP 4*

  Computing metrics to evaluate biological conservation
  - Moran’s I (computes if gene expression is localised, but it will not distringuish between 2 disting neighborhoods far away from each other)
  - Gini (compare if a gene is expressed just by a few cells, but does not see where those cell are, could still be sparse)
  - Custom function to compute spearman correlation between set of genes in close cells, far cells and random cells
  - Removing one cell type from RNA-seq before integration, and check if there is ATAC-seq in the neighbors of that cell type (there should not be, sign of overintegration)
