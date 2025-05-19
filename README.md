# Analyse different dimensionality reduction methods to see which one preserves biology best across unpaired modality 
## 2025 OZ Single Cell Hackathon 
### Authors: Harry Mueller, Maria Nucera, Tayla Albertini, Chen Zhu

- *STEP 1*

  Preprocessing scRNA-seq and scATAC-seq data

  Converting to anndata

  Adding cell type annotation
  
- *STEP 2*
  
  Integration with Harmony

  Testing different dimensionality reduction methods

- *STEP 3*

  Identify set of markers genes related to a cell type

  Computing metrics to evaluate biological conservation
  - Moran’s I (computes if gene expression is localised, but it will not distringuish between 2 disting neighborhoods far away from each other)
  - Geary’s C
  - Custom function to compute spearman correlation between set of genes in close cells, far cells and random cells
  - Removing one cell type from RNA-seq before integration, and check if there is ATAC-seq in the neighbors of that cell type (there should not be, sign of overintegration)
