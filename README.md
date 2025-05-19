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

  The idea is that those markers should be localised together in the embedding as they are spefic to a cell type

  The expression of those genes should not be sparse and not localised in 2 ore more neighborhood distant between each other
  
  Computing metrics to evaluate biological conservation
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


   [Code]((https://github.com/harrymueller/hackathon_multiomic_embeddings/blob/main/Biological_conservation_metrics%20(1).ipynb))



  - Removing one cell type from RNA-seq before integration, and check if there is ATAC-seq in the neighbors of that cell type (there should not be, sign of overintegration)
