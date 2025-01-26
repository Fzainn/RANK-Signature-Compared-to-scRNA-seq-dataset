# RANK-Signature-Compared-to-scRNA-seq-dataset
This document provides an overview of the analysis performed to evaluate the activity of the RANK signature in single-cell RNA sequencing (scRNA-seq) data from epithelial cells. The analysis includes calculating AUC scores for the RANK signature, visualizing the results, and comparing activity across cell clusters.

# Source Data
The scRNA-seq data used in this analysis is sourced from the CellxGene Data Portal:
* Collection : [Breast Cancer Single-Cell RNA-seq](https://cellxgene.cziscience.com/collections/4195ab4c-20bd-4cd3-8b3d-65601277e731)
* Dataset: Epithelial cells from breast cancer samples.
The data was downloaded in .h5ad format and processed using the Seurat and AUCell packages in R.

# Analysis Steps
1. Load RANK Signature Data: The RANK signature genes were loaded from a CSV file (filtered_data.csv).
2. Map Genes to Ensembl IDs: Gene symbols were mapped to Ensembl IDs using a GTF file (Homo_sapiens.GRCh38.113.gtf.gz), as RANK Signature in gene symbols format.
3. Subset scRNA-seq Data: the scRNA-seq data was subset to include only the RANK signature genes.
4. Calculate AUC Scores: AUC scores were calculated using the AUCell package to quantify the activity of the RANK signature in each cell.
5. Add AUC Scores to Seurat Object: The AUC scores were added to the Seurat object’s metadata for downstream analysis.
6. Visualize AUC Scores: A UMAP plot was generated to visualize the spatial distribution of RANK signature activity across cells and A histogram was created to show the distribution of AUC scores.
7. Compare AUC Scores Between Clusters: The mean AUC scores were calculated for each cluster to identify clusters with high RANK signature activity.
8. Save Results: The updated Seurat object with AUC scores was saved as subset_seurat_with_RANK_signature_AUC.rds and The UMAP plot and histogram were saved in a single PDF file (RANK_Signature_Analysis_Plots.pdf).

# Files Included 
1. Seurat Object:
   * File: subset_seurat_with_RANK_signature_AUC.rds
   * Description: Contains the scRNA-seq data with RANK signature AUC scores added to the metadata.

2. PDF File:
   * File: RANK_Signature_Analysis_Plots.pdf
   * Description: A histogram showing the distribution of RANK signature AUC scores and A UMAP plot showing the spatial distribution of RANK signature activity across cells.

3. Cluster Comparison Results:
   * Description: A table summarizing the mean AUC scores for each cluster (printed in the R script).
  
# Reproducing the Analysis 
To reproduce the analysis, follow these steps: 

1. Install Required Packages:
   * Seurat
   * AUCell
   * ggplot2
   * patchwork
   * rtracklayer
   * reticulate
2. Download the Data:
   * RANK signature genes from [here](https://www.ncbi.nlm.nih.gov/gds?linkname=pubmed_gds&from_uid=34004159)
   * Download the scRNA-seq data for epithelial cells from the [CellxGene Data Portal](https://cellxgene.cziscience.com/collections/4195ab4c-20bd-4cd3-8b3d-65601277e731) Save the data as scRNA-Epithelial.h5ad.
   * Annotation file for [Homo_sapiens.GRCh38.dna.alt.fa.gz](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/)





