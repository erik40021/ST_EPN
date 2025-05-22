# Source code for Schüftan et al.

This repository contains all code developed for the study "Integrated dual-omics analysis elucidates transcriptional heterogeneity and hypoxia-driven spatial organisation in supratentorial ependymoma" (unrevised, unpublished).

### 1. Single-nucleus analyses  
- *01_processing.R:* quality control, filtering, integration and clustering of raw single-nucleus data
- *02_clustering_analysis.R:* annotation, subtype variance calculation and subsetting of malignant nuclei
- *03_CNA_inference.R:* inference of copy-number alterations (CNAs)
- *04_NMF_programs.R:* calculation and clustering of NMF consensus programs
- *05_program_evaluation.R:* evaluation of consensus programs regarding sample-wise expression and clinical favourability

### 2. Spatial transcriptomic analyses  
- *01_processing.R:* quality control, filtering and transformation of spatial samples
- *02_NMF_programs.R:* calculation and clustering of NMF consensus programs adapted to spatial data
- *03_CNA_inference.R:* inference of copy-number alterations (CNAs)
- *04_spatial_complexity.R:* calculation and visualisation of spatial coherence and spatial complexity scores, and definition of structural zones
- *05_spatial_associations.R:* calculation and visualisation of spatial association scores, definition of structural zones, and evaluation of spatial zones regarding immune programs
- *06_higher_order_model.R:* unbiased construction of higher-order graph model based on association scores and visualisation in Cytoscape

### 3. Immune analyses
- *01_processing.R:* integration and clustering of immune single-nucleus data
- *02_clustering_analysis.R:* annotation of clusters and transfer of NMF consensus programs
- *03_cell_cell_interactions.R:* cell-cell interaction analysis between hypoxic and non-hypoxic tumour nuclei and immune nuclei

### 4. Methylation analyses
- *01_clustering.R:* multi-iteration clustering of methylation data from 49 patients, integrated with 230 reference methylomes
- *02_CDKN2AB.R:* automated detection of coy-number alteration at the CDKN2A/B locus

### Utility functions
Collection of all utility functions needed to run the analyses above.



## General notes
Raw and processed data used in this study can be downloaded under EGA repository (coming soon).

NMF analysis was mainly based on Gavish et al. 23 (see https://doi.org/10.1038/s41586-023-06130-4 for additional reference).

Spatial analyses were conceptually inspired by Greenwald et al. 24 (see https://doi.org/10.1016/j.cell.2024.03.029 for additional reference).

R version used was 4.4.1.


