# Source code for Schüftan et al. 2025

This repository contains all code developed for the study "Integrated dual-omics analysis elucidates transcriptional heterogeneity and hypoxia-driven spatial organisation in supratentorial ependymoma".

### Single-nucleus analyses  

- **01_processing.R:** quality control, filtering, integration and clustering of raw single-nucleus data
- **02_clustering_analysis.R:** annotation, subtype variance calculation and subsetting of malignant nuclei
- **03_CNA_inference.R:** inference of copy-number alterations (CNAs)
- **04_NMF_programs.R:** calculation and clustering of NMF consensus programs (see Gavish et al. https://doi.org/10.1038/s41586-023-06130-4 for reference)
- **05_program_evaluation.R:** evaluation of consensus programs regarding sample-wise expression and clinical favourability

### Spatial transcriptomic analyses  

- **01_processing.R:** quality control, filtering and transformation of spatial samples
- **02_NMF_programs.R:** calculation and clustering of NMF consensus programs adapted to spatial data
- **03_CNA_inference.R:** inference of copy-number alterations (CNAs)
- **04_spatial_complexity.R:** calculation and visualisation of spatial coherence and spatial complexity scores, and definition of structural zones
- **05_spatial_associations.R:** calculation and visualisation of spatial association scores, definition of structural zones, and evaluation of spatial zones regarding immune programs
- **06_higher_order_model.R:** unbiased construction of higher-order graph model based on association scores and visualisation in Cytoscape

### Immune analyses

- **01_processing.R:** integration and clustering of immune single-nucleus data
- **02_clustering_analysis.R:** annotation of clusters and transfer of NMF consensus programs
- **03_cell_cell_interactions.R:** cell-cell interaction analysis between hypoxic and non-hypoxic tumour nuclei and immune nuclei


### Utility functions
Collection of all utility functions needed to run the analyses above.


