# ST_EPN overview
Source code Schüftan et al. 2025

This repository contains code developed in "Integrated dual-omics analysis elucidates transcriptional heterogeneity and hypoxia-driven spatial organisation in supratentorial ependymoma".

### Single-nucleus transcriptomic analyses  
Code used for the analysis of MACSima multiplex immunofluorescence cell segmentation data.  

- **01_Seurat_from_MACSima.R:** Importing cell segmentation data into R and generating respective Seurat objects. Each row in the input data corresponds to one cell, columns contain individual protein intensities and cell coordinates.  
- **02_Seurat_MACSima_analysis.R:** Classifying each cell into tumor/non-tumor based on whether it is located inside or outside a tumor cell area (and not whether it is a tumor cell itself). The input data format is identical to that of '01_Seurat_from_MACSima.R' but only contains cells in tumor cell areas based on previous gating in the MACSiQ View software. Lymphoid, myeloid and CD204+ cells get quantified according to marker expression.  
- **03_immune_quantification.R:** Comparison of the proportions of lymphoid, myeloid and CD204+ cells between tumor and non-tumor cell areas.  
- **marker_gating_f

### Spatial transcriptomic analyses  
- **01_Seurat_from_MACSima.R:** Importing cell segmentation data into R and generating respective Seurat objects. Each row in the input data corresponds to one cell, columns contain individual protein intensities and cell coordinates.  
