# SCLC_proteomics

**This repository contains the custom R code used for the analyses reported in:**

**Szeitz, B., Megyesfalvi Z. et al. (2022). In-depth proteomic analysis reveals unique subtype-specific signatures in human small cell lung cancer.** 


**List of the scripts:**

*Helper functions:*

- **Install_load_packages.R**: load all required packages at the beginning of the script (if package is missing, then it is installed beforehand).
- **Utility_functions.R**: custom R functions used in the scripts.
- **Color_list.R**: specifying the colors for heatmap annotations.

*Analyses:*

- **Proteomics_post-processing.Rmd**: post-processing of the proteomic data (normalization, batch effect correction, filters).
- **Selection_of_secreted_proteins.Rmd**: using secretome databases to delineate secreted proteins among the identified proteins.
- **RNA-Seq_data_retrieval_and_processing.Rmd**: cleaning the transcriptomic data by George et al. [1], as well as processing of CCLE transcriptomic data [2] (normalization, differential expression analysis).
- **Differential_expression_with_culture_type.Rmd**: differential expression analysis with culture type.
- **Differential_expression_with_subtype.Rmd**: differential expression analysis with subtype.
- **NE_and_EMT_scores.Rmd**: calculate neuroendocrine and epithelial-mesenchymal transition scores.
- **Consensus_clustering_and_sources_of_variation.Rmd**: perform consensus clustering and PCA/PVCA analyses.
- **Pre-ranked_GSEA.Rmd**: perform pre-ranked gene set enrichment analyses.
- **sPLSDA.Rmd**: perform sparse partial least squares discriminant analysis.

*Figures:*

- **Figure_1.Rmd**
- **Figure_2_and_S4.Rmd**
- **Figure_3_and_4ab.Rmd**
- **Figure_4c_and_S8.Rmd**
- **Figure_S1.Rmd**
- **Figure_S2.Rmd**
- **Figure_S3.Rmd**
- **Figure_S9.Rmd**

See **Sessioninfo.txt** for details about the working environment in which the scripts were ran.

**References:**

[1] George J, Lim JS, Jang SJ, et al. Comprehensive genomic profiles of small cell lung cancer. Nature 2015; 524: 47-53

[2] Ghandi M, Huang FW, Jané-Valbuena J, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. Nature 2019; 569: 503-508
