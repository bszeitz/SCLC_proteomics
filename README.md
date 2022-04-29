# SCLC_proteomics

**This repository contains the custom R code used for the analyses reported in:**

Szeitz, B., Megyesfalvi, Z., Woldmar N., Valkó, Z., Schwendenwein, A., Bárány, N., Paku, S., László, V., Kiss, H., Bugyik, E., Lang, C., Szász, A.M., Pizzatti, L., Bogos, K., Hoda, M.A., Hoetzenecker, H., Marko-Varga, G., Horvatovich, P., Döme, B., Schelch, K., Rezeli, M. (2022). **In-depth proteomic analysis reveals unique subtype-specific signatures in human small cell lung cancer.** Manuscript submitted for publication. 

## Study summary

We report herein the distinct proteomic signatures of small cell lung cancer subtypes (SCLC-A/N/P/Y).

**Methods:**

In total, 26 patient-derived cell lines with predetermined transcriptomic subtype were subjected to label-free shotgun proteomics. Both the cell pellet (CP) and the cell media (CM) were analyzed. The nLC-MS/MS analysis was performed on Ultimate 3000 RSLCnano system coupled to Q Exactive HF-X with data-dependent acquisition, followed by database search in Proteome Discoverer 2.4.

**Details of the analyses performed in R:**

1. Proteomic data processing (**Proteomics_post-processing.Rmd**)

The data processing was performed separately for CP and CM samples. This includes: log2-transformation of protein intensities, median normalization and batch effect correction. Sample replicates were subsequently averaged and only proteins with min. 80% valid values across samples were kept. Missing values were imputed using imputation based on normal distribution in Perseus v.1.6.

2. Selection of secreted proteins (**Selection_of_secreted_proteins.Rmd**)

Three secretome databases were used to delineate secreted proteins among the identified proteins (relevant for the CM data). Proteins mentioned in at least two secretome databases were considered secreted in later analyses.

3. RNA-Seq data retrieval and processing (**RNA-Seq_data_retrieval_and_processing.Rmd**)

The FPKM and Z-scored values from George et al. [1] were retrieved, as well as the CCLE transcriptomic data (RSEM values) [2] were processed with limma (normalization, differential expression analysis between subtypes).

4. Differential expression analyses of the proteomic data (**Differential_expression_with_culture_type.Rmd**, **Differential_expression_with_subtype.Rmd**)

Differential expression analyses with subtypes and cell culture types (adherent, semi-adherent, in suspension) were conducted, as well as pathway enrichment analyses using the differentially expressed proteins.

5. Neuroendocrine and epithelial-mesenchymal transition scores (**NE_and_EMT_scores.Rmd**)

Neuroendocrine and epithelial-mesenchymal transition scores were calculated for each cell line based on the expression of curated marker proteins.

6. Consensus clustering and sources of variation (**Consensus_clustering_and_sources_of_variation.Rmd**)

Unsupervised consensus clustering on the CP dataset was performed, using proteins with SD>1.25. The main sources of variation in the CP and CM datasets were delineated using PCA and PVCA.

7. Pre-ranked gene set enrichment analysis (**Pre-ranked_GSEA.Rmd**)

The pGSEA was performed for all six subtype comparisons using the pre-ranked list of genes/proteins (*log2FC multiplied by -log10(p-value)*) from the proteomic CP and CCLE transcriptomic datasets. The results were further filtered for subtype-specific gene set enrichment.

8. Sparse partial least squares discriminant analysis (**sPLSDA.Rmd**)

The sPLS-DA was conducted using a combined CP+CM proteomic dataset to select proteins that show unique expression profile and are potential IHC/blood-based biomarker candidates for the subtypes.


*Scripts for manuscript figures:*

- **Figure_1.Rmd**
- **Figure_2_and_S4.Rmd**
- **Figure_3_and_4ab.Rmd**
- **Figure_4c_and_S8.Rmd**
- **Figure_S1.Rmd**
- **Figure_S2.Rmd**
- **Figure_S3.Rmd**
- **Figure_S9.Rmd**

*Helper scripts:*

- **Load_packages.R**: load all required packages at the beginning of the script (if package is missing, then it needs to be installed beforehand).
- **Utility_functions.R**: custom R functions used in the scripts.
- **Color_list.R**: specifying the colors for heatmap annotations.

See **Sessioninfo.txt** for details about the working environment in which the scripts were ran.

**Data availability:**

The proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifiers PXD029805 (cell pellet data) and PXD029821 (culture media data). 

**References:**

[1] George J, Lim JS, Jang SJ, et al. Comprehensive genomic profiles of small cell lung cancer. Nature 2015; 524: 47-53

[2] Ghandi M, Huang FW, Jané-Valbuena J, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. Nature 2019; 569: 503-508
