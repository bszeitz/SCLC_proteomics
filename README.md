# SCLC_proteomics

**This repository contains the custom R code used for the analyses reported in:**

**Szeitz B, Megyesfalvi Z, Woldmar N, et al. In‐depth proteomic analysis reveals unique subtype‐specific signatures in human small‐cell lung cancer. Clin Transl Med. 2022;12:e1060. 10.1002/ctm2.1060**

## Study summary

We report herein the distinct proteomic signatures of small cell lung cancer subtypes (SCLC-A/N/P/Y).

**Methods:**

In total, 26 patient-derived cell lines with predetermined transcriptomic subtype were subjected to label-free shotgun proteomics. Both the cell pellet (CP) and the cell media (CM) were analyzed. The nLC-MS/MS analysis was performed on Ultimate 3000 RSLCnano system coupled to Q Exactive HF-X with data-dependent and data-independent acquisition, followed by database search in Proteome Discoverer 2.4.

**Short description of the analyses performed in R:**

1. Proteomic data processing (**Proteomics_post-processing.Rmd**)

The data processing was performed separately for CP and CM samples. This includes: log2-transformation of protein intensities, median normalization and batch effect correction. Sample replicates were subsequently averaged and only proteins with min. 80% valid values across samples were kept. Missing values were imputed using imputation based on normal distribution in Perseus v.1.6.

2. Annotation of proteins (**Select_secreted_surface_plasma_druggable_proteins.Rmd**)

Previously established databases were used to annotate proteins as secreted [1,2,3], cell-surface [4,5], human blood plasma [6], and druggable [7] proteins.

3. Correlation with recently published proteomic dataset (**Correlation_with_ProCan-DepMapSanger_dataset.Rmd**)

The protein expression table from Gonçalves et al. [8] was retrieved to check correlation of protein abundances across the 18 commonly measured cell lines.

4. RNA-Seq data retrieval and processing (**RNA-Seq_data_retrieval_and_processing.Rmd**)

The FPKM and Z-scored values from George et al. [9] were retrieved and processed for ssGSEA, as well as the CCLE transcriptomic data (RSEM values) [10] were processed with limma (normalization, differential expression analysis between subtypes).

5. Differential expression analyses of the proteomic data (**Differential_expression_with_culture_type.Rmd**, **Differential_expression_with_chemo_status.Rmd**, **Differential_expression_with_subtype.Rmd**)

Differential expression analyses with subtypes, chemotherapy status (post-chemo vs chemo-naive) and cell culture types (adherent vs in suspension) were conducted, followed by pathway enrichment analyses using the DE proteins.

6. Neuroendocrine and epithelial-mesenchymal transition scores (**NE_and_EMT_scores.Rmd**)

Neuroendocrine and epithelial-mesenchymal transition scores were calculated for each cell line based on the relative expression of curated marker proteins.

7. Consensus clustering and sources of variation (**Consensus_clustering_and_sources_of_variation.Rmd**)

Unsupervised consensus clustering on the CP dataset was performed, using proteins with SD > 1.25. The main sources of variation in the CP and CM datasets were delineated using PCA and PVCA.

8. Pre-ranked gene set enrichment analysis (**Pre-ranked_GSEA.Rmd**)

The pGSEA was performed for all six subtype comparisons using the pre-ranked list of genes/proteins (*log2FC multiplied by -log10(p-value)*) from the proteomic CP and CCLE transcriptomic datasets. The results were further filtered for subtype-specific gene sets.

9. Single-sample gene set enrichment analysis (**Create_gct_for_ssGSEA_and_results.Rmd**)

A gct file with subtype-specific gene sets was prepared for ssGSEA. The ssGSEA was ran using the script from Krug et al. [11]. The ssGSEA results were then loaded in here and were processed further.

10. Sparse partial least squares discriminant analysis (**sPLSDA_CP.Rmd**, **sPLSDA_CM.Rmd** and **sPLSDA_results.Rmd**)

The sPLS-DA was conducted separately for the CP and CM proteomic dataset (sPLSDA_CP.Rmd and sPLSDA_CM.Rmd). The results are summarized in sPLSDA_results.Rmd. 
After selecting proteins that are potential IHC/blood-based biomarker candidates for the subtypes, the tissue transcriptomic dataset [9] was utilized to confirm our findings.

11. Drug response data retrieval and analysis (**Cancerrxgene_dataset_drug_response.Rmd**)

Drug response data [12] was cleaned, and selected drugs were analyzed. The goal was to see whether the subtypes show any significant differences in their drug response.


*Scripts for manuscript figures (others were directly exported from the scripts above):*

- **Figure_1_molecular_heterogeneity.Rmd**
- **Figure_2_culture_type.Rmd**
- **Figure_3_and_S4_unsupervised_clustering.Rmd**
- **Figure_4ace_and_5a_cnetplots.Rmd**
- **Figure_S1_batch_correction.Rmd**
- **Figure_S2_known_markers.Rmd**
- **Figure_S3_NE_EMT.Rmd**

*Helper scripts:*

- **Load_packages.R**: load all required packages at the beginning of the script (if package is missing, then it needs to be installed beforehand).
- **Utility_functions.R**: custom R functions used in the scripts.
- **Color_list.R**: specifying the colors for heatmap annotations.

See **Sessioninfo.txt** for details about the working environment in which the scripts were ran.

## Data availability

The proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifiers PXD029805 (cell pellet data) and PXD029821 (culture media data). 

## Important references

[1] Uhlén M, Karlsson MJ, Hober A, et al. The human secretome. Sci Signal 2019; 12: eaaz0274

[2] Chen G, Chen J, Liu H, et al. Comprehensive Identification and Characterization of Human Secretome Based on Integrative Proteomic and Transcriptomic Data. Front Cell Dev Biol 2019; 7: 299

[3] Meiken J, Walker G, Cooper CR, Min XJ. MetazSecKB: the human and animal secretome and subcellular proteome knowledgebase. Database (Oxford) 2015; 2015: bav077

[4] Hu Z, Yuan J, Long M, et al. The Cancer Surfaceome Atlas integrates genomic, functional and drug response data to identify actionable targets. Nat Cancer 2021; 2: 1406-1422

[5] Bausch-Fluck D, Goldmann U, Müller S, et al. The in silico human surfaceome. Proc Natl Acad Sci. 2018; 115: E10988-E10997

[6] https://www.proteinatlas.org/humanproteome/blood+protein

[7] https://www.proteinatlas.org/humanproteome/tissue/druggable

[8] Gonçalves E, Poulos RC, Cai Z, et al. Pan-cancer proteomic map of 949 human cell lines. Cancer Cell 2022; 40: P835-849.E8

[9] George J, Lim JS, Jang SJ, et al. Comprehensive genomic profiles of small cell lung cancer. Nature 2015; 524: 47-53

[10] Ghandi M, Huang FW, Jané-Valbuena J, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. Nature 2019; 569: 503-508

[11] https://github.com/broadinstitute/ssGSEA2.0

[12] Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Res 2013; 41: D955-D961.
