# integrative Glioblastoma Subtype classifier (iGlioSub)
Full paper accessible at: https://biodatamining.biomedcentral.com/articles/10.1186/s13040-021-00273-8

Glioblastoma (GBM) is the most aggressive and prevalent form of primary brain tumor, with a median survival of 15 months. Advancements in multi-omics profiling combined with computational algorithms have unraveled the existence of three GBM molecular subtypes (Classical, Mesenchymal, and Proneural) with clinical relevance. However, due to the costs of high-throughput profiling techniques, GBM molecular subtyping is not currently employed in clinical settings. Using machine learning, we constructed classifiers that require a minimum number of genomic features to overcome the limitations of efficiently classify GBM into clinically-relevant subtypes. These classification systems can be carried out using gene expression, DNA methylation (DNAm), or integration between these two types of data. The integrative GBM Subtype classifier (iGlioSub), constructed with only ten features per subtype, DNAm shows a superior performance (kappa=0.829) to stratifystratifying patients than gene expression (kappa=0.68) and DNAm-based classifiers (kappa=0.82). Importantly, the integrative classifiers constructed with only ten features per subtype have the highest classification efficiency (kappa=0.9). Also, expanding the understanding of the molecular differences between the GBM subtypes, this study reveals that each subtype presents unique DNAm patterns and gene pathway activation. These findings open an opportunity to explore novel subtype-specific complementary therapies. The data presented in this studyiGlioSub provides the basis to design cost-effective techniques to stratify GBM patients in routine pathology laboratories for clinical trials, which will significantly accelerate the discovery of more efficient GBM subtype-specific treatment approaches. 

**Order of the scripts**

1- Merge 27k-450k-GEO and Correct Batch Effect.R

2- Main Script - Create GeneExp, DNAm and Integrative Panels and Validate.R

The other two scripts (Error Rate and AUC Radar Plots.R and Create GSEA and GREAT Plots.R) must be used after the main script. Please remind that to create the GSEA and GREAT plots, the enriched pathways must be obtained using Metascape and GREAT tools.
