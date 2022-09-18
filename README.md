## DeepTarget: Deep mechanism of action of cancer drugs by integrating large-scale genetic and drug screens
 Goal: Our aim is to identify drug targets by integrating large-scale drug and genetic screens.
 
### Prerequisites
RStudio, R
Majorly Packages required: fgsea; data.table; parallel; ggplot2; pROC; tidyverse;

### Custom Function required for the pipeline:

> Step0_Write_Functions.Rmd

In this file, we defined all function which are required for the analysis and validation of our pipeline. 
This step is very important and imported to every file to generate figures.

> Step0B_Data_and_Libraries.Rmd

In this file, we have pre-loaded and preprocessed all the datasets (Other than Validation data) and libraries required for analysis, main Figure and supplementary figure.


### Data Curation File: 
Goto >> Tools >> Data_curation_code

#### Requred Files: For detailed method section, visit to Manuscript
We collected the viability screens after CRISPR-Cas9 and drug treatment from the DepMap database: https://depmap.org/portal/. 

---
1) Step_1: KnowTraget_prediction:

> Step_0A_data_curation_and saving.Rmd
* Required Data: 
  * OnTarget_v2.RDS
* Required File: 
  * Step0_Write_Functions.Rmd

In this Step, we first mined large-scale drug response screens (PRISM) for more than ~1500 cancer drugs and genome-wide CRISPR-Cas9 knockout viability profiles (CRISPR-KO essentiality, AVANA) from DepMap, which were commonly performed across 371 cancer cell lines. Integrating these two screens, we calculated a (drug, target) pair-wise similarity score (Pearson Correlation) between viability after drug treatment and that observed for each gene CRISPR-KO termed it the Drug-KO Similarity score (DKS score). This score can range from -1 to 1. Using this DKS Score, we curated data KnownTarget_predition which will be majorly used the futher analysis and validation of our Pipeline. 

The result of the correlation is saved in the Data folder:
> drugVScrispr_corr_features_list.RDS
> KnowTraget_prediction.RDS

Variable defined in the:
| Variable Name | Defination |
| -----------| ----------- | 
| MaxtargetName | If a drug have multiple target, what is the best target which have best correlation. |
| Maxcorr | The respective correlation stregth of the drug target pair. |

---
2) Step_2: KnowTraget_prediction: Extending the data of KnowTraget_prediction by including RNAexpression data.
> Step_0B_Data_curation_Interaction.Rmd
* Required Data: 
  * OnTarget_v2.RDS 
  * KnowTraget_prediction.RDS

* Required File: 
  * Step0_Write_Functions.Rmd

To this end, we first tested if indeed the correlation between the drug response and viability after the primary target CRISPR-KO (DKS score) decreases in cell lines where the primary target is not expressed. For this we computed the interaction between the expression data of the Cellines and DKS score and included these information which will be used when we will perform secondary target analysis and figure generation. 

The result of the correlation is saved in the Data folder:
> KnowTraget_prediction.RDS
---

3) Step_3: KnowTraget_prediction: Extending the data of KnowTraget_prediction by including Mutation data.
> Step_0C_Data_curation_Mutation_interaction.Rmd
* Required Data: 
  * OnTarget_v2.RDS 
  * KnowTraget_prediction.RDS

* Required File: 
  * Step0_Write_Functions.Rmd

We next identify whether a given drug more specifically targets the mutant or the wild-type form of its known target protein. For this we computed the interaction between the mutation data of the Cellines and DKS score and included these information which will be used when we will perform Mutation Analysis, Validation and figures.

The result of the correlation is saved in the Data folder:
> KnowTraget_prediction.RDS (Final)
---

4) Step_4: Creating correlation matrix for secondary target.
> create_CorrMat_secondary.Rmd
* Required Data: 
  * OnTarget_v2.RDS 
  * KnowTraget_prediction.RDS

* Required File: 
  * Step0_Write_Functions.Rmd

In this step, we repeated our target identification process (Step 1) in cell lines where their primary targets are not expressed to produce a secondary-DKS Score for such drugs. This saved dataset is used in secondary target analysis and figure generation.

The result of the correlation is saved in the Data folder:
> drugVScrispr_corr_features_list_secondary.RDS
---

5) Step_5: Pathway Enrichment
> create_CorrMat_secondary.Rmd
* Required library/package: 
  * library(fgsea)
  
* Required File: 
  * Step0_Write_Functions.Rmd
  * Step0B_Data_and_Libraries.Rmd (required dataset already loaded)

In this step, We are runing the pathway enrichment test across all the pathway (n=800) for all the drug (drugs from KnownTraget_prediction file) and save it in dataset gsea_enrichment_all_drugs, which will be used in Pathway enrichment test analysis and supplemntary figure generation.

Run gene enrichment function, it will take some time, ~20 mins. The result of the correlation is saved in the Data folder:
> moa_pathways.RDS

> gsea_enrichment_all_drugs.RDS
---

### Figures

To generate the figures:
> Goto >> Tools 
---
* Figure 1: Validation of our model: Overview of our pipeline and performance across eight gold-standard datasets.
> Figure_1_ValidationPrimaryTarget_AUC.RMD

In this step we are tested DeepTarget in eight independent gold standard datasets comprising high-confidence drug-target pairs collected from diverse sources (detailed defination in Manuscriprt, Method section). These datasets include drug-target pairs that have 
1) Clinical resistance mutation in the target from COSMIC (COSMIC resistance, N=16).
2) OncoKB (oncoKB resistance, N=28).
3) FDA Approval for a target mutation (FDA mutation-approval, N=86).
4) High-confidence as per the scientific advisory board of ChemicalProbes.org (SAB, N=24).
5) Multiple independent reports as per BioGrid (Biogrid Highly Cited, N=28).
6) Pharmacologically active status as per DrugBank, i.e., the drugs interacts directly with the target as part of its mechanism of action and are inhibitors (DrugBank Active Inhibitors, N=90).
7) DrugBank Active Antagonists (N=52).
8) Highly selective inhibitors based on their binding profile (SelecChem selective inhibitors, N=142).

We first defined the function to create randomize controls samples for testing. and tested the AUC for each of these 8 dataset. We have also ploted the data distribution (Box plot) and AUC plot for indiviisual dataset, which is shown in Figure 1 B -D.  

---





* Figure 2: Mutation Analysis and Validation

> Figure_2_Mutation_analysis.RMD

* Figure 3: Secondary Target Analysis and Validation

> Figure_3_and supp_12_secondary_target_and_validation.RMD

* Figure 4: Applications: Where we 

> Figure_3_and supp_12_secondary_target_and_validation.RMD




