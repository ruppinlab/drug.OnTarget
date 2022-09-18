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

#### Requred Files: For detailed method section, visit to Manuscript
We collected the viability screens after CRISPR-Cas9 and drug treatment from the DepMap database: https://depmap.org/portal/. 

---
1) KnowTraget_prediction: Goto >> Tools >> Data_curation_code

> Step_0A_data_curation_and saving.Rmd

* Required Data: OnTarget_v2.RDS
* Required File: Step0_Write_Functions.Rmd ## all function required for the project is defined here.

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
Goto >> Tools >> Data_curation_code

* Required Data: 
  * OnTarget_v2.RDS 
  * KnowTraget_prediction.RDS

* Required File: Step0_Write_Functions.Rmd

To this end, we first tested if indeed the correlation between the drug response and viability after the primary target CRISPR-KO (DKS score) decreases in cell lines where the primary target is not expressed. For this we computed the interaction between the expression data of the Cellines and DKS score and included these information which will be used when we will perform secondary target analysis and figure generation. 

The result of the correlation is saved in the Data folder:
> KnowTraget_prediction.RDS
---

3) Step_3: KnowTraget_prediction: Extending the data of KnowTraget_prediction by including Mutation data.
Goto >> Tools >> Data_curation_code

* Required Data: 
  * OnTarget_v2.RDS 
  * KnowTraget_prediction.RDS

* Required File: Step0_Write_Functions.Rmd

We next identify whether a given drug more specifically targets the mutant or the wild-type form of its known target protein. For this we computed the interaction between the mutation data of the Cellines and DKS score and included these information which will be used when we will perform Mutation Analysis, Validation and figures.

The result of the correlation is saved in the Data folder:
> KnowTraget_prediction.RDS (Final)
---

### Figures

To generate the figures:
> Goto >> Tools 

* Figure 1: Validation of our model



> Figure_1_ValidationPrimaryTarget_AUC.RMD

* Figure 2: Mutation Analysis and Validation

> Figure_2_Mutation_analysis.RMD

* Figure 3: Secondary Target Analysis and Validation

> Figure_3_and supp_12_secondary_target_and_validation.RMD

* Figure 4: Applications: Where we 

> Figure_3_and supp_12_secondary_target_and_validation.RMD




