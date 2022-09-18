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
1) KnowTraget_prediction:
Goto >> Tools >> Data_curation_code

* Required Data: OnTarget_v2.Rdata (Data folder)
* Required File: Step0_Write_Functions.Rmd ## all function required for the project is defined here.

Step_0A_data_curation_and saving : In this Step, we compute the correlation between the crispr KO and the Drug response. Using the correlation strenth (Pearson Correlation Rho), the p value and rank, we curated data KnownTarget_predition which will be majorly used the futher analysis and validation of our Pipeline. 

The result of the correlation is saved in the Data folder:
>> drugVScrispr_corr_features_list.RDS
---



### Figures

To generate the figures:

Figure 1: Validation of our model



* Goto >> Tools>> Figure_1_ValidationPrimaryTarget_AUC.RMD

Figure 2: Mutation Analysis and Validation

* Goto >> Tools>> Figure_2_Mutation_analysis.RMD

Figure 3: Secondary Target Analysis and Validation

* Goto >> Tools>> Figure_3_and supp_12_secondary_target_and_validation.RMD

Figure 4: Applications: Where we 

* Goto >> Tools>> Figure_3_and supp_12_secondary_target_and_validation.RMD




