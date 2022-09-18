## DeepTarget: Deep mechanism of action of cancer drugs by integrating large-scale genetic and drug screens
 Goal: Our aim is to identify drug targets by integrating large-scale drug and genetic screens.
 
### Prerequisites
RStudio, R
Majorly Packages required: fgsea; data.table; parallel; ggplot2; pROC; tidyverse;

### Data Curation File:
---
1) KnowTraget_prediction:
Goto >> Tools >> Data_curation_code

Step_0A_data_curation_and saving : In this Step, we compute the correlation between the crispr KO and the Drug response. Using the correlation strenth, the p value and rank, we created the variables such as 
* MaxtargetName: if a drug have multiple target, what is the best target which have best correlation.
* Maxcorr: The respective correlation stregth of the drug target pair.

---

We have initially downaloaded the 

---



### Figures

To generate the figures:

Figure 1: Validation of our model

Goto >> Tools>> Figure_1_ValidationPrimaryTarget_AUC.RMD

Figure 2: Mutation Analysis and Validation

Goto >> Tools>> Figure_2_Mutation_analysis.RMD

Figure 3: Secondary Target Analysis and Validation

Goto >> Tools>> Figure_3_and supp_12_secondary_target_and_validation.RMD

Figure 4: Applications: Where we 

Goto >> Tools>> Figure_3_and supp_12_secondary_target_and_validation.RMD




