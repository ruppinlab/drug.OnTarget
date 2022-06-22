# Find all the clinically approved drugs
DOI_cancer=onTarget$drugCategory$name[onTarget$drugCategory$drug_category=='targeted cancer']
FDAapproved_drugs=read.csv('2.Data/fda_approved_maybe_allDrugs.csv')
DOI_cancer_approved=na.omit(DOI_cancer[match(tolower(DOI_cancer), tolower(approved_drugs$name))])
onTarget$PredvsKnown_scores[onTarget$PredvsKnown_scores$CommonDrugName %in% DOI_cancer_approved,]
onTarget$PredvsKnown_scores