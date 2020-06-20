# Step 1: Validation via showing that CRISPR based model Likelihood Based model of prediction accuracy 
source('/Users/sinhas8/Project_OffTarget/Root_Step0.R')
source('/Users/sinhas8/myCustom_functions.R')
dim(avana_raw_matchedresp)
dim(resp_matchedavana_raw)

# Filter 1: We are only interested inhibitors
sum(grepl('inhibitor', mapping$moa, ignore.case = T))/length(grepl('inhibitor', mapping$moa, ignore.case = T))
mapping_inhibitors=mapping[grep('inhibitor', mapping$moa, ignore.case = T),]
resp_matchedavana_raw_inhibitors=resp_matchedavana_raw[grep('inhibitor', mapping$moa, ignore.case = T),]

# Filter 2: Drugs With multiple Target
mapping_inhibitors_MT=mapping_inhibitors[grep(',',mapping_inhibitors$target),]
resp_matchedavana_raw_inhibitors_MT=resp_matchedavana_raw_inhibitors[grep(',',mapping_inhibitors$target),]

table(mapping_inhibitors_MT$phase)

LBM_model<-function(id=1){
  Drug_Name=mapping_inhibitors_MT$broad_id[id]
  Targets=mapping_inhibitors_MT$target[id]
  Targets=unlist(strsplit(as.character(Targets), ',\\ ')[[1]])
  sapply(Targets)
  
}

