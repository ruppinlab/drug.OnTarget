# Categorizing Drugs
max_corr_Drug=colMaxs(cmat)
second_max_corr_Drug=apply(cmat, 2, function(x) max(x[x!=max(x)]) )
diff2Penultimate=max_corr_Drug - second_max_corr_Drug
plot(hist(max_corr_Drug - second_max_corr_Drug, 100))
plot(density(max_corr_Drug - second_max_corr_Drug))

# We propose here that single gene target drugs would have a max - sec_max to be greater than 0.05
singleTarget_Drugs=which(diff2Penultimate>0.05)
ateastOne_HighCorr_with_crispr_Drugs 
ateastOne_HighCorr_with_crispr_Drugs=which(colMaxs(cmat)>0.3)
polygenic_Drugs=ateastOne_HighCorr_with_crispr_Drugs[!ateastOne_HighCorr_with_crispr_Drugs %in% singleTarget_Drugs]
HighCorr_with_annotatedcrispr_Drugs = which(Annotated_Target_corr_max>0.3)
# High Concordant with at least one crispr are likely to be corrected
overlap<-function(list1, list2){
  sum(list1 %in% list2)*2/
  length(c(list1, list2))
}
overlap(HighCorr_with_annotatedcrispr_Drugs, singleTarget_Drugs)
###########################################################################
# Perfect annotation
###########################################################################
perfectly_predicted=which(Annotated_Target_corr_rank_mean== max(Annotated_Target_corr_rank_mean, na.rm = T) )
# Semi-Perfect annotation
semi_perfectly_predicted=which(Annotated_Target_corr_rank_mean>max(Annotated_Target_corr_rank_mean, na.rm = T)-100)
overlap(semi_perfectly_predicted, singleTarget_Drugs)
overlap(perfectly_predicted, polygenic_Drugs)

