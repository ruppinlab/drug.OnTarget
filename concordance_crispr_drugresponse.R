# Concordance btw drugs and crispr
##########################################################################
# Non-oncology drugs
##########################################################################
source('/Users/sinhas8/Project_OffTarget/Root_Step0.R')
concordance<-function(targetGene=mapping$target[1]){
  targetGene=as.character(targetGene)
  genetic_screen_profile=unlist(avana_raw_matchedresp[targetGene,])
  drug_screen_profile=unlist(resp_matchedavana_raw[match(targetGene,mapping$target),])
  unlist(cor.test(genetic_screen_profile, drug_screen_profile)[c(3,4)])
}

corr_score_btw_screengs=mclapply(1:nrow(mapping), function(x) err_handle(concordance(mapping$target[x])), mc.cores=4)
corr_score_btw_screengs=do.call(rbind, corr_score_btw_screengs)
sum(fdrcorr(corr_score_btw_screengs[,1])<0.2 & corr_score_btw_screengs[,2]>0.1, na.rm = T)


mapping$name[grep('03758309',mapping$name, ignore.case = T)]
corr_score_btw_screengs[grep('1776',mapping$name, ignore.case = T),]
