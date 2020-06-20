## Example 2: get the signatures (i.e. level 5 data) for the drug vorinostat
source('/Users/sinhas8/myCustom_functions.R')
setwd('/Users/sinhas8/Project_OffTarget/2.Data/PhaseI/')
load("siginfo.RData")
library(cmapR)
library(data.table)
require(matrixStats)
table(siginfo$pert_type)
yaptaz_lincs1000=siginfo[which(siginfo$pert_iname=='WWTR1' | siginfo$pert_iname=='YAP1' ),]
nrow(yaptaz_lincs1000)/18
write.csv(yaptaz_lincs1000, '/Users/sinhas8/Project_YAPTAZ/2.Data/yaptaz_lincs1000.csv' )
###########################################################################
# Functions to get a CGS
###########################################################################
# unweighted_consensus_crispr_signature
unweighted_consensus_crispr_signature<-function(geneName='BRAF'){
  sig.ids <- siginfo[pert_iname==geneName & pert_type=='trt_xpr', sig_id]
  sigs <- parse_gctx("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                     cid=sig.ids)
  rowMedians(apply(sigs@mat, 2, scale))
}

weighted_consensus_crispr_signature<-function(geneName='BRAF'){
  sig.ids <- siginfo[pert_iname==geneName & pert_type=='trt_xpr', sig_id]
  sigs <- parse_gctx("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                     cid=sig.ids)
  cellLine_names=sapply(colnames(sigs@mat), function(x) strsplit(x, '_')[[1]][2])
}
# Drug Signature
consensus_drug_signature_modf<-function(drugName_ID='metformin', phase='1'){
  sig.ids <- siginfo[pert_iname==drugName, sig_id]
  
  if(phase==1){
    sigs <- parse_gctx("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                       cid=sig.ids)
  } else { 
    sigs <- parse_gctx("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                     cid=sig.ids)
    }
  by_hours=sapply(sigs@cid, function(x) strsplit(strsplit(x, '_')[[1]][3], ':')[[1]][1])
  myhead(split(data.frame(t(apply(sigs@mat, 2, scale))), by_hours)[[1]])
  lapply(split(data.frame(t(apply(sigs@mat, 2, scale))), by_hours), function(x) rowMedians(t(x)))
}

###########################################################################
# Cmat Preprocessing
###########################################################################
cmat=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/corrMat_matricForm.RDS')
colnames(cmat)=rownames(prism)
rownames(cmat)=rownames(crispr)
cmat_rank=apply(cmat, 2, function(x) rank(x) )

###########################################################################
# Genes with KO avaiable
###########################################################################
GOI <- unique(siginfo[pert_type=='trt_xpr', pert_iname])
top_targets_id=apply(cmat_rank, 2, function(x) which(x==1))
Respective_Targets=rownames(cmat_rank)[top_targets_id]
top_targets=unique(Respective_Targets)
overlap_genes=top_targets[na.omit(match(GOI, top_targets))]
df2use=data.frame(geneName=overlap_genes[na.omit(match(Respective_Targets, overlap_genes))],
                  drugName=colnames(cmat)[!is.na(match(Respective_Targets, overlap_genes))])
df2use$drugname_common=depmap_19Q3$drug_target_map$name[match(df2use$drugName, 
                                  depmap_19Q3$drug_target_map$broad_id_trimmed)]
df2use$annotated_Target=depmap_19Q3$drug_target_map$target[match(df2use$drugName, 
                                                               depmap_19Q3$drug_target_map$broad_id_trimmed)]
df2use$pred_corr= sapply(1:nrow(df2use), function(x) grepl(df2use$geneName[x], df2use$annotated_Target[x]) )
head(df2use[!sapply(1:nrow(df2use), function(x) grepl(df2use$geneName[x], df2use$annotated_Target[x]) ), c(1,4)])
###########################################################################
# Mapping Drug Name
###########################################################################
all_drugs_P2=unique(siginfo[pert_type=='trt_cp', pert_iname])
match(df2use$drugname_common, all_drugs_P2)
# Test Phase I as well
# all_drugs_P1=unique(siginfo[pert_type=='trt_cp', pert_iname])
# match(df2use$drugname_common, all_drugs_P1)
###########################################################################
# Get cgs
###########################################################################
csg_genes_cko=sapply(as.character(df2use$geneName), unweighted_consensus_crispr_signature)
csg_drugs_trt=lapply(as.character(df2use$drugname_common),
                     function(x) err_handle(consensus_drug_signature(x)))
csg_drugs_trt_P1=lapply(as.character(df2use$drugname_common),
                     function(x) err_handle(consensus_drug_signature(x, phase=1)))
cor_ko_vs_trt_P1=do.call(rbind, sapply(1:ncol(csg_genes_cko), function(x)
  err_handle(unlist(cor.test(csg_drugs_trt_P1[[x]][[1]], csg_genes_cko[,x])[c(3,4)])) ))
df2use[which(cor_ko_vs_trt[,2]>0.1),]
###########################################################################
# Phase 2
###########################################################################
csg_genes_cko=sapply(as.character(df2use$geneName), unweighted_consensus_crispr_signature)
csg_drugs_trt=lapply(as.character(df2use$drugname_common),
                     function(x) err_handle(consensus_drug_signature(x)))
cor_ko_vs_trt=do.call(rbind, sapply(1:ncol(csg_genes_cko), function(x)
  err_handle(unlist(cor.test(csg_drugs_trt[[x]][[1]], csg_genes_cko[,x])[c(3,4)])) ))
df2use[which(cor_ko_vs_trt[,2]>0.1),]

###########################################################################
# 
###########################################################################
plot(hist(colMaxs(cmat), breaks = 30))
sum(colMaxs(cmat)>0.3 & Annotated_Target_corr_max>0.3, na.rm = T)

