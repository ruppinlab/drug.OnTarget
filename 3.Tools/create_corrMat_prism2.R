# Create Corrmat_crispr_Secondary
secondary=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/secondary_screen_processed.RDS')
onTarget$secondary_prism=secondary$secondary_screen
onTarget=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v3.RDS')
###########################################################################
# Step 0: match cols and rows
###########################################################################
secondary_prism_matched=secondary$secondary_screen[,!is.na(match(secondary$cellLine_Info, colnames(onTarget$avana)))]
avana_matched=onTarget$avana[,match(colnames(secondary_prism_matched), colnames(onTarget$avana))]
###########################################################################
# Step 0: REplicating create_CorrMatrix for shRNA
###########################################################################
use_cores=detectCores()
# use_cores=4
avana_matched_list=lapply(split(avana_matched, c(rep(1:use_cores, each=nrow(avana_matched)/use_cores),
                                 rep(use_cores, nrow(avana_matched)%%use_cores))),
                  matrix, ncol=ncol(avana_matched))
calling_CORtest <- function(avana_matched_list_mat){
  lapply(1:nrow(avana_matched_list_mat), function(y) sapply(1:nrow(secondary_prism_matched), function(x)
    cor(avana_matched_list_mat[y,],
        secondary_prism_matched[x,],
        use = "pairwise.complete.obs",
        quick = 1,
        weights.x = NULL,
        weights.y = NULL)
  ))
}
cmat=mclapply(1:length(avana_matched_list), function(x) err_handle(calling_CORtest(avana_matched_list[[x]])) ,
              mc.cores = detectCores())
# cmat_backup=cmat
cmat=unlist(cmat, recursive = F)
cmat=do.call(rbind, cmat)
rownames(cmat)=rownames(onTarget$avana)
colnames(cmat)=rownames(secondary_prism_matched)
saveRDS(cmat, '/Users/sinhas8/Project_OffTarget/2.Data/corrMat_secondary_crispr.RDS')
# Pick shRNA corr matrix from here
