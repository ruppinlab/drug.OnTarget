# Create a corrMat with shRNA dataset
require(parallel); require(tictoc); require(WGCNA); require(ggplot2); require(cowplot)
source('/Users/sinhas8/myCustom_functions.R')

shRNA=onTarget$achilles[,intersect(colnames(onTarget$achilles), colnames(onTarget$drug_prism))]
prism=onTarget$drug_prism[,intersect(colnames(onTarget$achilles), colnames(onTarget$drug_prism))]
###########################################################################
# Step 0: REplicating create_CorrMatrix for shRNA
###########################################################################
use_cores=detectCores()
# use_cores=4
shRNA_list=lapply(split(shRNA, c(rep(1:use_cores, each=nrow(shRNA)/use_cores),
                                   rep(use_cores, nrow(shRNA)%%use_cores))),
                   matrix, ncol=ncol(shRNA))

calling_CORtest <- function(shRNA_mat){
  lapply(1:nrow(shRNA_mat), function(y) sapply(1:nrow(prism), function(x)
    cor(shRNA_mat[y,],
        prism[x,],
        use = "pairwise.complete.obs",
        quick = 1,
        weights.x = NULL,
        weights.y = NULL)
  ))
}
tic()
cmat=mclapply(1:length(shRNA_list), function(x) err_handle(calling_CORtest(shRNA_list[[x]])) ,
            mc.cores = detectCores())
cmat=unlist(cmat, recursive = F)
cmat=do.call(rbind, cmat)
toc()
saveRDS(cmat, '/Users/sinhas8/Project_OffTarget/2.Data/corrMat_shRNA.RDS')
# Pick shRNA corr matrix from here
