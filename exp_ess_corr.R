# Expression vs Essentiality
COI=intersect(colnames(depmap$expression), colnames(depmap$avana))
GOI=intersect(rownames(depmap$expression), rownames(depmap$avana))

expvsEss_corr=mclapply(1:length(GOI), function(x) unlist(cor.test(depmap$expression[GOI[x], COI], depmap$avana[GOI[x], COI])[c(3, 4)]) )
expvsEss_corr=do.call(rbind, expvsEss_corr)
rownames(expvsEss_corr)=GOI
expvsEss_corr=data.frame(expvsEss_corr)

COSMIC=read.csv('/Users/sinhas8/APA_Adriana/COSMIC.tsv', sep='\t')
sum(expvsEss_corr$estimate.cor > 0.3, na.rm = T)

