# Validation A9 - COSMICv91 Drug Resistance file
require(data.table)
setwd('/Users/sinhas8/Project_OffTarget/')
CosmicResistanceMutations=read.table('CosmicResistanceMutations.tsv', sep='\t')
colnames(CosmicResistanceMutations)=as.character(unlist(CosmicResistanceMutations[1,]))
head(sort(table(CosmicResistanceMutations$`Pubmed Id`), decreasing = T))

CosmicResistanceMutations$`Drug Name`