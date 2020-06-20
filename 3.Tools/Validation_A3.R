# Validation in another dataset - From Miquel dataset
gs2=read.csv('/Users/sinhas8/Downloads/drugbank_chembldrug.tsv', sep='\t')
gs2=gs2[gs2$direction == -1,]
gs2$inchikey
