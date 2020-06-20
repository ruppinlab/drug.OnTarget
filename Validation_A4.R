# Validation A4: Via protein-drug affinity
require(data.table)
require(tictoc)
stitch=fread('/Users/sinhas8/Project_OffTarget/2.Data/9606.protein_chemical.links.detailed.v5.0.tsv',
             sep='\t')
geneMappping=fread('/Users/sinhas8/Downloads/ensp_to_genesymbol.txt', sep='\t')
stitch_proteinMapping=match(substring(stitch$protein, 6), geneMappping$Protein.stable.ID)
chemical_alias=fread('/Users/sinhas8/Project_OffTarget/2.Data/chemical.aliases.v5.0.tsv', sep='\t', nrows=10)

chemical_aliasv2=fread('/Users/sinhas8/Project_OffTarget/2.Data/chemical.aliases.v5.0.tsv',
                     sep='\t', drop=colnames(chemical_alias)[c(1,2,4)])

p53_drugs_effSize_delThis=read.csv('/Users/sinhas8/p53_drugs_effSize_delThis.csv')
head(onTarget$drugCategory[match(p53_drugs_effSize_delThis$X, onTarget$drugCategory$name),])

load('/Users/sinhas8/Downloads/drugs.of.interest.RData')
