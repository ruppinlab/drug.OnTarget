# read P53 mutation status from William Hahm NatGenetic paper
require('xlsx')
p53_status=readxl::read_xlsx('/Users/sinhas8/Downloads/NIHMS1501645-supplement-3.xlsx', skip = 1)
COI=na.omit(onTarget$annotation$depMapID[match(p53_status$CCLE_Cell_Line_Name, onTarget$annotation$CCLE_ID)])
p53_status=p53_status[!is.na(match(p53_status$CCLE_Cell_Line_Name, onTarget$annotation$CCLE_ID)),]
p53_status$Kras_mutation_status=onTarget$mutations_matrix['KRAS',COI]

p53_status$CellLines_ofDesiredp53andKRAS_Status=(p53_status$Kras_mutation_status==1 &
  p53_status$Genetic_and_Functional_p53_Status=='Wild-type_but_Non-functional')
p53_status=p53_status[order(p53_status$CellLines_ofDesiredp53andKRAS_Status, decreasing = T),]
p53_status$Kras_mutation_status=factor(p53_status$CellLines_ofDesiredp53andKRAS_Status, labels = c('WT', 'Mut'))
df2send$`Rank_Viability difference (CRISPR-shRNA) [Corrected]`=df2send$`Viability difference (CRISPR-shRNA) [Corrected]`
p53_status$`Viability difference (CRISPR-shRNA) [Corrected]`=df2send$`Viability difference (CRISPR-shRNA) [Corrected]`[
  match(p53_status$CCLE_Cell_Line_Name, rownames(df2send))]
p53_status$`Rank_Viability difference (CRISPR-shRNA) [Corrected]`=rank(p53_status$`Viability difference (CRISPR-shRNA) [Corrected]`, na.last = 'keep')
p53_status$`Rank_Viability difference (CRISPR-shRNA) [Corrected]`[p53_status$CellLines_ofDesiredp53andKRAS_Status]
p53_status[,c('Viability difference (CRISPR-shRNA) [Corrected]','CCLE_Cell_Line_Name')]

p53_status$Candidate_for_Experiment = p53_status$`Viability difference (CRISPR-shRNA) [Corrected]` > 0 &
  p53_status$CellLines_ofDesiredp53andKRAS_Status 

write.csv(p53_status, '/Users/sinhas8/Project_CRISPR/2.Data/p53_status_Hahn.csv')
