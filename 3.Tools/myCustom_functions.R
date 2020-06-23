###Following are custom functions I've made to avoid repetitively writing same things again - after three years of my PhD :p
myhead<-function(x){
  x[1:min(5, nrow(x)), 1:min(5, ncol(x))]
}
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
test_topbottomQuantile<-function(var1=UCD_score_M1, var2_tobetiled=unlist(Exp_metab[1,]), which_tail='g', numofQuantiles=3){
  length(var1); length(var2_tobetiled)
  var2_tobetiled=xtile(var2_tobetiled, numofQuantiles)
  ##Only keep top/bottom quantile
  var1= var1[var2_tobetiled== min(var2_tobetiled) | var2_tobetiled== max(var2_tobetiled)]
  var2_tobetiled= var2_tobetiled[var2_tobetiled== min(var2_tobetiled) | var2_tobetiled== max(var2_tobetiled)]
  c(sig=err_handle(wilcox.test(var1 ~ factor(var2_tobetiled), alternative=which_tail)$p.value),
    eff_size=err_handle(diff(aggregate(var1, by=list(var2_tobetiled), mean)[,2]) ) )
}
hypergeometric_test_for_twolists<-function(test_list, base_list, global, lowertail=FALSE) {
  #If lowertail=FALSE - we calculate the probability for enrichment
  length(base_list)
  adj_base_list=global[na.omit(match(base_list, global))]
  Matched_list=test_list[!is.na(match(test_list, adj_base_list))]
  phyper(length(Matched_list)-1, length(adj_base_list), length(global)- length(adj_base_list), length(test_list), lower.tail=lowertail)
}
fdrcorr<-function(test_list){p.adjust(test_list, method = 'fdr')}

# dftemp=data.frame(prob$samples, prob$types, prob$stage)
# sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUSC',], function(x) length(x))[1:3,3])/sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUSC',], function(x) length(x))[,3])
# sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUAD',], function(x) length(x))[1:4,3])/sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUAD',], function(x) length(x))[,3])
#split(aggregate(mat$X ~ mat$stage_trimmed+mat$hist, data=dftemp, function(x) length(x)))

#Subsetting a set of columns
colSubset<-function(mat, column_Names){
  mat[,na.omit(match(column_Names, colnames(mat)))]
}
#Subsetting a set of rows
rowSubset<-function(mat, row_Names){
  mat[na.omit(match(row_Names, rownames(mat))),]
}

vectorSubset<-function(vec, Names){
  vec[!is.na(match(names(vec), Names))]
}


installORload<-function(packages){
  package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}
# Scale 0-1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# Scale 0-1 cancer type specifically
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}

# Create a random df unique to the input df
randomize_DF_with_unique_Pairs<-function(df2randomize=Positive_Set_DrugGene_Pairs){
  df2randomize_collapsed=paste(df2randomize[,1], df2randomize[,2], sep = '_')
  fixedCol=df2randomize[,2]
  shuffledCol=sample(df2randomize[,1])
  new_df=paste(shuffledCol, fixedCol, sep = '_')
  while(sum(df2randomize_collapsed==new_df)>0){
    new_df[df2randomize_collapsed==new_df]=
      paste(sample(df2randomize[,1], sum(df2randomize_collapsed==new_df)), fixedCol[df2randomize_collapsed==new_df], sep = '_')
  }
  df2return=do.call(rbind, strsplit(new_df, '_'))
  df2return=data.frame(df2return[,1], df2return[,2])
  colnames(df2return)=colnames(df2randomize)
  df2return
}

randomize_DF_with_unique_Pairs_V2<-function(df2randomize=Positive_Set_DrugGene_Pairs){
  df2randomize_collapsed=paste(df2randomize[,1], df2randomize[,2], sep = '_')
  fixedCol=df2randomize[,2]
  shuffledCol=sample(rownames(onTarget$corrMat), nrow(df2randomize))
  new_df=paste(shuffledCol, fixedCol, sep = '_')
  while(sum(df2randomize_collapsed==new_df)>0){
    new_df[df2randomize_collapsed==new_df]=
      paste(sample(rownames(onTarget$corrMat), sum(df2randomize_collapsed==new_df)), fixedCol[df2randomize_collapsed==new_df], sep = '_')
  }
  df2return=do.call(rbind, strsplit(new_df, '_'))
  df2return=data.frame(df2return[,1], df2return[,2])
  colnames(df2return)=colnames(df2randomize)
  df2return
}