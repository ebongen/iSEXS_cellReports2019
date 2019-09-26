# Preprocess gse26163
# December 2017


#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("00_tools/plot_funx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/4_sortedCells/cd4_Tcells/gse26163/")

#..................#
##### Download #####
#..................#
gse26163 = getGEOData("GSE26163")
gse26163 = gse26163$originalData$GSE26163

#...............#
##### Pheno #####
#...............#
# Tidy pheno
gse26163$rawPheno = gse26163$pheno
gse26163$pheno = cleanUpPheno(gse26163$rawPheno, T)

# Remove redundant columns
gse26163$pheno$`age:ch1` = NULL
gse26163$pheno$`body mass index (kg/m2):ch1` = NULL
gse26163$pheno$`cell type:ch1` = NULL
gse26163$pheno$`day:ch1` = NULL
gse26163$pheno$`person:ch1` = NULL
gse26163$pheno$`Sex:ch1` = NULL

#.............#
##### Expr ####
#.............#
# No obvious batch effect
boxplot(gse26163$expr, main = "gse26163")

# -2.25, 15
# Raise all values above zero
range(gse26163$expr, na.rm = T)
gse26163$expr = gse26163$expr - min(gse26163$expr, na.rm = T) + 1

#....................#
##### Impute Sex #####
#....................#
# Purpose: 
#   Make sure each sample's labeled sex agrees with its 
#   gene expression
#
# Observiation: 
#   three samples have labeled sex that doesn't agree with their gene expression
#      Two CD4+ T cell samples may have been swapped
#      One NK cell sample labled "male" has female-like expression

# Results: 
#   I overwrote the labeled sex with the expression-based sex labels
#   Since I don't have the power to look at age in this cohort anyways
#   It doesn't matter whether the labeled ages are wrong

"XIST" %in% gse26163$keys # T
"RPS4Y1" %in% gse26163$keys # T
"KDM5D" %in% gse26163$keys # F

table(gse26163$pheno$sex) # 31 females, 37 males

# Two females have male-like expression
# One male has female-like expression
imputedSex = imputeSex(gse26163)
table(imputedSex, gse26163$pheno$sex)

quickViolin(gse26163, "RPS4Y1", "sex")
quickViolin(gse26163, "XIST", "sex")
quickScatter(gse26163, "RPS4Y1", "XIST", "sex")

# A male and a female CD4 T cell sample may have swapped
# One NK sample has the wrong labeled sex
View(gse26163$pheno[gse26163$pheno$sex != imputedSex,])

# Remove samples with incorrect sex labels
badSamples = rownames(gse26163$pheno)[gse26163$pheno$sex != imputedSex]
for(mySamp in badSamples){
  gse26163 = removeOneSample(gse26163, mySamp)
}

# Now it separates cleanly
quickViolin(gse26163, "RPS4Y1", "sex")

### Create class vector
gse26163$class = createClassVector(gse26163$pheno$sex, "female", gse26163$pheno)

#.............................#
##### CD4 T cell subset   #####
#.............................#
## Create a Dataset object of CD4 T cell samples
phenoCD4 = subset(gse26163$pheno, cell_type == "CD4+ lymphocytes")
gse26163_CD4 = subsetGEMFromPheno(gse26163, phenoCD4)
gse26163_CD4$formattedName = "GSE26163_CD4"

### Examine how variable the score is in males and females
# Sex chromosome score is very tight and consistent 
# No obvious sex difference in the variability of autosomal score
phenoCD4$xyScore = calculateScore(sexMetaObj$filterResults$xy, gse26163_CD4)
phenoCD4$autoScore = calculateScore(sexMetaObj$filterResults$autosomeOnly, gse26163_CD4)

ggplot(data = phenoCD4, aes(x=person, y = xyScore)) + geom_violin(aes(col = sex)) + geom_jitter(aes(col = sex)) + 
  scale_color_manual(values = c(femColor, malColor)) + ggtitle("CD4 - XY Score")

ggplot(data = phenoCD4, aes(x=person, y = autoScore)) + geom_violin(aes(col = sex)) + geom_jitter(aes(col = sex), width = 0.1) + 
  scale_color_manual(values = c(femColor, malColor)) + ggtitle("CD4 - Autosomal Score")

### Only keep first sample from each person
for(thisPerson in unique(gse26163_CD4$pheno$person)){
  if(table(gse26163_CD4$pheno$person)[thisPerson] > 1){
    myDays = gse26163_CD4$pheno$day[gse26163_CD4$pheno$person == thisPerson]
    myDays = as.numeric(strsplitVector(myDays, "day ", 2))
    minDay = min(myDays)
    minDay = paste("day", minDay)
    
    samplesToRemove = rownames(gse26163_CD4$pheno)[gse26163_CD4$pheno$person == thisPerson & gse26163_CD4$pheno$day != minDay]
    for(mySamp in samplesToRemove){
      gse26163_CD4 = removeOneSample(gse26163_CD4, mySamp)
      }
  }
}

rocPlot(sexMetaObj$filterResults$autosomeOnly, gse26163_CD4) # AUC = 0.73
rocPlot(sexMetaObj$filterResults$xy, gse26163_CD4) # AUC = 1

#.........................#
##### NK Cells Subset #####
#.........................#

### Create an NK cell only Dataset object
phenoNK = subset(gse26163$pheno, cell_type == "NK lymphocytes")
gse26163_NK = subsetGEMFromPheno(gse26163, phenoNK)
gse26163_NK$formattedName = "GSE26163_NK"

### Examine how variable the score is in males and females
# Sex chromosome score is very tight and consistent 
# No obvious sex difference in the variability of autosomal score
phenoNK$xyScore = calculateScore(sexMetaObj$filterResults$xy, gse26163_NK)
phenoNK$autoScore = calculateScore(sexMetaObj$filterResults$autosomeOnly, gse26163_NK)

ggplot(data = phenoNK, aes(x=person, y = xyScore)) + geom_violin(aes(col = sex)) + geom_jitter(aes(col = sex)) + 
  scale_color_manual(values = c(femColor, malColor)) + ggtitle("NK - XY Score")

ggplot(data = phenoNK, aes(x=person, y = autoScore)) + geom_violin(aes(col = sex)) + geom_jitter(aes(col = sex), width = 0.1) + 
  scale_color_manual(values = c(femColor, malColor)) + ggtitle("NK - Autosomal Score")


### Only keep first sample from each person
for(thisPerson in unique(gse26163_NK$pheno$person)){
  if(table(gse26163_NK$pheno$person)[thisPerson] > 1){
    myDays = gse26163_NK$pheno$day[gse26163_NK$pheno$person == thisPerson]
    myDays = as.numeric(strsplitVector(myDays, "day ", 2))
    minDay = min(myDays)
    minDay = paste("day", minDay)
    
    samplesToRemove = rownames(gse26163_NK$pheno)[gse26163_NK$pheno$person == thisPerson & gse26163_NK$pheno$day != minDay]
    for(mySamp in samplesToRemove){
      gse26163_NK = removeOneSample(gse26163_NK, mySamp)
    }
  }
}

# Only 15 people left, only has NK samples from 15 people
length(gse26163_NK$pheno$person) == length(unique(gse26163_NK$pheno$person))

rocPlot(sexMetaObj$filterResults$autosomeOnly, gse26163_NK) # AUC = 0.79
rocPlot(sexMetaObj$filterResults$xy, gse26163_NK) # AUC = 1


#.............#
##### Save ####
#.............#
save(gse26163, gse26163_CD4, gse26163_NK, file = "gse26163_cd4_nk.RData")
