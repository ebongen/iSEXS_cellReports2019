# November 19th 2017 
# Is ISEXS score diff in 1st degree SLE relatives? 


# Results: 
# ANA- Healthy Controls have higher XY-ISEXS than everyone else
#   ANA+ HC, first degre relative +/- ANA, and SLE
# SLE has lower autosomal ISEXS than everyone else
#   All the so-called healthy people look the same


# Conclusions: 
# If ANA- HC have higher XY score, they could have more X-escape which
#   according to that paper could reflect lower immune activation. Because
#   T cells reduced their X-escape once they were stimulated. So, it could
#   reflect low-level inflammation in all the other groups. 
# Something higher expressed in males could start being higher expresed
#  in SLE-prone females 
#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(GEOquery)
library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("00_tools/plot_funx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

load("1_metaAnaly/sexMetaObj.RData")
isexs = c(sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$posGeneNames,
          sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$negGeneNames)
isexsXY = c(sexMetaObj$filterResults$xy$posGeneNames, sexMetaObj$filterResults$xy$negGeneNames)
isexsAuto = c(sexMetaObj$filterResults$autosomeOnly$posGeneNames, sexMetaObj$filterResults$autosomeOnly$negGeneNames)

setwd("0_datasets/5_autoImmunity/SLE/blood/gse24706/")

#..................#
##### Download #####
#..................#

gse24706 = getGEOData("GSE24706")
gse24706 = gse24706$originalData$GSE24706

#...............#
##### Pheno #####
#...............#
gse24706$rawPheno = gse24706$pheno
gse24706$pheno = cleanUpPheno(gse24706$rawPheno, T)

# Give ANA status a prettier name
gse24706$pheno$ANAstatus = gse24706$pheno$`ana status:ch1`
gse24706$pheno$`ana status:ch1` = NULL

# Make disease status prettier
gse24706$pheno$diseaseState = gse24706$pheno$`patient status:ch1`
gse24706$pheno$`patient status:ch1` = NULL

# Sex column
imputedSex = imputeSex(gse24706)
gse24706$pheno$sex = imputedSex

"XIST" %in% gse24706$keys# T
"RPS4Y1" %in% gse24706$keys# T
"KDM5D" %in% gse24706$keys #F
quickViolin(gse24706, "XIST", "sex") # Weird wide XIST spread
quickViolin(gse24706, "RPS4Y1", "sex") # VERY clear RPS4Y1 split

# Make a clean group column 
disease = as.character(gse24706$pheno$diseaseState)
disease = strsplitVector(myVector = disease, mySplit = "[(]", whichPart = 2)
disease = strsplitVector(myVector = disease, mySplit = "[)]", whichPart = 1)
group = paste(disease, gse24706$pheno$ANAstatus, sep = "")
gse24706$pheno$group = group




#..............#
##### Expr #####
#..............#
# Illlumina beadchip
# gpl6884

# No obvious batch effect
# Majority of genes are squished between 6.5 and 7
# I don't like the distribution's compressedness
boxplot(gse24706$expr, main = "gse24706")
min(gse24706$expr) # 6.27
summary(gse24706$expr[,1])# 6.595, 6.911

#..............#
##### Keys #####
#..............#
# 50% of probes are NA
# This isn't unheard of, but not the best

length(gse24706$keys)
sum(is.na(gse24706$keys))/length(gse24706$keys) # 54% NA

#..............#
##### Save #####
#..............#

save(gse24706, gse24706_fem, gse24706_fem_clean, file = "gse24706.RData")

#...............#
##### ISEXS #####
#...............#
# Purpose: 
# Does ISEXS change with ANA status? 
# Does ISEXS change with FDR vs HC status? 

phenoFem = subset(gse24706$pheno, sex == "female")
gse24706_fem = subsetGEMFromPheno(gse24706, phenoFem)

# HC ANA- have the highest scores, seemes higher than
#   - HC ANA+
#   - FDR ANA+/-
#   - SLE
pdf("gse24706_xy_firstDegRelatives.pdf", width = 12)
violinPlot(sexMetaObj$filterResults$xy, gse24706_fem, "group") + ggtitle("XY-ISEXS")
dev.off()

# SLE is lower than everyone else
# All the FDR and HC seem about the same
pdf("gse24706_autosome_firstDegRela.pdf", width = 12)
violinPlot(sexMetaObj$filterResults$autosomeOnly, gse24706_fem, "group") + ggtitle("Autosome ISEXS")
dev.off()

xyScore = calculateScore(sexMetaObj$filterResults$xy, datasetObject = gse24706_fem)


## HC ANA- is clearly higher score than everybody 
t.test(xyScore[gse24706_fem$pheno$group== "HC ANA-"], xyScore[gse24706_fem$pheno$group == "HC ANA+"]) # 0.016
t.test(xyScore[gse24706_fem$pheno$group== "HC ANA-"], xyScore[gse24706_fem$pheno$group == "FDR ANA+"]) # 0.0014
t.test(xyScore[gse24706_fem$pheno$group== "HC ANA-"], xyScore[gse24706_fem$pheno$group == "FDR ANA-"]) # 0.053
t.test(xyScore[gse24706_fem$pheno$group== "HC ANA-"], xyScore[gse24706_fem$pheno$group == "SLE ANA+"]) # 0.0039

## SLE has no sig diff between ANA+ HC, or any of the FDR
t.test(xyScore[gse24706_fem$pheno$group== "FDR ANA+"], xyScore[gse24706_fem$pheno$group == "SLE ANA+"]) # 0.16
t.test(xyScore[gse24706_fem$pheno$group== "FDR ANA-"], xyScore[gse24706_fem$pheno$group == "SLE ANA+"]) # 0.6
t.test(xyScore[gse24706_fem$pheno$group== "HC ANA+"], xyScore[gse24706_fem$pheno$group == "SLE ANA+"]) # 0.2

#.....................................#
##### PCA - all Samples all genes #####
#.....................................#

group2 = paste(gse24706$pheno$group, gse24706$pheno$sex)
gse24706$pheno$group2 = group2


isexsAll = extractDataFromGEM(gse24706, c(sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$posGeneNames, 
                               sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$negGeneNames))
isexsAll = t(isexsAll)

isexsAll_pca = prcomp(isexsAll, center = T, scale. = T)

# One ANA+ healthy control is a major outlier
quickPCAplot(gse24706, isexsAll_pca, "group", PCx=1, PCy=2)

#...........................................#
##### PCA - Females only, only XY genes #####
#...........................................#
# Purpose: What does PCA look like?
# Answer: 
#  One of the females is a crazy outlier
#  It's an overlapping cloud with some signal in it
#  But, it's too confusing, best to leave it be 

gse24706_fem_clean = removeOneSample(gse24706_fem, "GSM608733")

femXYgenes = extractDataFromGEM(gse24706_fem_clean, isexsXY)
femXYgenes = t(femXYgenes)

pr.out = prcomp(femXYgenes, center = T, scale. = T)

biplot(pr.out, scale = 0)
quickPCAplot(gse24706_fem_clean, pr.out, "group", PCx=1, PCy=2)

#..................................#
##### Is it immune activation? #####
#..................................#
## Premise:
# If autosomal/XY ISEXS decreases with immune activation
# we'll see the same pattern in infection!

## Results: 
# XY-ISEXS:
#   - Slight decrease in females corresponds with symptom peak?
#   - Major increase in males after 100 hours
# Autosomal-ISEXS:
#   - Separates before and after infection, but mushed together
#     during symptom peak/onset or whatever

## Conclusion: 
#   - Yes, ISEXS XY and Autosomal changes with immune activation
#   - Exactly why is unclear

load("/labs/khatrilab/ebongen/viralChallenge_clean/0_data/DREAM_challenge/fullSynapseData.RData")

dee3 = fullSynapseData$DEE3_H1N1

phenoClean = subset(dee3$pheno, SHEDDING_SC1 == 1)
phenoClean = subset(phenoClean, SYMPTOMATIC_SC2 == 1)


dee3Clean = subsetGEMFromPheno(dee3, phenoClean)
dee3Clean$pheno$isexsXY = calculateScore(sexMetaObj$filterResults$xy, dee3Clean)
dee3Clean$pheno$isexsAuto = calculateScore(sexMetaObj$filterResults$autosomeOnly, dee3Clean)


ggplot(data = dee3Clean$pheno, aes(x=TIMEHOURS, y=isexsXY, col = GENDER)) +
  geom_point(size=2) + geom_smooth() +ggtitle("DEE3 XY ISEXS")


ggplot(data = dee3Clean$pheno, aes(x=TIMEHOURS, y=isexsAuto, col = GENDER)) +
  geom_point(size=2) + geom_smooth() +ggtitle("DEE3 autosomal ISEXS")

#...........................#
##### Cell Proportions? #####
#...........................#
# Questions: 
# Does immunoStates work in this dataset? 
# does xyScore correlate with CD4+ T cell or monocyte levels? 

myMeta = list(originalData = list(gse24706_fem =gse24706_fem))
test = immunoStatesDecov(myMeta)
