# 4/17/16 Download and pre-process GSE53195
# Re-done locally July 2017

# Background
# GSE53195 is a study of a Brisbane family cohort of parents and children
# many of the kids are twins, some kids are their siblings
# total samples is in the 800s, but 60 people fit within the 18-40 age range I need for discovery

# Interesting variables: 
# good number of kids bfrom 10-18, could get a decent look at puberty
# good number of adults in their 40s, but gets sparce in 50s and after

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/0_discovery/gse53195/")

#..................#
##### Download #####
#..................#

# Download dataset
gse53195 = getGEOData("GSE53195")
gse53195 = gse53195$originalData$GSE53195

#...............#
##### Pheno #####
#...............#

# Clean up pheno
gse53195$rawPheno = gse53195$pheno
gse53195$pheno = cleanUpPheno(gse53195$pheno, removeChar = T)

# SEx
head(gse53195$pheno$sex) # looks good

# Age
head(gse53195$pheno$age) # looks good


#..............................#
##### Divvy into Discovery #####
#..............................#
# New pheno with only samples between 18 and 40
phenoDisc = subset(gse53195$pheno, age <=40)
phenoDisc = subset(phenoDisc, age >=18)

table(phenoDisc$age) # mostly 18 and early 20s
table(phenoDisc$sex) # 34F, 26M

# Create a new dataset object with only 60 discovery samples
gse53195_discovery = subsetGEMFromPheno(gem = gse53195, newPheno = phenoDisc)

#..............#
##### expr #####
#..............#
# No batch effect
# Correct range and looks normalized
boxplot(gse53195_discovery$expr, main = "gse53195")

# No negative values
min(gse53195_discovery$expr, na.rm = T) #min is 5.59
#...............#
##### class #####
#...............#
# Female (F) is case(1)
gse53195_discovery$class = createClassVector(gse53195_discovery$pheno$sex, casesAre = "female", pheno = gse53195_discovery$pheno)

#................................#
##### Keys and formattedName #####
#................................#
# There's an aweful lot of NAs, but it looks alright
gse53195$keys[1:100]
length(gse53195$keys) #47323 probes
length(na.omit(gse53195$keys)) #27561 probes aren't NA
length(unique(na.omit(gse53195$keys))) #18057 unique genes

# They're all Australians of European descent
gse53195_discovery$formattedName = "GSE53195_Euro_WB"

#................................#
##### Doublecheck Sex Labels #####
#................................#
# Purpose: Make sure that males and females cluster approriately by XIST, RPS4Y1, and KDM5D expression
# Perfect agreement with recorded sex and imputed sex

"RPS4Y1" %in% gse53195_discovery$keys # True
"XIST" %in% gse53195_discovery$keys # True
"KDM5D" %in% gse53195_discovery$keys # False

# Clearly separates by RPS4Y1 and XIST 
quickScatter(gse53195_discovery, geneX = "RPS4Y1", geneY = "XIST", columnName = "sex")

#Good separation by both genes
quickViolin(gse53195_discovery, "XIST", "sex")
quickViolin(gse53195_discovery, "RPS4Y1", "sex")

# Find imputed Sex
imputedSex = imputeSex(gse53195_discovery, femGenes = "XIST", malGenes = NULL) # NAs in RPS4Y1, so can't use for classification
table(gse53195_discovery$pheno$sex, imputedSex)

# Perfect classification based on XIST expression, so don't have to remove any samples

#..............#
##### Save #####
#..............#
checkDataObject(gse53195_discovery, "Dataset") # True
checkDataObject(gse53195, "Dataset") # True

save(gse53195, gse53195_discovery, file = "./gse53195.RData")
