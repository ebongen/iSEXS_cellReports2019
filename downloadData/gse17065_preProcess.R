# 5/15/16 gse17065 Moroccan Dataset
# Re-done Juyl 2017 locally

# Purpose:
# Download, preprocess, and integrate age into gse17065

# Methods: 
# Three samples had their age written as a range between 50-55, so I put their ages as 53
# Median of expr was at 0, so I raised all expr values above zero
# ELANE was listed as its old name "ELA2", so I manually renamed it here - 1/14/17

#...................#
##### Load Stuff ####
#...................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/0_discovery/gse17065/")

library(MetaIntegrator)
source("/labs/khatrilab/ebongen/sexDifferences2.0/00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

# Load table with ages
gse17065_age = read.csv("GSE17065_Age.csv", header=T, row.names = 1)

#..................#
##### Download #####
#..................#
# originally done with MetaIntegrator_private
gse17065 = getGEOData("GSE17065","gse17065_Morocco_WB")
gse17065 = gse17065$originalData$GSE17065
checkDataObject(gse17065, "Dataset") # True

#...............#
##### Pheno #####
#...............#
# Clean up pheno
gse17065$rawPheno = gse17065$pheno
gse17065$pheno = cleanUpPheno(gse17065$pheno, removeChar = T)

# Check sex
head(gse17065$pheno$sex)

# Integrate age
age = gse17065_age$ageSimple
names(age) = rownames(gse17065_age)
gse17065$pheno$age = age[rownames(gse17065$pheno)]
head(gse17065$pheno$age) # looks good!



#..............#
##### expr #####
#..............#
# Ranges from -3 to 9, which is okay
# It's Illumina, so we can't easily go from raw data
# This is good enough
# A good deal of variation, but doesn't look like batch effect
pdf("gse17065_expr.pdf", width = 10, height = 10)
boxplot(gse17065$expr, main = "gse17065")
dev.off()

# Make all the expr values be above zero
summary(gse17065$expr) # all the medians are at 0
min(gse17065$expr) #-3.38 is the smallest one

# move above 0, so that there won't be any negative values
gse17065$expr = gse17065$expr - min(gse17065$expr) + 1

#........................#
##### Class and Keys #####
#........................#
# female is cases (1)
# male is controls (0)
gse17065$class = createClassVector(gse17065$pheno$sex, casesAre = "female", pheno = gse17065$pheno)

## Keys
# Looking good!
gse17065$keys[1:100]
length(gse17065$keys) # 48803
length(na.omit(gse17065$keys)) #24850 genes aren't NA


#.........................#
##### Divvy Discovery #####
#.........................#
# Count useful samples
youngPheno = subset(gse17065$pheno, age <=40)
youngPheno = subset(youngPheno, age >=18)
nrow(youngPheno) # 120
table(youngPheno$sex) # 71 female 49 male

hist(youngPheno$age) # relatively even spread between 18 and 40

gse17065_discovery = subsetGEMFromPheno(gse17065, youngPheno)

# Take a look at expr when there's fewer samples
# Eek this is ugly!
boxplot(gse17065_discovery$expr, main = "gse17065 Discovery")

#................................#
##### Doublecheck Sex Labels #####
#................................#
# Purpose: Make sure that males and females cluster approriately by XIST, RPS4Y1, and KDM5D expression

"RPS4Y1" %in% gse17065_discovery$keys # True
"XIST" %in% gse17065_discovery$keys # True
"KDM5D" %in% gse17065_discovery$keys # False

# Clearly separates by RPS4Y1 and XIST 
quickScatter(gse17065_discovery, geneX = "RPS4Y1", geneY = "XIST", columnName = "sex")

# A few mis-labeled samples
quickViolin(gse17065_discovery, "XIST", "sex")

# Find imputed Sex
imputedSex = imputeSex(gse17065_discovery)
table(gse17065_discovery$pheno$sex, imputedSex)

# Remove samples where imputed sample does not match labeled sex
mislabeledSamples = rownames(gse17065_discovery$pheno)[which(gse17065_discovery$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse17065_discovery = removeOneSample(gse17065_discovery, mySample)
}

# Removed mislabeled samples
quickViolin(gse17065_discovery, "XIST", "sex")

#..............#
##### Save #####
#..............#
checkDataObject(gse17065_discovery, "Dataset") # True
checkDataObject(gse17065, "Dataset") # True

save(gse17065_discovery, gse17065, file="gse17156_morocco.RData")
