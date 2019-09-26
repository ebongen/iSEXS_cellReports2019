# Download and Pre-process gse19151
# First done May 22 2016
# Re-done locally with updated MetaIntegrator July 2017
# Contact: Erika Bongen, ebongen@stanford.edu

# Background: 
# gse19151 is a study on thromboembolism with 53 healthy controls between ages [18,40]
# I'm using it as a discoery dataset in my meta-analysis comparing women and men
# Some of the healthy controls were taken from gse17156 (viral challenge RSV, HRV, H3N2)
# So, using this dataset means viral challenge datasets can't be validation

# Dataset attributes: 
# Whole blood
# Affy gpl571
# Age, raceEthnicity, sex

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/0_discovery/gse19151/")
#..................#
##### Download #####
#..................#
gse19151 = getGEOData("GSE19151")
gse19151 = gse19151$originalData$GSE19151


## Expr
# Log2 transformed? Yes!
# Batch effect? No!
# Contains negatives? No!
boxplot(gse19151$expr, main = "gse19151")
min(gse19151$expr) # 2.79

# Expr looks good!

## Pheno  
gse19151$rawPheno = gse19151$pheno
gse19151$pheno = cleanUpPheno(gse19151$pheno, removeChar = T)

# Sex
head(gse19151$pheno$sex) # looks good

# Age
head(gse19151$pheno$age) # looks good

# raceEthnicity
head(gse19151$pheno$raceEthnicity) # looks good

# Reanalysis column 
# Says what GSMs each column maps to in gse17156
# (the reused controls from a viral challenge study)
relations = strsplitVector(gse19151$pheno$relation, ": ", 2)
gse19151$pheno$relation = relations

## Class 
gse19151$class = createClassVector(gse19151$pheno$sex, "female", gse19151$pheno)

## Keys
gse19151$keys[1:100] # looks good

## formattedName
# Let's leave it be
gse19151$formattedName =  "GSE19151_UK_WB"

checkDataObject(gse19151, "Dataset") # True

#............................#
##### Divvy into healthy #####
#............................#
# Create verion of pheno with only healthy people ages [18, 40]
phenoYoung = subset(gse19151$pheno, source_name_ch1 == "blood_healthy control")
phenoYoung = subset(phenoYoung, age <= 40)
table(phenoYoung$age) # 18 to 38
table(phenoYoung$sex) # 26 female, 27 male

# Subset entire Dataset object
gse19151_discovery = subsetGEMFromPheno(gse19151, phenoYoung)
checkDataObject(gse19151_discovery, "Dataset") # True

#..........................#
##### Check Sex Labels #####
#..........................#

"RPS4Y1" %in% gse19151_discovery$keys # True
"XIST" %in% gse19151_discovery$keys # True
"KDM5D" %in% gse19151_discovery$keys # True

# Clearly separates by RPS4Y1 and XIST 
quickScatter(gse19151_discovery, geneX = "RPS4Y1", geneY = "XIST", columnName = "sex")

# Looks like two swapped samples
quickViolin(gse19151_discovery, "XIST", "sex")

# Find imputed Sex
imputedSex = imputeSex(gse19151_discovery)
table(gse19151_discovery$pheno$sex, imputedSex)

# Remove samples where imputed sample does not match labeled sex
mislabeledSamples = rownames(gse19151_discovery$pheno)[which(gse19151_discovery$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse19151_discovery = removeOneSample(gse19151_discovery, mySample)
}

# Removed mislabeled samples
quickViolin(gse19151_discovery, "XIST", "sex")


#..............#
##### Save #####
#..............#
checkDataObject(gse19151_discovery, "Dataset") # True
checkDataObject(gse19151, "Dataset") # True

save(gse19151, gse19151_discovery, file = "gse19151_thromboEmbolism.RData")
