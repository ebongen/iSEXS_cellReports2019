# 4/17/16 Preprocessing gse37069

# Background
# WB from Trauma patients at various points of recovery
# This dataset has 32 healthy controls within my age range
# it has repeat samples from GSE11375 and GSE36809

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/1_validation/gse37069/")
#..................#
##### Download #####
#..................#
gse37069 = getGEOData("GSE37069")
gse37069 = gse37069$originalData$GSE37069

#...............#
##### Pheno #####
#...............#
gse37069$rawPheno = gse37069$pheno
gse37069$pheno = cleanUpPheno(gse37069$pheno, removeChar = T)

View(gse37069$pheno)

# Age
head(gse37069$pheno$age) # looks good

# Sex
head(gse37069$pheno$sex) # looks good

# Group 
group = as.character(gse37069$pheno$source_name_ch1)
group[which(group != "Control")] = "traumaPatient"
group[which(group == "Control")] = "healthy"
gse37069$pheno$group = as.factor(group)

#......................................#
##### Class +Keys + Formatted Name #####
#......................................#
## Class
gse37069$class = createClassVector(gse37069$pheno$sex, "female", gse37069$pheno)

## Keys
# look good!
gse37069$keys[1:100]
length(na.omit(gse37069$keys)) # 35986 non-NA keys

## Formatted Name
gse37069$formattedName = "GSE37069_USA_WB"

#..........................#
##### Divvy by healthy #####
#..........................#
phenoDisc = subset(gse37069$pheno, group == "healthy")
phenoDisc = subset(phenoDisc, age <= 40)
table(phenoDisc$age) # between 18 and 40, well spread out
table(phenoDisc$sex) # 14 female, 18 male

gse37069_validation = subsetGEMFromPheno(gse37069, phenoDisc)

#..............#
##### expr #####
#..............#
# Correct range, 
# Log2 transformed? Yup!
# Batch effect? Nope!
# Above zero? Yup!
boxplot(gse37069_validation$expr, main = "gse37069")
min(gse37069_validation$expr) # 2.19, looks good


#..........................#
##### Check Sex Labels #####
#..........................#
# Purpose: Make sure the sex labels match expresion of known genes
# Results: All discovery samples have correct sex labels, according to RPS4Y1 expression


"XIST" %in% gse37069_validation$keys # True
"RPS4Y1" %in% gse37069_validation$keys # TRUE
"KDM5D" %in% gse37069_validation$keys # True

# Separates clearly by sex
# two females look like they might be males
quickViolin(gse37069_validation, "XIST", "sex")
quickViolin(gse37069_validation, "RPS4Y1", "sex")
quickScatter(gse37069_validation, "RPS4Y1", "XIST", "sex")


imputedSex = imputeSex(gse37069_validation)
all(imputedSex == gse37069_validation$pheno$sex) # False

# Remove samples where imputed sample does not match labeled sex
mislabeledSamples = rownames(gse37069_validation$pheno)[which(gse37069_validation$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse37069_validation = removeOneSample(gse37069_validation, mySample)
}

# Looks good!
quickViolin(gse37069_validation, "XIST", "sex")

# Control 19297865 has two samples run on her blood
# GSM909661 is a duplicate sample, so I'll remove it
gse37069_validation = removeOneSample(gse37069_validation, "GSM909661")

#..............#
##### Save #####
#..............#
checkDataObject(gse37069, "Dataset")
checkDataObject(gse37069_validation, "Dataset")

save(gse37069, gse37069_validation, file = "gse37069.RData")
