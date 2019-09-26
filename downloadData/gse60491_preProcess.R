# Download and Preprocess GSE60491 
# (UCLA personality dataset)

# Background
# This will be a discovery dataset! 

## Dataset details
# Mostly undergrad aged people
# PBMCs
# 51% non-white
# Cigarrets per day
# Exercise per day
# Birth control yes/no
# anti-depressent yes/no
# other drug yes/no
# raceEthnicity (White vs notWhite)

#...................#
##### Load Stuff ####
#...................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/0_discovery/gse60491/")

#..................#
##### Download #####
#..................#

gse60491 = getGEOData("GSE60491", formattedNames = "GSE60491_USA_PBMC")
gse60491 = gse60491$originalData$GSE60491
#...............#
##### Pheno #####
#...............#
# Clean up ugly columns in pheno
gse60491$rawPheno = gse60491$pheno
gse60491$pheno = cleanUpPheno(gse60491$pheno, T)

# Age
gse60491$pheno$age # looks good!

# Sex
gse60491$pheno$male
sex = as.character(gse60491$pheno$male)
sex[which(sex == "0")] = "female"
sex[which(sex == "1")] = "male"
gse60491$pheno$sex = as.factor(sex)

# raceEthnicity
raceEthnicity = as.character(gse60491$pheno$caucasian)
raceEthnicity[which(raceEthnicity == "1")] = "White"
raceEthnicity[which(raceEthnicity == "0")] = "NonWhite"
gse60491$pheno$raceEthnicity = as.factor(raceEthnicity)

# Columns that should be numeric
numberCols = colnames(gse60491$pheno)[10:24]
for (myCol in numberCols){
  gse60491$pheno[[myCol]] = as.numeric(as.character(gse60491$pheno[[myCol]]))
  print(myCol)
  print(gse60491$pheno[[myCol]])
}
#..............#
##### Expr #####
#..............#
# Batch effect? Nope!
# right range? 6.3 to 14, Good!
# Negative values? Nope!
# Log2 transformed? Yup!
range(gse60491$expr)
boxplot(gse60491$expr, main = "gse60491")

#......................#
##### Class + Keys #####
#......................#
# Set class vector
gse60491$class = createClassVector(gse60491$pheno$sex, "female", gse60491$pheno)

# Sanity check keys, looks alright!
length(gse60491$keys) # 34568 keys
length(na.omit(gse60491$keys)) # 17896 aren't NA
gse60491$keys[1:100]

#..............................#
##### Divvy into Discovery #####
#..............................#
phenoDisc = subset(gse60491$pheno, age <= 40) # remove old people
phenoDisc = subset(phenoDisc, birthcontrol == 0)

table(phenoDisc$age) # mostly undergrads
table(phenoDisc$sex) # 62 women, 32 men

gse60491_discovery = subsetGEMFromPheno(gse60491, phenoDisc)


#................................#
##### Doublecheck Sex Labels #####
#................................#
# Purpose: Make sure that males and females cluster approriately by XIST, RPS4Y1, and KDM5D expression

"RPS4Y1" %in% gse60491_discovery$keys # True
"XIST" %in% gse60491_discovery$keys # True
"KDM5D" %in% gse60491_discovery$keys # False

# Clearly separates by RPS4Y1 and XIST 
quickScatter(gse60491_discovery, geneX = "RPS4Y1", geneY = "XIST", columnName = "sex")

# One male that looks pretty female, and one female that's caught in between
# My method removes the questionable male, and leaves the questionable female
quickViolin(gse60491_discovery, "XIST", "sex")
quickViolin(gse60491_discovery, "RPS4Y1", "sex")

# Find imputed Sex
imputedSex = imputeSex(gse60491_discovery)
table(gse60491_discovery$pheno$sex, imputedSex)

# Remove samples where imputed sample does not match labeled sex
mislabeledSamples = rownames(gse60491_discovery$pheno)[which(gse60491_discovery$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse60491_discovery = removeOneSample(gse60491_discovery, mySample)
}

# Removed mislabeled samples
quickViolin(gse60491_discovery, "XIST", "sex")
quickViolin(gse60491_discovery, "RPS4Y1", "sex")

#..............#
##### Save #####
#..............#
checkDataObject(gse60491, "Dataset") # True
checkDataObject(gse60491_discovery, "Dataset") # True

save(gse60491, gse60491_discovery, file = "gse60491_UCLA.RData")
