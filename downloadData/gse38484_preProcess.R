# 4/18/16 Preprocess GSE38484

# Background: 
# Schezophrenia dataset
# Expression values aren't normally distributed, they're strangely comppressed around ~8-10

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/1_validation/gse38484/")

#..................#
##### Download #####
#..................#

# Schizophrenia dataset
gse38484 = getGEOData("GSE38484")
gse38484 = gse38484$originalData$GSE38484

# Clean up pheno
gse38484$rawPheno = gse38484$pheno
gse38484$pheno = cleanUpPheno(gse38484$pheno, T)

# age looks good
gse38484$pheno$age

# Sex looks good
gse38484$pheno$sex

# Create class vector
gse38484$class = createClassVector(gse38484$pheno$sex, "female", gse38484$pheno)

#............................#
##### Validation Dataset #####
#............................#
# Create a Dataset object with only healthy controls age 18-40

myPheno = subset(gse38484$pheno, status == "CONTROL")
myPheno = subset(myPheno, age <=40)
dim(myPheno)

table(myPheno$sex) # 31 female, 25 males


# Create a new Dataset object
gse38484_validation = subsetGEMFromPheno(gse38484, myPheno)
checkDataObject(gse38484_validation, "Dataset")

# Create formatted name
gse38484_validation$formattedName = "GSE38484_Validation"

#....................#
##### Check Expr #####
#....................#
# No obvious batch effect
# Weirdly squished boxplots, but it'll be okay
boxplot(gse38484_validation$expr, main = "gse38484")

#..........................#
##### Check Sex Labels #####
#..........................#
# Purpose: Make sure the sex labels match expresion of known genes
# Results: All discovery samples have correct sex labels, according to RPS4Y1 expression


"XIST" %in% gse38484_validation$keys # True
"RPS4Y1" %in% gse38484_validation$keys # TRUE
"KDM5D" %in% gse38484_validation$keys # False

# Separates clearly by sex
# One male and one female might have swapped
quickViolin(gse38484_validation, "XIST", "sex")
quickViolin(gse38484_validation, "RPS4Y1", "sex")
quickScatter(gse38484_validation, "RPS4Y1", "XIST", "sex")


imputedSex = imputeSex(gse38484_validation)
all(imputedSex == gse38484_validation$pheno$sex) # False

# Remove samples where imputed sample does not match labeled sex
mislabeledSamples = rownames(gse38484_validation$pheno)[which(gse38484_validation$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse38484_validation = removeOneSample(gse38484_validation, mySample)
}

# Looks good!
quickViolin(gse38484_validation, "XIST", "sex")
quickViolin(gse38484_validation, "RPS4Y1", "sex")

#.....................................#
##### Create Healthy only Dataset #####
#.....................................#
# Purpose: healthy timecourse

phenoHealthy = subset(gse38484$pheno, status == "CONTROL")
gse38484_healthy = subsetGEMFromPheno(gse38484, phenoHealthy)


### Check sex labels
"XIST" %in% gse38484_healthy$keys # T
"RPS4Y1" %in% gse38484_healthy$keys # T
"KDM5D" %in% gse38484_healthy$keys # F

# There are some mislabeled samples
quickScatter(gse38484_healthy, "RPS4Y1", "XIST", "sex")

# Impute sex and remove samples
imputedSex = imputeSex(gse38484_healthy)
badSamples = rownames(gse38484_healthy$pheno)[which(imputedSex != gse38484_healthy$pheno$sex)]

# If you check the scatterplot again, this fixes it!
for(mySample in badSamples){
  gse38484_healthy = removeOneSample(gse38484_healthy, mySample)
}

#..............#
##### Save #####
#..............#
checkDataObject(gse38484, "Dataset") # True
checkDataObject(gse38484_validation, "Dataset") # True

save(gse38484, gse38484_validation, gse38484_healthy, file = "gse38484_validation.RData")
