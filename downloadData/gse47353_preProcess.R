# August 2017 Preprocess gse47353

# Background: 
# Study on flu vaccination
# 44 people within my age range
# Affy gpl6244

# Variables:
# timepoints before and after vaccination
# Age, Sex, raceEthnicity

# Notes: 
# 40% of keys are NAs, that's fine

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/0_discovery/gse47353/")
#..................#
##### Download #####
#..................#
gse47353 = getGEOData("GSE47353", formattedNames = "GSE47353_USA_PBMC")
gse47353 = gse47353$originalData$GSE47353

#...............#
##### Pheno #####
#...............#

gse47353$rawPheno = gse47353$pheno
gse47353$pheno = cleanUpPheno(gse47353$pheno, T)
View(gse47353$pheno)

# Age
gse47353$pheno$age # looks good

#Sex 
gse47353$pheno$sex # looks good

# raceEthnicity
# Weirdly misspelled
raceEthnicity = as.character(gse47353$pheno$ethnicty)
unique(raceEthnicity)
raceEthnicity[which(raceEthnicity == "Asain")] = "Asian"
gse47353$pheno$raceEthnicity = as.factor(raceEthnicity)

# Timepoint
timepoint = as.character(gse47353$pheno$sample_collection_time)
timepoint = strsplitVector(timepoint, "\\(day", 2)
timepoint = gsub(pattern = "\\)", replacement = "", x = timepoint)
gse47353$pheno$timepoint = as.numeric(timepoint)

# Vaccine responses
vaccineCol = colnames(gse47353$pheno)[19:28]
for (myCol in vaccineCol){
  gse47353$pheno[[myCol]] = as.numeric(as.character(gse47353$pheno[[myCol]]))
  #print(gse47353$pheno[[myCol]])
}

#......................#
##### Class + Keys #####
#......................#

# Set Class vector
gse47353$class = createClassVector(gse47353$pheno$sex, casesAre = "female", pheno = gse47353$pheno)

# Sanity check keys
gse47353$keys[1:100] # first 100 are all NAs
length(gse47353$keys) # 33297 probes
length(na.omit(gse47353$keys)) # 24218 probes are not NA, number changed in July 2017
length(na.omit(unique(gse47353$keys))) # 22402 unique genes, number changed in July 2017
sum(is.na(gse47353$keys)) / length(gse47353$keys) # 27% of keys are NAs

# Although 27% of keys are NA, it's alright

#................#
##### Subset #####
#................#

# subset useful baseline samples within age range
phenoDisc = subset(gse47353$pheno, timepoint == 0)
phenoDisc = subset(phenoDisc, age <= 40)
dim(phenoDisc) # 44 samples
table(phenoDisc$age) # mostly people in their 20s
table(phenoDisc$sex) # 25 female, 19 male
table(phenoDisc$ethnicty) # 30 White, 10 Asian, 1 Hispanic, 3 AA

gse47353_discovery = subsetGEMFromPheno(gem = gse47353, newPheno = phenoDisc)

#..............#
##### expr #####
#..............#
# Batch effect? No!
# Log2 transformed? Yes!
# All above zero? Yes!

# correct range and mostly normal looking
boxplot(gse47353_discovery$expr, main = "gse47353 Discovery")
min(gse47353$expr) #1.62 

#................................#
##### Doublecheck Sex Labels #####
#................................#
# Purpose: Make sure that males and females cluster approriately by XIST, RPS4Y1, and KDM5D expression

"RPS4Y1" %in% gse47353_discovery$keys # True
"XIST" %in% gse47353_discovery$keys # False
"KDM5D" %in% gse47353_discovery$keys # True

# Clearly separates by RPS4Y1 and KDM5D 
quickScatter(gse47353_discovery, geneX = "RPS4Y1", geneY = "KDM5D", columnName = "sex")

# A few mis-labeled samples
quickViolin(gse47353_discovery, "RPS4Y1", "sex")

# Find imputed Sex
# One female sample that looks male
imputedSex = imputeSex(gse47353_discovery)
table(gse47353_discovery$pheno$sex, imputedSex)

# Remove samples where imputed sample does not match labeled sex
mislabeledSamples = rownames(gse47353_discovery$pheno)[which(gse47353_discovery$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse47353_discovery = removeOneSample(gse47353_discovery, mySample)
}

# Removed mislabeled samples
quickViolin(gse47353_discovery, "RPS4Y1", "sex")


### Double check sex labels in full dataset
quickViolin(gse47353, 'RPS4Y1', 'sex') # 5 mislabeled samples

# Impute sex
imputedSex = imputeSex(gse47353)
table(gse47353$pheno$sex, imputedSex) # 3 'females' are males, 2 'males' are female

# Remove mislabeled samples
mislabeledSamples = rownames(gse47353$pheno)[which(gse47353$pheno$sex != imputedSex)]
length(mislabeledSamples) # 5 samples
for(mySample in mislabeledSamples){
  gse47353 = removeOneSample(gse47353, mySample)
}

# Now perfect separation!
quickViolin(gse47353, 'RPS4Y1', 'sex')


#..............#
##### Save #####
#..............#
checkDataObject(gse47353_discovery, "Dataset") # True
checkDataObject(gse47353, "Dataset") # True

save(gse47353, gse47353_discovery, file = "gse47353.RData")
