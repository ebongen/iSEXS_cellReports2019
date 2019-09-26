# 4/23/2018 - GSE75511 - WB, PBMC, and dried blood spot the same people

# Background:
#   - Cole from UCLA ran this study
#   - Looked at healthy people in Chicago
#   - Mostly White, but some minority people
#   - Compared gene expression from WB, PBMC, and dried blood spot

# Purpose: 
#   - Whether WB or PBMC separates males and females better
#     can tell us what cell types are important
#   - Aka how much granulocytes are important 

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

# Load libraries
library(MetaIntegrator)
library(ggplot2)
library(cowplot)

# Source code
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/4_sortedCells/WBvsPBMC/gse75511/")

#..................#
##### Download #####
#..................#
# download
gse75511 = getGEOData("GSE75511")
gse75511 = gse75511$originalData$GSE75511

# Clean up pheno
gse75511$rawPheno = gse75511$pheno
gse75511$pheno = cleanUpPheno(gse75511$rawPheno, T)

# Add subject ID
subjID = strsplitVector(gse75511$pheno$title, "-", 2)
head(subjID)
subjID = paste("subject", subjID, sep = '')
unique(table(subjID)) # 3 of each
gse75511$pheno$subjectID = subjID

# Add sex label
sex = gse75511$pheno$female
table(sex) # 3 missing
sex[sex == "missing"] = NA
sex = ifelse(sex == "1", yes = "female", no = "male")

sum(is.na(sex)) # 3 NAs
table(sex)
gse75511$pheno$sex = sex

#..............................#
##### Examine Demographics #####
#..............................#
# Purpose:
#   - If enough minorities, could be worth putting in Validation
#
# Results: 
#   - Mostly white people
#   - A dataset from Cole is in Discovery
#   - Hard to decide whether to put PBMC or WB in Validation
#
# Conclusion:
#   - Include in cell type analysis, not in Validation

# Look at subject level demographics
myCol = c("age", "bmi", "alcohol", "smoke", "female", "black", "asian",
          "otherminority", "white", "educattain", "subjectID")

demo = gse75511$pheno[,myCol]
demo = unique(demo)
dim(demo) # 83
nrow(gse75511$pheno)/3 # 83

table(demo$female, demo$age <=40) # 18 males vs 37 females

demoYoung = subset(demo, age <=40)
range(demoYoung$age) # 18-38
table(demoYoung$black) # 5 young black people

apply(demoYoung, 2, table)
#..............#
##### Expr #####
#..............#
# Purpose:
#   - Make sure expr looks good
#
# Results: 
#   - All positive
#   - No NAs
#   - Within expected range for log transformed


range(gse75511$expr) # 6-13

# Check for batch effect
# Very clean, no batch effect
png("boxplot.png", width=  1000)
boxplot(gse75511$expr, main = "GSE75511")
dev.off()

#..............#
##### Keys #####
#..............#
# Nearly 50% NAs
# I don't like it, but it isn't unusual

# Proportion of NAs
sum(is.na(gse75511$keys)) # 16k NAs
sum(is.na(gse75511$keys))/length(gse75511$keys) # 48% NA

# All but KDM5D present
c("XIST", "RPS4Y1", "KDM5D") %in% gse75511$keys

#....................#
##### Sex Labels #####
#....................#
# Purpose: 
#   - Add sex labl 

# Some mislabeling
imputedSex = imputeSex(gse75511)
table(imputedSex, gse75511$pheno$sex) # 7 discordant samples

# No samples with mixed labels
unique(table(gse75511$pheno$sex, gse75511$pheno$subjectID))

# Obvious mislabelings
# NA person is probably female
quickViolin(gse75511, "RPS4Y1", "sex")

# Identify bad samples
badSamples= na.omit(rownames(gse75511$pheno)[imputedSex != gse75511$pheno$sex])
length(badSamples) # 7 samples
badSamples = c(badSamples, rownames(gse75511$pheno)[is.na(gse75511$pheno$sex)])
length(badSamples) # 10 samples (added 3 NA)
badSamples

# Remove bad Samples
for(mySamp in badSamples){
  gse75511 = removeOneSample(gse75511, mySamp)
}


# Nice and clean!
quickViolin(gse75511, "RPS4Y1", "sex")


# Add class label
gse75511$class = createClassVector(gse75511$pheno$sex, "female", gse75511$pheno)

#................................#
##### Divvy by platform Type #####
#................................#
# Purpose: 
#   - Create separate dataset objects for WB, PBMC, and dried blood spot

table(gse75511$pheno$sampletype)

### Whole Blood
phenoWB = subset(gse75511$pheno, sampletype == "PAXgene")
dim(phenoWB)
table(phenoWB$sex, phenoWB$age <=40)
gse75511_wb = subsetGEMFromPheno(gse75511, phenoWB)
gse75511_wb$formattedName = "GSE75511 WB"

### PBMC
phenoPBMC = subset(gse75511$pheno, sampletype == "PBMC")
dim(phenoPBMC)
table(phenoPBMC$sex, phenoPBMC$age <=40)
gse75511_pbmc = subsetGEMFromPheno(gse75511, phenoPBMC)
gse75511_pbmc$formattedName = "GSE75511 PBMC"

### DBS
phenoDBS = subset(gse75511$pheno, sampletype == 'DBS')
dim(phenoDBS)
table(phenoDBS$sex, phenoDBS$age <=40)
gse75511_dbs = subsetGEMFromPheno(gse75511, phenoDBS)
gse75511_dbs$formattedName = "GSE75511 DBS"

#..............#
##### Save #####
#..............#
save(gse75511, gse75511_dbs, gse75511_pbmc, gse75511_wb, file = "gse75511.RData")
