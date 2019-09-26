# 4/20/2018 - Download and Preprocess GSE61821

# Background:
#   - Study from Singapore that looks at mild flu
#     in Singaporeans and severe flu from various
#     medical centers in Asia


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

setwd("0_datasets/1_validation/gse61821/")

#..................#
##### Download #####
#..................#
# Purpose: 
#   - Download gse61821 and clean up pheno


gse61821 = getGEOData("GSE61821")
gse61821 = gse61821$originalData$GSE61821

# Clean up pheno
gse61821$rawPheno = gse61821$pheno
gse61821$pheno = cleanUpPheno(gse61821$pheno, T)


# Add study column 
# Mild flu and other-febrile illness (OFI) came from Singapore's EDEN study
# Severe and moderate flu came from hospitilized pations at various centers in
# south east Asia (e.g. Singapore, Indonesia, Vietnam, Thailand). They're all inthe
# SEAICRN study
study = ifelse(gse61821$pheno$timepoint %in% c("Acute", "Follow Up"), yes = "EDEN", no = "SEAICRN")
table(study)
gse61821$pheno$study = study

# Impute sex
"XIST" %in% gse61821$keys# T
"RPS4Y1" %in% gse61821$keys #T
"KDM5D" %in% gse61821$keys # F

gse61821$pheno$sex = imputeSex(gse61821)
gse61821$pheno$imputedSex = rep(TRUE, nrow(gse61821$pheno))

#.........................................#
##### Preprocess Mild Cohort #####
#.........................................#
phenoEDEN = subset(gse61821$pheno, study == "EDEN")
gse61821_eden = subsetGEMFromPheno(gse61821, phenoEDEN)

# 310 samples
dim(gse61821_eden$pheno)

# Check expr
# No NAs, but relatively large negative values
# May need to raise above zero
range(gse61821_eden$expr) # -9 to 16
hist(apply(gse61821_eden$expr, 2, median)) # Not a median at zero, but close to


# Some messiness, but no obvious batch effect 
png("boxplot_eden.png", width = 1000)
boxplot(gse61821_eden$expr, main = "EDEN")
dev.off()

# Raise expr above zero
gse61821_eden$expr = gse61821_eden$expr + abs(min(gse61821_eden$expr)) + 1.01


# Subset to a healthy pheno
phenoHC = subset(gse61821_eden$pheno, timepoint == "Follow Up")
dim(phenoHC)
table(phenoHC$sex) # 51 females, 104 males
table(phenoHC$sex, phenoHC$virus_type)

# Subset by age
phenoHC = subset(phenoHC, age <=40)
range(phenoHC$age) # 18-37
table(phenoHC$sex) # 34 females, 80 males

# Very similar age distributions
ggplot(phenoHC, aes(x=age, col=sex))+geom_density()

# Create Dataset object 
gse61821_val = subsetGEMFromPheno(gse61821_eden, phenoHC)

# Add class vector
gse61821_val$class = createClassVector(gse61821_val$pheno$sex, "female", gse61821_val$pheno)
checkDataObject(gse61821_val, "Dataset") # True

#....................................#
##### Create SEAICRN Dataset Obj #####
#....................................#
# Purpose: 
#   - Create Dataset object of severe/moderate flu 
#     infected people
#   - before and after infeciton
#
# Results: 
#   - too few post-infection young people to be worth it
#     5 females, 9 males

phenoSEA = subset(gse61821$pheno, study == "SEAICRN")
dim(phenoSEA)
table(phenoSEA$severity)

gse61821_SEAICRN = subsetGEMFromPheno(gse61821, phenoSEA)

# Check expr
range(gse61821_SEAICRN$expr) # -2 to 7

# Messy, but no obvious batch effect
png("boxplot_seaicrn.png", width = 1000)
boxplot(gse61821_SEAICRN$expr, main= "SEAICRN")
dev.off()

# Raise expr above zero
gse61821_SEAICRN$expr = gse61821_SEAICRN$expr + abs(min(gse61821_SEAICRN$expr)) + 1.01

# 25 healthy young people
table(gse61821_SEAICRN$pheno$age <=40, gse61821_SEAICRN$pheno$timepoint)

# Create pheno of healthy young people
# Not worth having
phenoHC = subset(gse61821_SEAICRN$pheno, age <=40)
range(phenoHC$age) # 5-40
phenoHC = subset(phenoHC, age >=18)
phenoHC = subset(phenoHC, timepoint == "day_28")
table(phenoHC$sex) # 5 females vs 9 males

# Do males and females have the same age distributions?
# Completely different age distributions
ggplot(phenoHC, aes(x=age, col = sex)) + geom_density() + ggtitle("Healthy people who had severe flu")

# Create HC dataset obj
gse61821_SEAICRN_val = subsetGEMFromPheno(gse61821, phenoHC)

# Add class vector
gse61821_SEAICRN_val$pheno$class = createClassVector(gse61821_SEAICRN_val$pheno$sex, "female", gse61821_SEAICRN_val$pheno)
checkDataObject(gse61821_SEAICRN_val, "Dataset") # True

#..............#
##### Save #####
#..............#
save(gse61821_eden, gse61821_SEAICRN, gse61821_val, gse61821_SEAICRN_val, file = "gse61821.RData")
