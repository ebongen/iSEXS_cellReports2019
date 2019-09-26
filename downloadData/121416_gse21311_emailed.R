# 12/14/16 Create Dataset from Gibson's emailed data of gse21311

# Background:
# The preprocessed expresion on GEO has scrambled sex lables
# I emailed Gibson, and he doesn't know why
# But, their in house data looks fine, so he emailed that to me

# Purpose: 
# Their in house data is centered on zero, not sure if it's good for what I want to do 
# Put the data he emailed me into a GEM


# Update: January 2018 
# I forgot to check expr for batch effects and normalization
# The data they gave me was very messy and required normalizing
# 

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
library(dplyr)
source("~/Tools/Graphing Scripts/quickPlots.R")
source("00_tools/general_GEMfunx.R")


setwd("0_datasets/1_validation/gse21311/")

#........................#
##### Pre-Processing #####
#........................#
# Load in the data Gibson gave me
expr = read.csv("gibson_redCross_processedExpr.csv", stringsAsFactors = F)
sexLabels = read.delim("gibson_sexLabels.txt")
sampleAttrib = read.csv("mason10_gse21311_sampleAttributes.csv", stringsAsFactors = F)

# Create keys
keys = as.character(expr$TargetID)
names(keys) = rownames(expr)


# Pre-process expr
expr$TargetID = NULL
expr$ProbeID= NULL
expr$Accession = NULL

expr = as.matrix(expr)
class(expr) = "numeric"
rownames(expr)=  names(keys)


# Create the world's smallest pheno
pheno = data.frame(sampleID = colnames(expr), tissue = rep("leukocytes", ncol(expr)), source= "email with Greg Gibson")
rownames(pheno) = colnames(expr)
pheno$sex = sexLabels$sex
pheno$sex = ifelse(pheno$sex == "F", yes = "female", no = "male")

age = as.numeric(sampleAttrib$Age[which(sampleAttrib$Array == rownames(pheno))])
pheno$age = age
pheno$age_note = rep("Age from Table S2 Mason 2010")

# Create class
class = ifelse(pheno$sex == "female", yes = 1, no = 0)
names(class) = rownames(pheno)

gse21311_email <-list(expr = expr, pheno=pheno, keys = keys, class = class, formattedName = "GSE21311_email" )


checkDataObject(gse21311_email, "Dataset") # true


## Looking good! Only one male that seems female-ish
quickViolin(gse21311_email, "XIST", "sex")

#..............................#
##### Quality Control Expr #####
#..............................#
# Added January 2017

boxplot(gse21311_email$expr, main = "gse21311 - email")

exprPretty = normalize.quantiles(gse21311_email$expr)
boxplot(exprPretty, main = "Normalized gse21311")

dim(exprPretty) == dim(gse21311_email$expr)
colnames(exprPretty) = colnames(gse21311_email$expr)
rownames(exprPretty) = rownames(gse21311_email$expr)

# Raise above zero
exprPretty = exprPretty + abs(min(exprPretty, na.rm = T)) + 1
boxplot(exprPretty, main = "Pretty above zero")

gse21311_pretty = gse21311_email
gse21311_pretty$expr = exprPretty
checkDataObject(gse21311_pretty, "Dataset")

# Add annotations
gse21311_pretty$exp_comment = "Data from email normalized by myself"
gse21311_pretty$key_comment = "Data from email"

### Check iSEXS
# Cleaning up expr does not change hump in ladies
load("/labs/khatrilab/ebongen/sexDifferences2.0/1_metaAnaly/sexMetaObj.RData")

myData = gse21311_pretty$pheno
myData$autosomal = calculateScore(sexMetaObj$filterResults$autosomeOnly, gse21311_pretty)

ggplot(myData, aes(x=age, y = autosomal, col = sex)) + geom_smooth() + geom_point(size = 3)



### Check 5 flu genes
# 5 flu genes still change 
load("/labs/khatrilab/ebongen/viralChallenge_clean/1_metaAnaly/fluMeta.RData")
myData = gse21311_pretty$pheno
fluMeta$filterResults$forwardSearch_20fdr$posGeneNames = character()
myData$score = calculateScore(fluMeta$filterResults$forwardSearch_20fdr, gse21311_pretty)

ggplot(myData, aes(x=age, y = score)) + geom_smooth() + geom_point(size = 3)

### Does Deconvolution work? 
# Nope, deconvolution is still terrible
pretyDeconvo = MetaIntegrator::immunoStatesDecov(list(originalData = list(gse21311_pretty = gse21311_pretty)))
hist(pretyDeconvo$immunoStates$gse21311_pretty$Correlation)

save(gse21311_email, file = "gse21311_email_unnormalizedExpr.RData")
gse21311_email = gse21311_pretty
#......................#
##### Divvy by Age #####
#......................#
phenoYoung = subset(pheno, age <=40)
range(phenoYoung$age)

gse21311_email_validation = subsetGEMFromPheno(gse21311_email, phenoYoung)
gse21311_email_validation$formattedName = "GSE21311_validation"

checkDataObject(gse21311_email_validation, "Dataset") # True

phenoOld = subset(gse21311_email$pheno, age >= 50)
range(phenoOld$age)
table(phenoOld$sex)
gse21311_email_over50 = subsetGEMFromPheno(gse21311_email, phenoOld)
gse21311_email_over50$formattedName = "GSE21311_over50"

#..........................#
##### Check Sex Labels #####
#..........................#
# Purpose: Make sure the sex labels match expresion of known genes
# Results: All discovery samples have correct sex labels, according to RPS4Y1 expression


"XIST" %in% gse21311_email_validation$keys # True
"RPS4Y1" %in% gse21311_email_validation$keys # TRUE
"KDM5D" %in% gse21311_email_validation$keys # False

# VAlidation samples
# Separates clearly by sex
# One sample with high XIST and high RPS4Y1
quickViolin(gse21311_email_validation, "XIST", "sex")
quickViolin(gse21311_email_validation, "RPS4Y1", "sex")
quickScatter(gse21311_email_validation, "RPS4Y1", "XIST", "sex")

# Over 50 samples
# Separates clearly by sex
quickViolin(gse21311_email_over50, "XIST", "sex")
quickViolin(gse21311_email_over50, "RPS4Y1", "sex")
quickScatter(gse21311_email_over50, "RPS4Y1", "XIST", "sex")


imputedSex = imputeSex(gse21311_email_validation)
all(imputedSex == gse21311_email_validation$pheno$sex) # True

# No conflicting samples, so don't need to remove any!

#..............#
##### Save #####
#..............#

save(gse21311_email, gse21311_email_validation,gse21311_email_over50, file="gse21311_email.RData")
