# 4/24/2018 Preprocess Japanese eQTL dataset

# Background:
#   - only 7 males younger than 40
#   - Young people don't separate well by autosomal, but old people do
#   - This is weird because it's opposite of other cohorts we've looked at

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

setwd("0_datasets/2_age/postMenopause/gse53351/")
#..................#
##### Download #####
#..................#
# Get from GEO
gse53351 = getGEOData("GSE53351")
gse53351 = gse53351$originalData$GSE53351

# Clean up pheno
gse53351$rawPheno = gse53351$pheno
gse53351$pheno= cleanUpPheno(gse53351$rawPheno, T)

# Add subject ID
subjID = strsplitVector(gse53351$pheno$title, "\\.", 1)
length(unique(subjID)) == nrow(gse53351$pheno) # No duplicates
gse53351$pheno$subjectID = subjID

#..............#
##### Expr #####
#..............#
# Purpose: 
#   - Make sure expr is:
#        - all above zero
#        - log transformed
# 
# Results:
#   - No NAs
#   - Some negative values
#   - Log transformed

# Contains negative values
range(gse53351$expr) # -7 to 12

# Look for batch effect
# No obvious batch effect
png("boxplot.png", width = 1000)
boxplot(gse53351$expr, main = "GSE53351")
dev.off()


# Raise above zero
gse53351$expr = gse53351$expr + abs(min(gse53351$expr)) + 1.01

#..............#
##### Keys #####
#..............#

# Normal level of NAs
sum(is.na(gse53351$keys)) # 16k NAs
sum(is.na(gse53351$keys))/length(gse53351$keys) # 40% NA, that's normal

# All present
c("XIST", "KDM5D", "RPS4Y1") %in% gse53351$keys

#....................#
##### Sex Labels #####
#....................#
# Purpose:
#   - See if sex labels match gene expression


# Impute sex
imputedSex = imputeSex(gse53351)

# No discordent sex labels
table(imputedSex, gse53351$pheno$sex)

# Beautiful separation
quickViolin(gse53351, "RPS4Y1", "sex")

# Create class vector
gse53351$class = createClassVector(gse53351$pheno$sex, "female", gse53351$pheno)

#..............#
##### Save #####
#..............#

save(gse53351, file = "gse53351.RData")
