# 4/20/2018 - Download and Preprocess GSE19442

# Background: 
#  - Chaussabel 
#  - Black South Africans with TB
#  - Using latent TB as "healthy"

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

setwd("0_datasets/1_validation/gse19442/")

#..................#
##### Download #####
#..................#
# Download
gse19442 = getGEOData("GSE19442")
gse19442 = gse19442$originalData$GSE19442

# Clean up pheno
gse19442$rawPheno = gse19442$pheno
gse19442$pheno = cleanUpPheno(gse19442$rawPheno, removeChar = T)
View(gse19442$pheno)


# Sanity check expr
dim(gse19442$expr) # 48k probes, 51 samples
sum(is.na(gse19442$keys)) # 24k NAs
any(is.na(gse19442$expr)) # no NAs in expr
range(gse19442$expr, na.rm = T) # not log transformed
summary(gse19442$expr[,1])

# Clean up expr
gse19442$expr = gse19442$expr + abs(min(gse19442$expr)) + 1.01
min(gse19442$expr) # 1
range(gse19442$expr)
gse19442$expr = log2(gse19442$expr)
range(gse19442$expr) # Now it's the expected range

# Is there batch effect in expr?
# No obvious batch effect, weirdly squished boxplots but otherwise good
boxplot(gse19442$expr, main = "gse19442")

#....................#
##### Sex Labels #####
#....................#
# Check number of males vs females
table(gse19442$pheno$sex) # 25 females, 26 males
table(gse19442$pheno$sex, gse19442$pheno$illness) # 20 latent females, 11 latent males

# Check sex labels
"XIST" %in% gse19442$keys # T
"RPS4Y1" %in% gse19442$keys # T
"KDM5D" %in% gse19442$keys # F

# All match!
imputedSex = imputeSex(gse19442)
table(imputedSex, gse19442$pheno$sex)

# Create Class vector
gse19442$class = createClassVector(gse19442$pheno$sex, "female", gse19442$pheno)

#..........................#
##### Divvy by Illness #####
#..........................#
# Purpose: 
#   - Create a datset object that only has latent TB people
pheno_val = subset(gse19442$pheno, illness == "LATENT TB")
table(pheno_val$sex)
gse19442_val = subsetGEMFromPheno(gse19442, pheno_val)

#..............#
##### Save #####
#..............#

save(gse19442, gse19442_val, file = "gse19442_latentTB.RData")
