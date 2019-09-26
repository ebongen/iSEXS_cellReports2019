# November 28 2017 - Download and preprocess gse65133

# Purpose: 
# this cohort has 10 males and 10 females with immune cell
# frequencies measured via flow

#....................#
##### Load Stuff #####
#....................#


setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")


# Load stuff
library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")


setwd("0_datasets/4_sortedCells/cellProp/gse65133/")

# Load cell frequencies from Franecso (preprocessed from $pheno)
cellProp = read.delim("PBMCs-Fig3a-Flow-Cytometry.txt")
rownames(cellProp) = cellProp$ID

#..................#
##### Download #####
#..................#


gse65133 = getGEOData("GSE65133")
gse65133 = gse65133$originalData$GSE65133

#...............#
##### Pheno #####
#...............#
gse65133$rawPheno = gse65133$pheno
gse65133$pheno = cleanUpPheno(gse65133$rawPheno, T)

# set sample IDs
sampleID = strsplitVector(gse65133$pheno$title, ", ", 2)
gse65133$pheno$sampleID = sampleID

# Add in cell proportions
cellProp = cellProp[gse65133$pheno$sampleID,]
all(gse65133$pheno$sampleID == rownames(cellProp))
gse65133$pheno = as.data.frame(cbind(gse65133$pheno, cellProp))


# Impute sex
"XIST" %in% gse65133$keys # T
"RPS4Y1" %in% gse65133$keys # T
"KDM5D" %in% gse65133$keys # F


imputedSex = imputeSex(gse65133)
table(imputedSex)  # 10 males and 10 females
gse65133$pheno$sex = imputedSex
gse65133$pheno$sexImputed = rep(TRUE, nrow(gse65133$pheno))


#..............#
##### Expr #####
#..............#

# No obvious batch effect
#kindof squished but I've seen that before for ILlumina
boxplot(gse65133$expr, main = "gse65133")

#........................#
##### Class and Keys #####
#........................#
# Create class vector
gse65133$class = createClassVector(gse65133$pheno$sex, "female", gse65133$pheno)

#..............#
##### Save #####
#..............#

save(gse65133, file = "gse65133.RData")
