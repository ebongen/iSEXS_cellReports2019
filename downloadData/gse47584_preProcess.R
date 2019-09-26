# 8/25/16 Download and Pre-process gse47584 Klinefelter's
# Adapted from 7/11/2016 code

# Purpose: 
# get GSE42331 into the Khatri lab Dataset format


#...................#
##### Load Stuff ####
#...................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)

# This script has my cleanup pheno function
source("00_tools/general_GEMfunx.R")

# graphing scripts
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/3_humanPerturbations/gse47584_klinefelter/")
#..................#
##### Download #####
#..................#

gse47584 = getGEOData("GSE47584", "GSE47584_XXY")
gse47584 = gse47584$originalData$GSE47584

#...............#
##### Pheno #####
#...............#

gse47584$rawPheno = gse47584$pheno
gse47584$pheno = cleanUpPheno(gse47584$rawPheno, T)
View(gse47584$pheno)

# Create sex column
# No sex annoatated, probably all male
# both klinefelter's and "normal" have high expression of Y-chr gene RPS4Y1
quickViolin(gse47584, "RPS4Y1", "disease_state")

gse47584$pheno$sex = as.factor(rep("male", nrow(gse47584$pheno)))

# Create group column 
myGroup = as.character(gse47584$pheno$disease_state)
unique(myGroup)
myGroup[which(myGroup == "normal")] = "XY Male"
myGroup[which(myGroup == "Klinefelter's Syndrome")] = "XXY Male"
gse47584$pheno$group = as.factor(myGroup)

#..............#
##### Expr #####
#..............#
# GPL14550	Agilent-028004 SurePrint G3 Human GE 8x60K Microarray

# range between 2ish and 16ish
# medians matchup
# no obvious batch effect
boxplot(gse47584$expr, main = "gse47584 - XXY")

#........................#
##### class and keys #####
#........................#
# Class by chromosome number
# 1 = XXY
# 0 = XY
gse47584$class = createClassVector(gse47584$pheno$group, casesAre = "XXY Male", pheno = gse47584$pheno)

# Keys

# An acceptible amount of NAs
numNAprobes = sum(is.na(gse47584$keys)) #16284
numProbes = length(gse47584$keys) # 42405
numNAprobes/numProbes # 38% of probes are NA, this is normal

# probes look good
gse47584$keys[which(!is.na(gse47584$keys))[1:100]]

#......................#
##### Check Labels #####
#......................#

"XIST" %in% gse47584$keys # TRUE
"ZFX" %in% gse47584$keys# TRUE
"RPS4Y1" %in% gse47584$keys # TRUE
"KDM5D" %in% gse47584$keys# TRUE

pdf("gse47584_checkLabels.pdf")
quickViolin(gse47584, "XIST", "group")
quickViolin(gse47584, "RPS4Y1", "group")
quickViolin(gse47584, "ZFX", "group")
dev.off()
#..............#
##### Save #####
#..............#
save(gse47584, file = "gse47584_XXY.RData")
