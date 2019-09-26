# 5/27/2019 - Preporcess gse39088

# Purpose:
#  = One of the reviewers wanted to see what iSEXS looks like in SLE
#  = gse39088 has a good number of healthy and healthy SLE all female adults

# Results:
#   = expr is centered on 0
#   = But, it still recapitulates what I observed in gse49454


# Conclusion:
#  = We can use gse39088 in sup figure for iSEXS paper

#....................#
##### Load Stuff #####
#....................#


setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Source functions
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')
library(cowplot)


# Load Winn's SLE meta-objects
load('/labs/khatrilab/hayneswa/SLE/forErika/testingObjectPub.RData')

gse39088 = trainingObjectPub$originalData$WB...GSE39088

setwd("0_datasets/5_autoImmunity/SLE/blood/gse39088/")

#..........................#
###### Clean up Pheno ######
#..........................#
# Purpose:
#   = Make pheno easier to look at and standardized with my other GEMs

# Clean up pheno
gse39088$rawPheno = gse39088$pheno
gse39088$pheno = cleanUpPheno(gse39088$rawPheno, T)

# Check group
table(gse39088$pheno$group)

# Create a group column better for plots
group2 = as.character(gse39088$pheno$group)
group2[group2 == 'healthy'] = "Healthy"
group2 = paste("Female\n", group2, sep="")
gse39088$pheno$group2 = group2

#........................#
##### Check out Expr #####
#........................#
# Purpose;
#   = Make sure expr passes my qc
# Results:
#   = No NAs
#   = has negative values


# centered at mother frickin zero
range(gse39088$expr) # -8 to 6.79

gse39088$expr = gse39088$expr + abs(min(gse39088$expr)) +1.01

boxplot(gse39088$expr, main="GSE39088")


# Raising expr has no effect on p-values between the two groups
violinPlot(xySig, gse39088, "group") # reduced in SLE
violinPlot(autoSig, gse39088, "group") # reduced in SLE

#..............#
##### Save #####
#..............#

save(gse39088, file="gse39088.RData")


