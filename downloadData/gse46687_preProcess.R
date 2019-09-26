# 7/11/16 Download and Pre-process GSE46687 Turner's


# Purpose: 
# get GSE46687 into the Khatri lab Dataset format

# Results: 
# One Turner's Syndrome 45Xm sample had a non-normal distribution, so I removed it


#...................#
##### Load Stuff ####
#...................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)

# This script has my cleanup pheno function
source("00_tools/general_GEMfunx.R")

# graphing scripts
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/3_humanPerturbations/gse46687_turner/")
#..................#
##### Download #####
#..................#

gse46687 = getGEOData(gseVector = "GSE46687", formattedNames = "GSE46687_X0")
gse46687 = gse46687$originalData$GSE46687


#...............#
##### Pheno #####
#...............#
gse46687$rawPheno = gse46687$pheno
gse46687$pheno = cleanUpPheno(pheno = gse46687$pheno, removeChar = T)
View(gse46687$pheno)

# Group column 
myGroup = as.character(gse46687$pheno$karyotype)
unique(myGroup)
myGroup[which(myGroup == "46XX")] = "Female XX"
myGroup[which(myGroup == "45Xm")] = "Female X0"
myGroup[which(myGroup == "45Xp")] = "Female X0"
gse46687$pheno$group = as.factor(myGroup)

#..............#
##### Expr #####
#..............#
# Generally between -2 and 16
# Medians a bit jagged but mostly the same
# One Turner patient just has a weird sample
# No obvious batch effect
# Looks fine to me!
boxplot(gse46687$expr, main = "gse46687 - Turner's")


# Should I remove the weird sample?
# Yes, it has a different distribution than the other samples
hist(gse46687$expr[,11], main = "Weird Sample #11, Not Normally Dist")
hist(gse46687$expr[,12], main = "Normal Sample #12, Normally Dist")

# 45Xm sample = Turner's Syndrome with maternal X Chr
gse46687$pheno$title[11]

gse46687 = removeOneSample(gem = gse46687, sampleID = colnames(gse46687$expr)[11])
checkDataObject(gse46687, "Dataset")

#........................#
##### Class and Keys #####
#........................#
# Class isn't important here

#Keys 
# acceptable number of NAs
numNA = sum(is.na(gse46687$keys)) # 12139
numTotal = length(gse46687$keys) # 54675
numNA/numTotal # 23% NA, so that's fine!

# Keys look good!
gse46687$keys[which(!is.na(gse46687$keys)[1:200])]

#......................#
##### Check Labels #####
#......................#


"XIST" %in% gse46687$keys
"RPS4Y1" %in% gse46687$keys
"KDM5D" %in% gse46687$keys

pdf("gse46687_checkLabels.pdf")
quickViolin(gse46687, "XIST", "group")
quickViolin(gse46687, "RPS4Y1", "group")
quickViolin(gse46687, "KDM5D", "group")
dev.off()


#..............#
##### Save #####
#..............#
save(gse46687, file ="gse46687_Turner.RData")
