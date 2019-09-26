# Preprocess gse18323 -- Malaria Challenge study


# Background: 
# Malaria Challenge study
# 18-45 years old
# 33 males, 15 females

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
library(GEOquery)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/1_validation/gse18323_PBMC/")


#..................#
##### Download #####
#..................#
gse18323 = getGEOData("GSE18323")

gse18323_GPL570 = gse18323$originalData$GSE18323_GPL570
gse18323_gpl571 = gse18323$originalData$GSE18323_GPL571


#.............................#
##### I don't need gpl570 #####
#.............................#
# gpl570 only has 5 days post-challenge
# I'm not interested in it

# gpl570 is only 5 days post-challenge
gse18323_GPL570$rawPheno = gse18323_GPL570$pheno
gse18323_GPL570$pheno = cleanUpPheno(gse18323_GPL570$rawPheno, T)
View(gse18323_GPL570$pheno)

#......................#
##### gpl571 Pheno #####
#......................#
gse18323_gpl571$rawPheno = gse18323_gpl571$pheno
gse18323_gpl571$pheno = cleanUpPheno(gse18323_gpl571$rawPheno, T)
View(gse18323_gpl571$pheno)

# Create timepoint column
myTitle = as.character(gse18323_gpl571$pheno$title)

# Create timeopoint column
gse18323_gpl571$pheno$timepoint = c(as.numeric(strsplitVector(myTitle[1:203], "point T", 2)),
                                    rep(0,length(myTitle)-203))


# Create Subject ID column
subjID1 = strsplitVector(myTitle[1:203], mySplit = " at", 1) # first 203 subjects
subjID2 = strsplitVector(myTitle[204:length(myTitle)], mySplit = ", ", whichPart = 1)
subjID = strsplitVector(c(subjID1, subjID2), " ", whichPart = 2)
gse18323_gpl571$pheno$subjectID = subjID

#...........................#
##### Subset validation #####
#...........................#

# Create the Validation dataset
phenoVal = subset(gse18323_gpl571$pheno, timepoint == 0)
gse18323_validation = subsetGEMFromPheno(gse18323_gpl571, phenoVal)
checkDataObject(gse18323_validation, "Dataset") # True


#....................#
##### Check Expr #####
#....................#

# Looks perfect!
# No batch effect
# good range between 2 and 15
boxplot(gse18323_validation$expr, main = "gse18323 Validation")

# Good, all above zero
min(gse18323_validation$expr) # 2.66

# Check Keys, look good!
print(gse18323_validation$keys[1:20])

#....................#
##### Impute Sex #####
#....................#

"XIST" %in% gse18323_validation$keys # True
"RPS4Y1" %in% gse18323_validation$keys #True
"KDM5D" %in% gse18323_validation$keys # True

imputedSex = imputeSex(gse18323_validation)
gse18323_validation$pheno$imputedSex = imputedSex

# Perfect separation, so don't need to exclude any samples
quickViolin(gse18323_validation, "XIST", "imputedSex")
quickViolin(gse18323_validation, "RPS4Y1", "imputedSex")
quickViolin(gse18323_validation, "KDM5D", "imputedSex")


# Create class vector
gse18323_validation$class = createClassVector(gse18323_validation$pheno$imputedSex, "female", gse18323_validation$pheno)

#..............#
##### Save #####
#..............#

save(gse18323_validation, file = "gse18323_validation.RData")
save(gse18323_GPL570, gse18323_gpl571, gse18323_validation, file = "gse18323_full.RData")
