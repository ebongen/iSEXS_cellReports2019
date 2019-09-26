# 8/25/16 Download and Pre-process gse42331 Klinefelter's
# Updated May 2017 for correct sex labels (sex column is correct, not source_name_ch1)
# Adapted from 7/2016 code

# Purpose: 
# get GSE42331 into the Khatri lab Dataset format


#...................#
##### Load Stuff ####
#...................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/3_humanPerturbations/gse42331_klinefelter/")

#..................#
##### Download #####
#..................#

gse42331 = getGEOData("GSE42331", "GSE42331_XXY")
gse42331 = gse42331$originalData$GSE42331

#...............#
##### Pheno #####
#...............#
# Make pheno cleaner looking
gse42331$rawPheno = gse42331$pheno
gse42331$pheno = cleanUpPheno(pheno = gse42331$rawPheno, removeChar = T)


# Add group column 
# This original way leads a female and Klinefelters sample to swap labels
# because the source_name_ch1 column is wrong
# myGroup = as.character(gse42331$pheno$source_name_ch1)
# unique(myGroup)
# myGroup[which(myGroup == "whole blood (KS)")] = "XXY Male"
# myGroup[which(myGroup == "whole blood (male)")] = "XY Male"
# myGroup[which(myGroup == "whole blood (female)")] = "XX Female"
# gse42331$pheno$group = as.factor(myGroup)

myGroup = paste(as.character(gse42331$pheno$sex), as.character(gse42331$pheno$condition))
myGroup[which(myGroup == "female control")] = "XX Female"
myGroup[which(myGroup == "male control")] = "XY Male"
myGroup[which(myGroup == "male Klinefelter Syndrome")] = "XXY Male"
gse42331$pheno$group = myGroup



#..............#
##### Expr #####
#..............#
# GPL6244	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array
# There's one huge outlier that I set to NA
# Otherwise looks good, range between 2ish and 14ish
# medians mostly line up
# no obvious batch effect


# There is one insane outlier
boxplot(gse42331$expr, main = "gse42331 - Klinefelter's")

# Most samples are within the ranges of 2.2 and 14.5, except one
summary(gse42331$expr)


### Explore outlier

# maximum value is 8170000
max(gse42331$expr)

test = gse42331$expr
sum(test == max(test)) # only one

# which probe has the outlier? 
# probe 32236
for ( i in 1:nrow(test)){
  if(max(test) %in% test[i,]){
    print(i)
  }
}

gse42331$keys[32236] # gene PHF8

# values around 8.3 and 9.0 ish
# probably a misplaced decimal
test[32236,]
summary(test[32236,])

# Setting the outlier value to NA fixes the problem!
test[which(test==max(test))] = NA
boxplot(test)

# remove outlier value from actual expr instead of test object
gse42331$expr[which(gse42331$expr == max(gse42331$expr))] = NA
boxplot(gse42331$expr, main = "gse42331 - outlier removed")

#........................#
##### Class and Keys #####
#........................#
# Class doesn't matter here

# Keys look fine
gse42331$keys[1:200] # lots of NA's in the beginning
sum(is.na(gse42331$keys)) # 8896 NA's the first time I donwloaded, 9079 in May 2017
sum(is.na(gse42331$keys))/length(gse42331$keys) # 27% NA, which is normal

#..........................#
##### Check Sex Labels #####
#..........................#
# Definitely Y chromosome labels are correct
# No XIST, so hard to double check XX labels, but 
# XXY are more similar to XX


"XIST" %in% gse42331$keys # NOpe
"RPS4Y1" %in% gse42331$keys # True
"KDM5D" %in% gse42331$keys # TRUE
"EIF1AX" %in% gse42331$keys # True
"ZFX" %in% gse42331$keys # True

pdf("sexLabelPlots.pdf")
quickViolin(gse42331, "RPS4Y1", "group")
quickViolin(gse42331, "KDM5D", "group")

quickViolin(gse42331, "EIF1AX", "group")
quickViolin(gse42331, "ZFX", "group")
dev.off()

#..............#
##### Save #####
#..............#
save(gse42331, file = "gse42331_XXY.RData")
