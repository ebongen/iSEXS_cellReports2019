# 5/26/2019 - Preporcess gse65391

# Purpose: 
#   = A reviewer asked how the sex difference gene signatures look
#     in SLE. So, I'm using gse8884


# gse65391 Background:
#   - Massive SLE study with nearly 1000 samples from Virginia Pascual
#   - Turns out 

# iSEXS Results:
#  - Lymphocytes and Monocytes correlate in the right direction, but neighther are significant
#  - SLEDAI doesn't correlate with Auto-iSEXS or XY-iSEXS
#  - Could it be that 18 year olds aren't adult enough? Or pSLE is too different?

#....................#
##### Load Stuff #####
#....................#

setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Source functions
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')
library(cowplot)

setwd("0_datasets/5_autoImmunity/SLE/blood/gse65391/")

# Load Winn's SLE meta-objects
load('/labs/khatrilab/hayneswa/SLE/forErika/testingObjectPub.RData')

gse65391 = testingObjectPub$originalData$GSE65391


# Capitalizes the first character of each word in a a string
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

setwd("0_datasets/5_autoImmunity/SLE/blood/gse65391/")

#..........................#
###### Clean up Pheno ######
#..........................#
# Purpose:
#   = Make pheno easier to look at and standardized with my other GEMs

# Clean up pheno
gse65391$rawPheno = gse65391$pheno
gse65391$pheno = cleanUpPheno(gse65391$rawPheno, T)

# Manually create age column, cuz cleanUpPheno screwed up
age = as.character(gse65391$rawPheno$characteristics_ch1.13)
table(grepl(pattern = "age: ", x = age ))
age = strsplitVector(age, "age: ", 2)
age = as.numeric(age)
range(age) # 6-21
gse65391$pheno$age = age



# Create group column
group = paste(gse65391$pheno$sex, gse65391$pheno$`disease state`)
group = sapply(X = group, simpleCap)
names(group) = NULL
gse65391$pheno$group = group

# Make SLEDAI numeric
gse65391$pheno$sledai = as.numeric(as.character(gse65391$pheno$sledai))

# Make cell type percentages numeric
gse65391$pheno$lymphocyte_percent = as.numeric(as.character(gse65391$pheno$lymphocyte_percent))
gse65391$pheno$monocyte_percent = as.numeric(as.character(gse65391$pheno$monocyte_percent))

#....................#
##### Check expr #####
#....................#
# Purpose:
#   = Make sure expr passes QC
#
# Results:
#   = any NAs? Nope!
#   = Log transformed already
#   = No negative values
#   = Boxplot is a bit messy, but no obvious batch effect

# No NAs
any(is.na(gse65391$expr))

# range 1.9- 15.6
range(gse65391$expr)


# Check the boxplot 
jpeg("gse65391_boxplot.jpeg", width = 1000)
boxplot(gse65391$expr, main = "GSE65391")
dev.off()



#...........................................#
##### Quick Exploration of Healthy kids #####
#...........................................#
# Purpose:
#   = I've never been able to dig into what Auto-iSEXS looks like in childhood
#   = gse65391 actually has a decent age range of healthy kids
#
# Results:
#   = Not enough healthy young males to to draw conclusions there
#   = Enough healthy young females to suggest that Auto-iSEXS may be at a peak
#     at puberty years (~10-12yo) and then decreases (~15-16yo). But not enough
#     females to prove anything
#   = Also, not enough adults to understand where Auto-iSEXS in teenaged females compares
#
# Conclusions:
#   = Enough data to show that something interesting is happening
#   = Not enough to figure out what that is

# Quick check of healthy females
phenoHealthy = subset(gse65391$pheno, `disease state` == "Healthy")
phenoHealthyF = subset(phenoHealthy, sex == 'female')
gse65391_hf = subsetGEMFromPheno(gse65391, phenoHealthyF)
gse65391_hf$pheno$autoISEXS = calculateScore(autoSig, gse65391_hf)

# Interesting! Females seem to be consistently decreasing
# Maybe there's a peak at purberty, that goes down to Adult levels
ggplot(gse65391_hf$pheno, aes(x=age, y=autoISEXS)) + geom_smooth()+
  geom_point(size=2)

### Check out healthy boys and girls
# Boys seem to be increasing? from 9-15 yo 
# Not enough samples to be confident of anything
gse65391_h = subsetGEMFromPheno(gse65391, phenoHealthy)
gse65391_h$pheno$autoISEXS = calculateScore(autoSig, gse65391_h)

ggplot(gse65391_h$pheno, aes(x=age, y=autoISEXS, col=sex)) + geom_smooth()+
  geom_point(size=2) +scale_color_manual(values = c(femColor, malColor))

#...............................#
##### Subset to just Adults #####
#...............................#
# Purpose:
#   = I can't use people < 18yo in the iSEXS paper because child development
#     throws a wrench in the iSEXS score
#   = I need to figure out if the 45 adults in this cohort are enough to show that
#     Auto-iSEXS in SLE tracks with lymphocyte and monocyte levels

# Subset to just adults
phenoAdult = subset(gse65391$pheno, age >=18)
table(phenoAdult$group) # 
gse65391_adult = subsetGEMFromPheno(gse65391, phenoAdult)


# Calculate Auto-iSEXS and XY-iSEXS
gse65391_adult$pheno$autoISEXS = calculateScore(autoSig, gse65391_adult)
gse65391_adult$pheno$xyISEXS = calculateScore(xySig, gse65391_adult)

# Check Autosomal-iSEXS: Sig difference between healthy females and SLE females
# Teh males trend in the expected direction, but there's not a big enough sample size
violinPlot(autoSig, gse65391_adult, "group")

# No difference between SLE and healthy females
# Males are much lower, as expected
violinPlot(xySig, gse65391_adult, "group")


### Check correlations
ggplot(gse65391_adult$pheno, aes(x=lymphocyte_percent, y=autoISEXS, col=sex))+
  geom_smooth(method="lm") + geom_point(size=2) + scale_color_manual(values = c(femColor, malColor))
cor.test(gse65391_adult$pheno$lymphocyte_percent, gse65391_adult$pheno$autoISEXS) # not sig, right direction


ggplot(gse65391_adult$pheno, aes(x=monocyte_percent, y=autoISEXS, col=sex))+
  geom_smooth(method="lm") + geom_point(size=2) + scale_color_manual(values = c(femColor, malColor))
cor.test(gse65391_adult$pheno$monocyte_percent, gse65391_adult$pheno$autoISEXS) # not sig, right direction

ggplot(gse65391_adult$pheno, aes(x=sledai, y=autoISEXS))+
  geom_smooth(method="lm", col='purple') + geom_point(size=2, aes(col=sex)) + scale_color_manual(values = c(femColor, malColor))
cor.test(gse65391_adult$pheno$sledai, gse65391_adult$pheno$autoISEXS) # nearly sig, Neg correlation, I guess we expect


phenoSLEf= subset(gse65391_adult$pheno, group =='Female SLE')
cor.test(phenoSLEf$sledai, phenoSLEf$autoISEXS) # not significant when you look at females alone
cor.test(phenoSLEf$sledai, phenoSLEf$xyISEXS) # not significant when you look at females alone

#..............#
##### Save #####
#..............#

save(gse65391, gse65391_adult, file="gse65391.RData")
