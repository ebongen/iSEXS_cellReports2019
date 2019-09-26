# Preprocess SDY212

# Background: 
#   - Mark Davis HIPC study
#   - Vaccinated people >60 and people <30ish

#....................#
##### Load Stuff #####
#....................#
#
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Source functions
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

# Load packages
library(MetaIntegrator)
library(cowplot)

setwd("0_datasets/6_infection/vaccine/steven/")

# Load Steven's data
stevenData = readRDS(file = '/labs/khatrilab/ebongen/friends/steven/fluCombo.rdat')
length(stevenData) # 22 datasets!

setwd('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/2_age/postMenopause/sdy212/')

#.........................#
##### Check out Pheno #####
#.........................#
# Purpose; 
#   - make sure pheno formatted well
#   - fix Dataset obj formatting

sdy212 = stevenData$SDY212

# Some of the samples are in the wrong order
checkDataObject(sdy212, 'Dataset')

# pheno and expr have the same samples as class
all(rownames(sdy212$pheno) %in% names(sdy212$class)) # true
all(colnames(sdy212$expr) %in% names(sdy212$class)) # True

# Rearrange samples
sdy212$expr = sdy212$expr[,names(sdy212$class)]
sdy212$pheno = sdy212$pheno[names(sdy212$class),]
checkDataObject(sdy212, 'Dataset') # True!

# Add group label
sdy212$pheno$ageGroup = ifelse(sdy212$pheno$age > 40, yes = 'Older', no = 'Younger')
table(sdy212$pheno$ageGroup, sdy212$pheno$sex)
sdy212$pheno$group = paste(sdy212$pheno$ageGroup, sdy212$pheno$sex)

#..............#
##### Expr #####
#..............#
# Purpose: 
#   - Make sure Expr is formatted well
# Results:
#   - Has NAs
#   - range 6-15
#   - No batch effect
#
# Conclusion:
#   - We're good to go!

range(sdy212$expr, na.rm = T)

sum(is.na(sdy212$expr)) # 87 NAs

# No batch effect
png('sdy212_expr.png')
boxplot(sdy212$expr, main = 'SDY212')
dev.off()

#..........................#
##### Check Sex Labels #####
#..........................#
# Purpose: 
#   - Make sure sex labels match gene expr

# Both there!
'XIST' %in% sdy212$keys
'RPS4Y1' %in% sdy212$keys

# One older male and one older female are swapped
quickScatter(sdy212, 'RPS4Y1', 'XIST', 'sex')
quickViolin(sdy212, 'RPS4Y1', 'group')

# Remove swapped samples
imputedSex = imputeSex(sdy212)
badSamples = rownames(sdy212$pheno)[sdy212$pheno$sex != imputedSex]
length(badSamples) # 2 bad samples

for(mySamp in badSamples){
  sdy212 = removeOneSample(sdy212, mySamp)
}

#..............#
##### Save #####
#..............#

save(sdy212, file='sdy212.RData')
