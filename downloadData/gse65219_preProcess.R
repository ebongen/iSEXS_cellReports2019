# August 2017 Preprocessing gse65219

# Dataset Background: 
# Study concerning young people vs people in their 90's
# 
# 146 nonagenarians (103 females, 43 males) 
# 30 young controls (19-30 years of age, 21 females, 9 males)



#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/1_validation/gse65219_PBMC/")

#..................#
##### Download #####
#..................#
gse65219 = getGEOData("GSE65219")
gse65219 = gse65219$originalData$GSE65219

#...............#
##### Pheno #####
#...............#
# Pheno looks good!

gse65219$rawPheno = gse65219$pheno
gse65219$pheno = cleanUpPheno(gse65219$rawPheno, T)

View(gse65219$pheno)

# check Sex
head(gse65219$pheno$sex)

# Check age
head(gse65219$pheno$age)

#...............#
##### Class #####
#...............#


gse65219$class = createClassVector(gse65219$pheno$sex, "female", gse65219$pheno)

#...........................................#
##### Divvy into Old and Young datasets #####
#...........................................#

# Create dataset for validation with only young people
phenoYoung = subset(gse65219$pheno, age <=40)
range(phenoYoung$age) # 19 to 30
gse65219_validation = subsetGEMFromPheno(gse65219, phenoYoung)
gse65219_validation$formattedName = "GSE65219_validation"
checkDataObject(gse65219_validation, "Dataset") # True

# Create a dataset of only old people
phenoOld = subset(gse65219$pheno, age > 40)
range(phenoOld$age) # only 90's
gse65219_90s = subsetGEMFromPheno(gse65219, phenoOld)
gse65219_90s$formattedName = "GSE65219_nonagenarian"
checkDataObject(gse65219_90s, "Dataset") # True

#..............#
##### Expr #####
#..............#
# Check for bath effect and weird expression distributions

# Minimum is 6.67, so don't have to raise above zero
min(gse65219$expr)

# No major batch effect in validation dataset
# Just very squashed boxplots
boxplot(gse65219_validation$expr, main = "gse65219 Young")


# What about old people? They look good to me
# No major batch effect
boxplot(gse65219_90s$expr[,1:50], main = "First 50 old people")
boxplot(gse65219_90s$expr[,50:100], main = "Second 50 old people")

#.......................................#
##### Check Sex Labels - Validation #####
#.......................................#
# Sex labels agree completely in Validation/Young peopel set

"XIST" %in% gse65219$keys # True
"RPS4Y1" %in% gse65219$keys # True
"KDM5D" %in% gse65219$keys # True


# Separates clearly by sex
# Looks to be entirely correct
quickViolin(gse65219_validation, "XIST", "sex")
quickViolin(gse65219_validation, "RPS4Y1", "sex")
quickScatter(gse65219_validation, "RPS4Y1", "XIST", "sex")

# Entirely matches, don't need to change anything!
imputedSex = imputeSex(gse65219_validation)
all(imputedSex == gse65219_validation$pheno$sex) # True


#....................................#
##### Check Sex Labels - Elderly #####
#....................................#
# A couple samples are clearly sex labeled wrong, so my code removed them
# Two additional samples are a little weirdly in-between, and my code leaves them in
# I"m okay with this

"XIST" %in% gse65219$keys # True
"RPS4Y1" %in% gse65219$keys # True
"KDM5D" %in% gse65219$keys # True


# Looks like ~2 samples might have swapped labels
# Also one weird one
quickViolin(gse65219_90s, "XIST", "sex")
quickViolin(gse65219_90s, "RPS4Y1", "sex")
quickScatter(gse65219_90s, "RPS4Y1", "XIST", "sex")


imputedSex = imputeSex(gse65219_90s)
all(imputedSex == gse65219_90s$pheno$sex) # False

mislabeledSamples = rownames(gse65219_90s$pheno)[which(gse65219_90s$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse65219_90s = removeOneSample(gse65219_90s, mySample)
}

# Looks good!
# Still a couple samples that are weirdly in-between by RPS3Y1 and KDM5D expression
quickViolin(gse65219_90s, "XIST", "sex")
quickScatter(gse65219_90s, "RPS4Y1", "XIST", "sex")
quickScatter(gse65219_90s, "KDM5D", "XIST", "sex")
quickScatter(gse65219_90s, "RPS4Y1", "KDM5D", "sex")

#..............#
##### Save #####
#..............#

save(gse65219, gse65219_90s, gse65219_validation, file = "gse65219_validation.RData")
