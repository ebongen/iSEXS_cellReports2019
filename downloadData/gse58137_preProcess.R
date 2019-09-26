# 8/16/2017 Preprocess gse58137

# Background:
# Australian dataset on immuno-aging
# I'm gonna use it to see how iSEXS changes with age

# Notes:
# Use platform gpl6947 because gpl1 only has 10 males

#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/2_age/postMenopause/gse58137/")

#..................#
##### Download #####
#..................#
gse58137 = getGEOData("GSE58137")
gse58137gpl1 = gse58137$originalData$GSE58137_GPL10558
gse58137gpl6 = gse58137$originalData$GSE58137_GPL6947

# Tidy up Pheno
gse58137gpl1$rawPheno = gse58137gpl1$pheno
gse58137gpl6$rawPheno = gse58137gpl6$pheno
gse58137gpl1$pheno = cleanUpPheno(gse58137gpl1$rawPheno, T)
gse58137gpl6$pheno = cleanUpPheno(gse58137gpl6$rawPheno, T)

View(gse58137gpl1$pheno)

# Fix age
age1 = gse58137gpl1$rawPheno$characteristics_ch1.2
age1 = as.numeric(strsplitVector(age1, ": ", 2))
length(age1) #106
gse58137gpl1$pheno$age = age1
gse58137gpl1$pheno$ageUnits = NULL

age6 = gse58137gpl6$rawPheno$characteristics_ch1.2
age6 = as.numeric(strsplitVector(age6, ": ", 2))
length(age6) # 253
gse58137gpl6$pheno$age = age6
gse58137gpl6$pheno$ageUnits = NULL

# Fix Sex
gse58137gpl1$pheno$sex = ifelse(gse58137gpl1$pheno$sex == "females", yes = "female", no = "male")
table(gse58137gpl1$pheno$sex) # 96 females, 10 males

gse58137gpl6$pheno$sex = ifelse(gse58137gpl6$pheno$sex == "females", yes = "female", no = "male")
table(gse58137gpl6$pheno$sex) # 168 female, 85 male

# Set Class
gse58137gpl1$class = createClassVector(gse58137gpl1$pheno$sex, "female", gse58137gpl1$pheno)
gse58137gpl6$class = createClassVector(gse58137gpl6$pheno$sex, "female", gse58137gpl6$pheno)

# Check keys
# Looks good!
gse58137gpl1$keys[1:10]
gse58137gpl6$keys[1:10]
length(gse58137gpl6$keys) # 21394

# Check expr
# medians of zero, gotta raise it above zero
pdf("gse58137_gpl6_boxplot.pdf")
boxplot(gse58137gpl6$expr, main = "gse58137 gpl6")
dev.off()

# Move entire boxplot above zero
gse58137gpl1$expr  = gse58137gpl1$expr - min(gse58137gpl1$expr, na.rm = T) + 1
gse58137gpl6$expr  = gse58137gpl6$expr - min(gse58137gpl6$expr, na.rm = T) + 1

#.................................#
##### Check Sex Labels - gpl1 #####
#.................................#
# One person labeled as "female" has male-like expression
# I removed them

"XIST" %in% gse58137gpl1$keys # True
"RPS4Y1" %in% gse58137gpl1$keys # True
"KDM5D" %in% gse58137gpl1$keys # False


# One female is actually male
quickViolin(gse58137gpl1, "XIST", "sex")
quickViolin(gse58137gpl1, "RPS4Y1", "sex")
quickScatter(gse58137gpl1, "RPS4Y1", "XIST", "sex")

imputedSex = imputeSex(gse58137gpl1)
sum(imputedSex != gse58137gpl1$pheno$sex) #1 

mismatchedSamples = rownames(gse58137gpl1$pheno)[which(imputedSex != gse58137gpl1$pheno$sex)]
for(i in 1:length(mismatchedSamples)){
  gse58137gpl1 = removeOneSample(gse58137gpl1, mismatchedSamples[i])
}



#.................................#
##### Check Sex Labels - gpl6 #####
#.................................#
# Two labeled "males" are actually female

"XIST" %in% gse58137gpl6$keys # True
"RPS4Y1" %in% gse58137gpl6$keys # True
"KDM5D" %in% gse58137gpl6$keys # False

# Two males are actually female
quickViolin(gse58137gpl6, "XIST", "sex")
quickViolin(gse58137gpl6, "RPS4Y1", "sex")
quickScatter(gse58137gpl6, "RPS4Y1", "XIST", "sex")

# Impute sex and remove mismatched samples
imputedSex = imputeSex(gse58137gpl6)
sum(imputedSex != gse58137gpl6$pheno$sex) #2 

mismatchedSamples = rownames(gse58137gpl6$pheno)[which(imputedSex != gse58137gpl6$pheno$sex)]
for(i in 1:length(mismatchedSamples)){
  gse58137gpl6 = removeOneSample(gse58137gpl6, mismatchedSamples[i])
}

#.............................#
##### Create Young Subset #####
#.............................#
youngPheno = subset(gse58137gpl6$pheno, age <=40)
youngPheno = subset(youngPheno, age >=18)
table(youngPheno$sex) # 75 females, 24 males

# are ages significantly different between males and females?
# Yup, p = 0.01
maleAge = youngPheno$age[youngPheno$sex == "male"]
femaleAge = youngPheno$age[youngPheno$sex == "female"]
t.test(maleAge, femaleAge) # p = 0.01


gse58137gpl6_young = subsetGEMFromPheno(gem = gse58137gpl6, newPheno = youngPheno)

#..............#
##### Save #####
#..............#
save(gse58137gpl6, file = "gse58137_timecourse.RData")
save(gse58137gpl1, file = "gse58371_gpl1_extra.RData")
save(gse58137gpl6_young, file = "gse58137_young.RData")
