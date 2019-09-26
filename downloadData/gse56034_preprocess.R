# December 2017 - Download and preprocess gse56034

# Background: 
# gse56034 is magnetic bead purified CD14+ monocytes
# HC from Boston
# a TON of them!

### Duplicate notes: 
# Contains technical replicates
# Identified samples from same person and averaged expression, even if taken from different times
# My numbers of duplicates and their numbers of duplicates from supplemental methods don't agree
#    They're similar though, and I can't find a better way to identify duplicates
#    All duplicates have consistent sex labels and ages within a few years of each other
#    Most duplicates have the same age

#....................#
##### Load Stuff #####
#....................#


setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("~/Tools/Graphing Scripts/quickPlots.R")
source("00_tools/general_GEMfunx.R")
load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/4_sortedCells/monocytes/gse56034/")

#..................#
##### Download #####
#..................#
gse56034 = getGEOData("GSE56034")
gse56034 = gse56034$originalData$GSE56034

#...............#
##### Pheno #####
#...............#
# Clean up Pheno
gse56034$rawPheno = gse56034$pheno
gse56034$pheno = cleanUpPheno(gse56034$rawPheno, T)

### Create age column
gse56034$pheno$age = as.numeric(gse56034$pheno$`age_(yrs)`)


### Check sex labels for agreement with expr
"XIST" %in% gse56034$keys # False
"RPS4Y1" %in% gse56034$keys # True
"KDM5D" %in% gse56034$keys # True

imputedSex = imputeSex(gse56034)
sum(imputedSex != gse56034$pheno$sex) # All of them agree!

# Since sex passes inspection, we can add class vector
gse56034$class = createClassVector(gse56034$pheno$sex, "female", gse56034$pheno)

### Extract subject IDs
subjectID = strsplitVector(gse56034$pheno$title, "\\.", 1)
length(subjectID) # 485
length(unique(subjectID)) # 449 # about 35 duplicate samples
gse56034$pheno$subjectID = subjectID



#..................................#
##### Remove Duplicate samples #####
#..................................#
# Background: 
#  There are some duplicate samples from the same subject
#  I looked in the paper and these are probably technical replicates
#  that they did to deal with batch effect
#  All subject IDs are only assigned one sex, indciating that they respond to the same person

# Purpose: 
# Identify subjects that have duplicate samples
# Identify the GSMs for each subject
# Quickly look at correlations between each GSM (I'm curious, okay?!)
# Then, average the expression between the GSMs

### Identify samples with duplicate expresion
hasDupe = c()
for(mySamp in gse56034$pheno$subjectID){
  if(sum(gse56034$pheno$subjectID == mySamp) > 1){
    hasDupe = c(hasDupe, mySamp)
  }
}
hasDupe = unique(hasDupe)
hasDupe

View(gse56034$pheno[gse56034$pheno$subjectID %in% hasDupe,c("title", "subjectID", "age", "sex")])

### How well are the technical replicates correlated? 
dupeCor = list()
for (mySamp in hasDupe){
  gsm = rownames(gse56034$pheno)[gse56034$pheno$subjectID == mySamp]
  expr1 = na.omit(gse56034$expr)[,gsm[1]]
  myCor = c()
  for(i in 2:length(gsm)){
    myCor = c(myCor,cor(expr1, na.omit(gse56034$expr)[,gsm[i]] ))
  }
  dupeCor[[mySamp]] = myCor
}

### What about all samples?
allCor = cor(na.omit(gse56034$expr))
dim(allCor)
#View(allCor)

### While the duplicated sample include some with lower correlations than expected
# They are enriched for high correlations
# So, I"ll go ahead and collapse these samples together
sum(allCor > 0.99)/length(allCor) # duplicate samples 40%
sum(unlist(dupeCor) > 0.99)/length(unlist(dupeCor)) # all samples13%

### Collapse samples from the same sample together
gse56034 = averageSubjectSamples(gse56034, "subjectID")
checkDataObject(gse56034, "Dataset") # True

# Make age numeric
age = gse56034$pheno$age
ageNumeric = c()
for (i in 1:length(age)){
  if(grepl(pattern = "; ", x = age[i])){
    myAges =as.numeric(unlist(strsplit(x = age[i], split = "; ")))
    ageNumeric[i] = mean(myAges)
  } else{
    ageNumeric[i] = as.numeric(age[i])
  }
}

gse56034$pheno$age = ageNumeric

## Make class vector again
gse56034$class = createClassVector(gse56034$pheno$sex, "female", gse56034$pheno)



#..............#
##### Expr #####
#..............#
# No obvious batch effect. A few crazy negatibe outlier genes
boxplot(gse56034$expr, main = "gse56034")

# Exploring very low outlier (one probe -10ish)
which(gse56034$expr == min(gse56034$expr, na.rm = T))/nrow(gse56034$expr)
min(gse56034$expr[,30], na.rm = T)
View(gse56034$pheno[which(gse56034$pheno$subjectID == gse56034$pheno$subjectID[30]),])
gse56034$expr[gse56034$expr == min(gse56034$expr, na.rm = T)]
gse56034$keys[21572] # It's in UTY, an ISEXS gene higher expressed in males
# she'a s female so that makes sense

# IN order to avoid raising the entire matrix by 10, I'm setting that value to NA
gse56034$expr[which(gse56034$expr == min(gse56034$expr, na.rm = T))] = NA
min(gse56034$expr, na.rm = T) # New min in -5.411, much better

# Raise entire matrix so that min is 1
gse56034$expr = gse56034$expr - min(gse56034$expr, na.rm = T) + 1


#......................#
##### Divvy by age #####
#......................#
# Purpose: 
# If I divvy by decade, I can find age specific sex differences


# Fudge! Statisitically significant difference in age for males and females!
t.test(gse56034$pheno$age[gse56034$pheno$sex == "female"], gse56034$pheno$age[gse56034$pheno$sex == "male"])

# 18-54
range(gse56034$pheno$age)


table(gse56034$pheno$age <=30, gse56034$pheno$sex) # 186 females, 98 males
table(gse56034$pheno$age <=25, gse56034$pheno$sex) # 130 females, 73 males
table((gse56034$pheno$age > 25 & gse56034$pheno$age <=30), gse56034$pheno$sex) # 56 females, 25 males
table((gse56034$pheno$age > 30 & gse56034$pheno$age <=40), gse56034$pheno$sex) # 30 females, 34 males
table((gse56034$pheno$age > 40 & gse56034$pheno$age <=50), gse56034$pheno$sex) # 30 females, 55 males
table(gse56034$pheno$age >50, gse56034$pheno$sex) # 5 females, 11 males


### 18-25 year olds
pheno25 = subset(gse56034$pheno, age <= 25)
gse56034_25 = subsetGEMFromPheno(gse56034, pheno25)
t.test(gse56034_25$pheno$age[gse56034_25$pheno$sex == "female"], gse56034_25$pheno$age[gse56034_25$pheno$sex == "male"]) # p = 0.09
gse56034_25$formattedName = "GSE56034 18-25yo"
checkDataObject(gse56034_25, "Dataset")
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56034_25) # 0.58

### 26-30 year olds
pheno30 = subset(gse56034$pheno, age >25 & age <=30)
range(pheno30$age) # 26-30
t.test(pheno30$age[pheno30$sex == "female"], pheno30$age[pheno30$sex == "male"]) # p = 0.6
gse56034_30 = subsetGEMFromPheno(gse56034, pheno30)
gse56034_30$formattedName = "GSE56034 26-30yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56034_30) # 0.57

### 31-40yo
pheno40 = subset(gse56034$pheno, age <=40 & age > 30)
range(pheno40$age)
t.test(pheno40$age[pheno40$sex == "female"], pheno40$age[pheno40$sex == "male"]) # p = 0.59
hist(pheno40$age[pheno40$sex == "female"]) # Mostly close to 30 and slowly tapers down
hist(pheno40$age[pheno40$sex == "male"]) # More evenly spread, not all that different really
gse56034_40 = subsetGEMFromPheno(gse56034, pheno40)
gse56034_40$formattedName = "GSE56034 31-40yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56034_40) # 0.69

### Over 40
pheno54 = subset(gse56034$pheno, age > 40)
range(pheno54$age) # 41-54
t.test(pheno54$age[pheno54$sex == "female"], pheno54$age[pheno54$sex == "male"]) # p = 0.54
gse56034_54 = subsetGEMFromPheno(gse56034, pheno54)
gse56034_54$formattedName = "GSE56034 41-54yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56034_54) # 0.6

## subset into young
phenoYoung = subset(gse56034$pheno, age <=40)
range(phenoYoung$age) # 18-40, good!
t.test(phenoYoung$age[phenoYoung$sex == "female"], phenoYoung$age[phenoYoung$sex == "male"]) # p = 0.07
table(phenoYoung$sex) # 216 females, 132 males
gse56034_young =subsetGEMFromPheno(gse56034, phenoYoung)
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56034_young) # 0.6

#..............#
##### Save #####
#..............#

save(gse56034, gse56034_young, file = "gse56034_full.RData")
save(gse56034_25, gse56034_30, gse56034_40, gse56034_54, file = "gse56034_ageBuckets.RData")
