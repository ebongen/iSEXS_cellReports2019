# 5/28/2018 - Download and preprocess SLE cohorts for sex differences

# Background:
#   - I have 3 genes that are HC Males > HC Females
#   - But, They also seem to be higher in aSLE, but not pSLE
#   - I need large datasets of females only of aSLE and pSLE

#....................#
##### Load Stuff #####
#....................#

setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Source functions
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

# Load Meta object
load("1_metaAnaly/sexMetaObj.RData")

library(cowplot)

# 3 genes I'm interested in
forestPlot(valMetaObj, "CD99") # consistent down (expected cuz PAR1)
forestPlot(valMetaObj, "ZBED1") # consistent down (expected cuz PAR1)
forestPlot(valMetaObj, "CYBB") # Garbage, pity, this one's interesting!

setwd("0_datasets/5_autoImmunity/SLE/")

#.......................#
##### Winn's Object #####
#.......................#
# Purpose: 
#   - Winn gave me some of his SLE Dataset objects
#   - Identify useful Dataset objects from within 


#..................................#
##### Winn - trainingObjectPub #####
#..................................#
# Purpose: 
#   - Identify useful aSLE and pSLE cohorts within trainingObjectPub
#
# Adult SLE:
#   - gse50635:
#        - Niewold, Mayo Clinic
#        - Whole blood
#        - All adult
#        - All female
#        - 16 H, 33 SLE
#   - gse39088
#        - Lauwerys, Brussels 
#        - Whole blood
#        - All adult
#        - All female
#        - 34 H, 26 SLE
#   - gse17755
#        - Nishimoto sensei from Wakayama
#        - PBMC
#        - All adult
#        - All female
#        - 22 H, 22 SLE
#        - DO NOT INCLUDE: early work showed male vs female expression was strange
#   - gse22098 - Adult
#        - WB
#        - Chaussabel, Baylor
#        - All adult
#        - Contains aSLE and pSLE
#        - Females: 20 aSLE, 31 H, 12 pSLE
#        - Males: 8 aSLE, 11 H, 1 pSLE
# Pediatric SLE
#    - gse8650
#       - PBMC
#       - Chaussabel, Baylor
#       - All <18
#       - 3:1, female:male
#       - SLE: 33 females, 5 males
#       - HC: 15 females, 6 males
#
#    - gse22098 - Pediatric
#         - WB
#         - Chaussabel, Baylor
#         - But, how can I guarentee that there's no repeats?
#     - gse65391 - Pediatric
#         - Non-Chaussabel, Baylor
# Rejects: 
#   - gse11909
#       - I'm worried about duplicate samples across Chaussabel datasets
#       - So, I'll only use gse8650 from Chaussabel
load('/labs/khatrilab/hayneswa/SLE/forErika/trainingObjectPub.RData')

# gse50635
names(trainingObjectPub$originalData)


### gse22098 - Adults
# Sex labels not given, so I gotta check them
'XIST' %in% trainingObjectPub$originalData$WB...GSE22098..Adults$keys # T
'RPS4Y1' %in% trainingObjectPub$originalData$WB...GSE22098..Adults$keys # T
'KDM5D' %in% trainingObjectPub$originalData$WB...GSE22098..Adults$keys # F

# Log transform gse22098 
range(trainingObjectPub$originalData$WB...GSE22098..Adults$expr) # -38 to 54k! Not log transformed!
gse22098_adult = trainingObjectPub$originalData$WB...GSE22098..Adults
gse22098_adult$expr = gse22098_adult$expr + abs(min(gse22098_adult$expr)) + 1.01
gse22098_adult$expr = log2(gse22098_adult$expr)

# Looks like there's some batch effect that needs dealing with
boxplot(gse22098_adult$expr, main= 'gse22098')

# Quantile normalize
test = preprocessCore::normalize.quantiles(gse22098_adult$expr)
boxplot(test, main = 'normalized expr')
rownames(test) = rownames(gse22098_adult$expr)
colnames(test) = colnames(gse22098_adult$expr)
gse22098_adult$expr = test
checkDataObject(gse22098_adult, "Dataset") # T

### It's gonna be hard to differentiate aSLE males and females
# I'll let imputeSex make the cutoff and kick out the aSLE patients with 
# midling expression of XIST and RPS4Y1
quickViolin(gse22098_adult, 'XIST', 'group') # Hard to separate people
quickViolin(gse22098_adult, 'RPS4Y1', 'group') # obvious males in H and pSLE, aSLE is weird
quickScatter(gse22098_adult, 'RPS4Y1', 'XIST', 'group')

gse22098_adult$pheno$imputedSex = imputeSex(gse22098_adult)
quickScatter(gse22098_adult, 'RPS4Y1', 'XIST', 'imputedSex')
table(gse22098_adult$pheno$imputedSex, gse22098_adult$pheno$group)


#.................................#
##### Winn - testingObjectPub #####
#.................................#
# GSE65391 - Useful non-Chaussabel, but also from Baylor pSLE cohort
#   - GPL10558
# Lots of annotation:
#   - gse49454
#   - gse88884
#
# PBMCs:
#   - GSE50772
#        - Genentech study confirming IFN signature
#        - University of Michigan patients
#        - All adults, m
#   - GSE61635


# Load Winn's object
load('/labs/khatrilab/hayneswa/SLE/forErika/testingObjectPub.RData')

names(testingObjectPub$originalData)

### gse65391
View(testingObjectPub$originalData$GSE65391$pheno[1:5,])
pheno = cleanUpPheno(testingObjectPub$originalData$GSE65391$pheno, T)
age = as.character(testingObjectPub$originalData$GSE65391$pheno$characteristics_ch1.13)
age = strsplitVector(age, ': ', 2)
head(age)
age = as.numeric(age)
hist(age)

for(myDataset in testingObjectPub$originalData){
  if(grepl(pattern = "PBMC", x = myDataset$formattedName)){
    print(myDataset$formattedName)
    print(table(myDataset$class))
  }
  
}

### gse50772
# Probably enough numbers of males to comare HC males vs HC females vs SLE females
'RPS4Y1' %in% testingObjectPub$originalData$GSE50772$keys
quickViolin(testingObjectPub$originalData$GSE50772, 'RPS4Y1', 'disease status')
quickViolin(testingObjectPub$originalData$GSE50772, 'CD99', 'disease status') # No diff
imputedSex = imputeSex(testingObjectPub$originalData$GSE50772)
group = paste(imputedSex, testingObjectPub$originalData$GSE50772$pheno$`disease status`)
testingObjectPub$originalData$GSE50772$pheno$group = group

# But, a bit unclear
quickViolin(testingObjectPub$originalData$GSE50772, 'CD99', 'group') # subset of SLE with high CD99

### gse61635
# Lilly in Indy, also University of Miami?
View(testingObjectPub$originalData$GSE61635$pheno)

#.............................#
##### Datasets to Analyze #####
#.............................#

# TrainingPub object - Adults
#   - gse50635: - WB
#   - gse39088 - WB
#   - gse22098 - Adult - WB
#
# Training Pub object - Peds
#    - gse8650
#    - gse22098 - Pediatric
#
# Testing Pub object - Adults
#   - gse49454 - WB
#   - gse88884 - WB
#   - GSE50772 (if too much trouble, ditch)
#   - GSE61635 (if too much trouble, ditch)
#
# Testing Pub object - Peds
#   - GSE65391

sleDatasets = list()
sleDatasets_males = list()


#..................#
##### gse50635 - WB #####
#..................#
#   - gse50635:
#        - Niewold, Mayo Clinic
#        - Whole blood
#        - All adult
#        - All female
#        - 16 H, 33 SLE
#        - All in 40s and 50s

# Get dataset object
sleDatasets$gse50635 = trainingObjectPub$originalData$WB...GSE50635

# Check out pheno
sleDatasets$gse50635$rawPheno = sleDatasets$gse50635$pheno
sleDatasets$gse50635$pheno = cleanUpPheno(sleDatasets$gse50635$rawPheno, T)
View(sleDatasets$gse50635$pheno)

# Check expr
#   - No NAs
#   - All above zero
#   - No obvious batch effect
#   - some wibbliness, but not too bad
range(sleDatasets$gse50635$expr) #1-14
boxplot(sleDatasets$gse50635$expr, main ='gse50635')

# Check sex lables
#   - All labeled as female
#   - All have low expression of KDM5D and RPS4Y1
table(sleDatasets$gse50635$pheno$sex) # all female
'RPS4Y1' %in% sleDatasets$gse50635$keys # T
'KDM5D' %in% sleDatasets$gse50635$keys # T
quickScatter(sleDatasets$gse50635, 'RPS4Y1', 'KDM5D', 'subject_type')

# Tidy group labels
sleDatasets$gse50635$pheno$group[sleDatasets$gse50635$pheno$group == 'H'] = "healthy"

# Check age
#   - A good chunk are >50
hist(sleDatasets$gse50635$pheno$age)

# Check class
table(sleDatasets$gse50635$class, sleDatasets$gse50635$pheno$group)

#..................#
##### gse39088 - WB #####
#..................#
#   - gse39088
#        - Lauwerys, Brussels 
#        - Whole blood
#        - All adult
#        - All female
#        - 34 H, 26 SLE

# get dataset
sleDatasets$gse39088 = trainingObjectPub$originalData$WB...GSE39088

# Clean up pheno
sleDatasets$gse39088$rawPheno = sleDatasets$gse39088$pheno
sleDatasets$gse39088$pheno = cleanUpPheno(sleDatasets$gse39088$rawPheno, T)

# Expr 
#   - Contains negative values
#   - No NAs
#   - Centered around zero
#   - So, I raised it above zero
range(sleDatasets$gse39088$expr) # -8 to 6
boxplot(sleDatasets$gse39088$expr, main = 'gse39088')
sleDatasets$gse39088$expr = sleDatasets$gse39088$expr + abs(min(sleDatasets$gse39088$expr)) + 1.01

# Passes
checkDataObject(sleDatasets$gse39088, "Dataset")

# Check sex labels
#   - One outlier by XIST measurement, but 
#     no males, so hard to define
#   - I'll leave it in
table(sleDatasets$gse39088$pheno$sex) # all female
'XIST' %in% sleDatasets$gse39088$keys
'RPS4Y1' %in% sleDatasets$gse39088$keys
quickScatter(sleDatasets$gse39088, 'RPS4Y1', 'XIST', 'group')


# Tidy group labels
table(sleDatasets$gse39088$pheno$group)
sleDatasets$gse39088$pheno$group[sleDatasets$gse39088$pheno$group == 'healthy'] = "healthy"

# Check age
#   - Moslty younger than 40
hist(sleDatasets$gse39088$pheno$age)
range(sleDatasets$gse39088$pheno$age) # 19 to 50

# Class
table(sleDatasets$gse39088$class, sleDatasets$gse39088$pheno$group)
#........................#
##### gse22098 Adult #####
#........................#
#   - gse22098 - Adult
#        - WB
#        - Chaussabel, Baylor
#        - All adult
#        - Contains aSLE and pSLE
#        - Females: 20 aSLE, 31 H, 12 pSLE
#        - Males: 8 aSLE, 11 H, 1 pSLE


sleDatasets$gse22098_adult = trainingObjectPub$originalData$WB...GSE22098..Adults

# Pheno
sleDatasets$gse22098_adult$rawPheno = sleDatasets$gse22098_adult$pheno
sleDatasets$gse22098_adult$pheno = cleanUpPheno(sleDatasets$gse22098_adult$rawPheno, T)


# Expr 
#   - Not log2 transformed, so I raised it above zero and did that
#   - There's definitely batch effect
range(sleDatasets$gse22098_adult$expr) # -38 to 54k
sleDatasets$gse22098_adult$expr = sleDatasets$gse22098_adult$expr + abs(min(sleDatasets$gse22098_adult$expr)) + 1.01
sleDatasets$gse22098_adult$expr = log2(sleDatasets$gse22098_adult$expr)
range(sleDatasets$gse22098_adult$expr) # 0.14 - 15
boxplot(sleDatasets$gse22098_adult$expr, main = 'gse22098')

### Expr batch effect
#    - Batch effect seems to correspond to group
illnessColors = sleDatasets$gse22098_adult$pheno$illness
table(illnessColors)
illnessColors[illnessColors == 'ASLE'] = 'grey'
illnessColors[illnessColors == 'pSLE'] = 'pink'
illnessColors[illnessColors == 'PSLE'] = 'pink'
illnessColors[illnessColors == 'Still'] = 'blue'
illnessColors[illnessColors == 'Strep and Staph'] = 'red'
boxplot(sleDatasets$gse22098_adult$expr, col = illnessColors, main = 'gse22098')


### Quantile normalize expr
test = preprocessCore::normalize.quantiles(sleDatasets$gse22098_adult$expr)
boxplot(test, main = 'normalized expr')
rownames(test) = rownames(sleDatasets$gse22098_adult$expr)
colnames(test) = colnames(sleDatasets$gse22098_adult$expr)
sleDatasets$gse22098_adult$expr = test
checkDataObject(sleDatasets$gse22098_adult, "Dataset") # T

### Sex Labels
# There are some males
table(sleDatasets$gse22098_adult$pheno$sex, sleDatasets$gse22098_adult$pheno$group)
'XIST' %in% sleDatasets$gse22098_adult$keys # T
'RPS4Y1' %in% sleDatasets$gse22098_adult$keys # T
'KDM5D' %in% sleDatasets$gse22098_adult$keys # F

# adult SLE has some weird looking people
quickScatter(sleDatasets$gse22098_adult, 'RPS4Y1', 'XIST', 'group')

# At leaste one 'male' has 'female-like' expression
quickScatter(sleDatasets$gse22098_adult, 'RPS4Y1', 'XIST', 'sex')
quickScatter(gse22098_adult, 'RPS4Y1', 'XIST','imputedSex')

# Remove samples whose gene expression doesn't match label
imputedSex = imputeSex(sleDatasets$gse22098_adult)
badSamples = rownames(sleDatasets$gse22098_adult$pheno)[imputedSex != sleDatasets$gse22098_adult$pheno$sex]
for(mySample in badSamples){
  sleDatasets$gse22098_adult = removeOneSample(sleDatasets$gse22098_adult, mySample)
}


# Remove pediatric SLE patients
phenoClean = subset(sleDatasets$gse22098_adult$pheno, group != "PSLE" )
sleDatasets$gse22098_adult = subsetGEMFromPheno(sleDatasets$gse22098_adult, phenoClean)

# Check age 
hist(sleDatasets$gse22098_adult$pheno$age) # half below 40

# Check Class
table(sleDatasets$gse22098_adult$class, sleDatasets$gse22098_adult$pheno$group)

# Group 
table(sleDatasets$gse22098_adult$pheno$group)
sleDatasets$gse22098_adult$pheno$group[sleDatasets$gse22098_adult$pheno$group == 'H'] = 'healthy'

# Save dataset obj that has males
# Males - 4 SLE vs 4 HC
table(sleDatasets$gse22098_adult$pheno$sex, sleDatasets$gse22098_adult$pheno$group)
sleDatasets_males$gse22098_adults = sleDatasets$gse22098_adult

# Remove Males
phenoClean = subset(sleDatasets$gse22098_adult$pheno, sex == 'female')
sleDatasets$gse22098_adult = subsetGEMFromPheno(sleDatasets$gse22098_adult, phenoClean)


#..................#
##### gse49454 - WB #####
#..................#
# Background: 
#   - Chaussabel, Baylor
#   - Longitudinal study of different IFN-assoc. gene modules
#   - Related modules to disease severity

# Load my version of gse49454
# expr and sex labels are already checked
load('blood/gse49454/gse49454.RData')

### Subject IDs
# There are 62 patients, according to their paper
patientID = strsplitVector(gse49454$pheno$title, 'V', 1)
length(unique(patientID)) # 63!
patientID
unique(table(patientID)) # max 6 visits, 20 controls?
unique(gse49454$pheno$group[patientID == "C3"]) # C3 refers to healthy controls
patientID[patientID == 'C3'] = as.character(gse49454$pheno$title)[patientID =='C3']
length(unique(patientID)) # 82 --> 62 patients + 20 healthy controls
gse49454$pheno$subjectID = patientID

### Visit IDs
visitID = strsplitVector(gse49454$pheno$title, 'V', 2)
table(visitID)
visitID[grepl(pattern = 'S', x = visitID)] = '1'
gse49454$pheno$visitID = as.numeric(visitID)

### Group 
table(gse49454$pheno$group)
group = as.character(gse49454$pheno$group)
group[group == 'Healthy control of SLE'] = 'healthy'
gse49454$pheno$group = group

### Class
table(gse49454$class) # all zeros
table(gse49454$pheno$group)
gse49454$class = createClassVector(gse49454$pheno$group, 'SLE', gse49454$pheno)

### Remove duplicate subjects
samplesToUse =c()
for(mySubject in unique(gse49454$pheno$subjectID)){
  if(sum(gse49454$pheno$subjectID == mySubject) == 1){
    samplesToUse = c(samplesToUse, rownames(gse49454$pheno)[gse49454$pheno$subjectID == mySubject])
  } else{
    miniPheno = subset(gse49454$pheno, subjectID == mySubject)
    samplesToUse = c(samplesToUse, rownames(miniPheno)[miniPheno$visitID == min(miniPheno$visitID)])
  }
}
length(samplesToUse) == length(unique(gse49454$pheno$subjectID)) # true, 82 samples

# Get down to the first visit for each sample
pheno1visit = gse49454$pheno[samplesToUse,]
table(pheno1visit$group)# 20 healthy, 62 SLE
sleDatasets_males$gse49454 = subsetGEMFromPheno(gse49454, pheno1visit)


### remove males
quickViolin(gse49454, 'RPS4Y1', 'sex') # clean separation
phenoClean = subset(pheno1visit, sex == 'female')
sleDatasets$gse49454 = subsetGEMFromPheno(gse49454, phenoClean)
table(sleDatasets$gse49454$class) # 17 healthy, 53 SLE

#..................#
##### gse88884  - WB Clinical Trial #####
#..................#
# Background:
#   - Clinical trial of some kind of Ab


sleDatasets$gse88884 = testingObjectPub$originalData$GSE88884

### Clean pheno
# Make the column names match what cleanUpPheno expects
sleDatasets$gse88884$rawPheno = sleDatasets$gse88884$pheno
myColnames= colnames(sleDatasets$gse88884$rawPheno)
myColnames[grepl(pattern = "Sample_", x = myColnames)] = strsplitVector(myColnames[grepl(pattern = 'Sample_', x=myColnames)], 'Sample_', 2)
colnames(sleDatasets$gse88884$rawPheno) = myColnames

sleDatasets$gse88884$pheno = cleanUpPheno(sleDatasets$gse88884$rawPheno, T)
table(sleDatasets$gse88884$pheno$time) # all at baseline
table(sleDatasets$gse88884$pheno$group) # 60 HC, 1760 SLE

### Expr
#   - No NAs
#   - No negatives
#   - No obvious batch effect
range(sleDatasets$gse88884$expr) # 0.8-15
any(is.na(sleDatasets$gse88884$expr)) # false

# No obvious batch effect
png('gse88884_expr.png', width = 1500)
boxplot(sleDatasets$gse88884$expr, main = 'gse88884')
dev.off()

### Give better expr rownames
probeIDs = paste('probe', names(sleDatasets$gse88884$keys), sep='')
names(sleDatasets$gse88884$keys) = probeIDs
rownames(sleDatasets$gse88884$expr) = probeIDs
checkDataObject(sleDatasets$gse88884, 'Dataset') # True

### Sex labels
# Some samples with unknown sex labeled as '---'
# change '--' to 'unknown'
# I thought this woudl fix quickViolin, but it didn't 
table(sleDatasets$gse88884$pheno$sex, sleDatasets$gse88884$pheno$group)
sleDatasets$gse88884$pheno$sex = as.character(sleDatasets$gse88884$pheno$sex)
sleDatasets$gse88884$pheno$sex[sleDatasets$gse88884$pheno$sex == '--'] = 'unknown'

'XIST' %in% sleDatasets$gse88884$keys #T
'RPS4Y1' %in% sleDatasets$gse88884$keys #T
'KDM5D' %in% sleDatasets$gse88884$keys #T

# quickViolin isn't working for 'sex' label
#quickViolin(myDataset = sleDatasets$gse88884,myGene =  'RPS4Y1',columnName =  'sex')
quickViolin(myDataset = sleDatasets$gse88884,myGene =  'RPS4Y1',columnName =  'group') # clearly males in both healthy and SLE
quickViolin(myDataset = sleDatasets$gse88884,myGene =  'IFI44',columnName =  'group') # IFN gene up as expected

### Manual sex plot
sexGenes = getSampleLevelGeneData(sleDatasets$gse88884, c('XIST', 'RPS4Y1', 'KDM5D'))
dim(sexGenes)
rownames(sexGenes)

myData = sleDatasets$gse88884$pheno
myData$XIST = unlist(sexGenes['XIST',])
myData$RPS4Y1 = unlist(sexGenes['RPS4Y1',])
myData$KDM5D = unlist(sexGenes['KDM5D',])

# about 5 samples with a given sex label are scrambled 
ggplot(myData, aes(x=sex, y=RPS4Y1)) + geom_jitter(width=0.1)


### impute sex
# 1) replace unknown with imputed sex
# 2) remove samples whose imputed sex does not match labeled sex

# Impute sex
imputedSex = imputeSex(sleDatasets$gse88884)
sex = sleDatasets$gse88884$pheno$sex
table(imputedSex, sex) # matches mostly, 7 mismatches

# Replace unknowns with imputed sex
sex[sex == 'unknown'] = imputedSex[sex=='unknown']
sleDatasets$gse88884$pheno$sex = sex

# identify and remove mislabeled samples
badSamples = rownames(sleDatasets$gse88884$pheno)[sex !=imputedSex]
for(mySample in badSamples){
  sleDatasets$gse88884 = removeOneSample(sleDatasets$gse88884, mySample)
}

### Standardize group
table(sleDatasets$gse88884$pheno$group)
group = as.character(sleDatasets$gse88884$pheno$group)
group[group == "Normal"] = 'healthy'
sleDatasets$gse88884$pheno$group = group
table(sleDatasets$gse88884$pheno$group, sleDatasets$gse88884$pheno$sex)


### Save males
sleDatasets_males$gse88884 = sleDatasets$gse88884

### Remove males
phenoClean = subset(sleDatasets$gse88884$pheno, sex == 'female')
sleDatasets$gse88884 = subsetGEMFromPheno(sleDatasets$gse88884, phenoClean)

#..................#
##### GSE50772 - PBMC #####
#..................#
# Purpose: 
#   - Need clean PBMC dataset to see of CD99 is up in SLE PBMCs
# Background:
#   - Dataset seeing if IFN signature tracked with disease severity
#   - Looks to be adults 
#   - Half of controls are males
#
# Results:
#   - CD99 significantly up in SLE PBMCs
gse50772 = testingObjectPub$originalData$GSE50772

gse50772$rawPheno = gse50772$pheno
gse50772$pheno = cleanUpPheno(gse50772$rawPheno, T)
table(gse50772$pheno$group) # 11 female controls, 59 female SLE


### Check expr
#   - No NAs
#   - some negative values
#   - No obvious batch effect
range(gse50772$expr) # -3 to 16
boxplot(gse50772$expr, main = 'gse50772')
gse50772$expr = gse50772$expr + abs(min(gse50772$expr)) + 1.01

### Check sex labels
'XIST' %in% gse50772$keys # T
'RPS4Y1' %in% gse50772$keys # T
'KDM5D' %in% gse50772$keys # T

# Sex labels are correct
quickViolin(gse50772, 'RPS4Y1', 'group')
quickViolin(gse50772, 'XIST', 'group')

# imputeSex 
imputedSex = imputeSex(gse50772)
table(imputedSex, gse50772$pheno$group)
gse50772$pheno$sex = imputedSex

### CD99?
quickViolin(gse50772, 'CD99', 'group')

CD99 = unlist(getSampleLevelGeneData(gse50772, 'CD99'))

# with raising expr p = 0.042
# without raising expr, p = 0.042
t.test(CD99[gse50772$pheno$group == 'female Control'],
       CD99[gse50772$pheno$group == 'female SLE'])


phenoFem = subset(gse50772$pheno, sex == 'female')
gse50772_fem = subsetGEMFromPheno(gse50772, phenoFem)

sleDatasets$gse50772 = gse50772_fem

#..................#
##### GSE61635 - WB #####
#..................# 
# Purpose: 
#   - Is CD99 consistently up in SLE PBMCs??
#   - This cohort is actually WB, I emailed the authors

gse61635 = testingObjectPub$originalData$GSE61635

# Multiple visits
gse61635$rawPheno = gse61635$pheno
gse61635$pheno = cleanUpPheno(gse61635$rawPheno, T)

# Get subjectID
#subjectID = strsplitVector(gse61635$pheno$, 'V', 1)


### Impute sex
gse61635$pheno$sex = imputeSex(gse61635)
gse61635$pheno$group = paste(gse61635$pheno$sex, gse61635$pheno$diagnosis)


### Quick check for genes
quickViolin(gse61635, 'RPS4Y1', 'diagnosis')
quickViolin(gse61635, 'CD99', 'group') # CD99 is reduced, so looks like

gse61635_cells = MetaIntegrator::immunoStatesDecov(list(originalData=list(gse61635=gse61635)))
hist(gse61635_cells$immunoStates$gse61635$Correlation) # pretty high correlations
hist(gse61635_cells$immunoStates$gse61635$neutrophil) # 20-50% neutrophil

gse61635_cells$originalData$gse61635$pheno$granulocyte = gse61635_cells$immunoStates$gse61635$granulocyte

ggplot(gse61635_cells$originalData$gse61635$pheno, aes(x=group, y=granulocyte)) + geom_violin(fill='grey')+
  geom_jitter(size=2, width=0.1)


hist(gse24706_cells$immunoStates$gse24706_fem_clean$neutrophil, main='neut') # < 15% neutrophil


### Check with sexMetaObj
sexMetaObj_cells = MetaIntegrator::immunoStatesDecov(sexMetaObj)

granulocytes = c()
cohort = c()

for(i in 1:length(sexMetaObj$originalData)){
  myTitle = sexMetaObj$originalData[[i]]$formattedName
  #hist(sexMetaObj_cells$immunoStates[[i]]$granulocyte, main = myTitle)
  granulocytes = c(granulocytes, sexMetaObj_cells$immunoStates[[i]]$granulocyte)
  cohort = c(cohort, rep(myTitle, length(sexMetaObj$originalData[[i]]$class)))
  
}

test = as.data.frame(cbind(granulocytes, cohort))
test$granulocytes = as.numeric(as.character(test$granulocytes))

ggplot(test, aes(x=cohort, y=granulocytes)) + geom_violin(fill='grey') + 
  geom_jitter(width=0.1)


# PBMCs
#   - gse60491 --> neut --> Most < 5%, with long tail to 25%
#   - gse47353 --> neut --> most < 6%, looks like nomral dist around 3%
#   - gse21862 --> neut --> Most < 0.5%, with long tail to 3.5%
# WB
#   - gse53195 --> neut --> normalish distr around 10%, range 0-20%
#   - gse19151 --> neut --> Most between 5% and 25%, roughly normal, spread well
#   - gse17065 --> neut --> normalish around 10%, most below 15%


#..................#
##### GSE81622 - PBMC #####
#..................#
# Purpose: 
#   - SEe if CD99 is up in PBMCs in SLE
#   - Weak directionality, but not convincing

gse81622 = getGEOData('GSE81622')
gse81622 = gse81622$originalData$GSE81622

### Clean up pheno
gse81622$rawPheno = gse81622$pheno
gse81622$pheno =cleanUpPheno(gse81622$rawPheno, T)

# Group labels
group = strsplitVector(gse81622$pheno$title, '-', 1)
table(group)
gse81622$pheno$group = group
table(group, gse81622$pheno$sex)


# create group2
gse81622$pheno$group2 = paste(gse81622$pheno$sex, gse81622$pheno$group)
table(gse81622$pheno$group2)

# Create class vector
myClass = ifelse(gse81622$pheno$group == 'normal control', yes = 0, no = 1)
table(myClass, gse81622$pheno$group)
names(myClass) = rownames(gse81622$pheno)
gse81622$class = myClass

# Check expr
#   - No negatives
#   - no NAs
#   - No obvious batch effect
range(gse81622$expr) # 6-14
boxplot(gse81622$expr, main = 'gse81622')

### Sex labels
'XIST' %in% gse81622$keys # t
'RPS4Y1' %in% gse81622$keys # t
'KDM5D' %in% gse81622$keys # f

# ~ 4 mislabeled samples
quickViolin(gse81622, 'XIST', 'sex')
quickViolin(gse81622, 'RPS4Y1', 'sex')
quickScatter(gse81622, 'RPS4Y1', 'XIST', 'sex')

# Remove 5 bad samples
imputedSex = imputeSex(gse81622)
badSamples = rownames(gse81622$pheno)[gse81622$pheno$sex != imputedSex]
badSamples # 5 bad samples
for(mySample in badSamples){
  gse81622 = removeOneSample(gse81622, mySample)
}

### Sanity check genes
quickViolin(gse81622, 'MX1', 'group2') # up in SLE
quickViolin(gse81622, 'IFI44', 'group2') # up in SLE


### CD99
quickViolin(gse81622, 'CD99', 'group2')
CD99 = unlist(getSampleLevelGeneData(gse81622, "CD99"))
t.test(CD99[gse81622$pheno$group2 == 'female normal control'], 
       CD99[gse81622$pheno$group2 == "female SLE patient"])


violinPlot(filter3gene, gse81622, 'group2')

### subset to only female
phenoFem = subset(gse81622$pheno, sex == 'female')
gse81622_fem = subsetGEMFromPheno(gse81622, phenoFem)
sleDatasets$gse81622 = gse81622_fem

#...........................................#
##### Look at 3 genes in all female SLE #####
#...........................................#
# Questions:
#   - What do ZBED1, CYBB, and CD99 look like in all female SLE?
#        = CD99: unexpectedly down in SLE
#        = ZBED1: It's a mess
#        = CYBB: Up as expected
#   - How many datasets have a sig score difference?
#        = 3/5 have a sig score difference
#   - What does this mean?
# 
# Conclusions:
#   - I'd hoped that CD99 would be up on SLE blood
#   - But, it is messy and unclear
#   - Even after looking at an additional CD99 dataset, 
#     there isn't clear incrase in CD99

### Meta-analyze female SLE

# Prep and run meta-analysis
sleFem= list(originalData = sleDatasets)
checkDataObject(sleFem, 'Meta', "Pre-Analysis") # True
sleFem = runMetaAnalysis(sleFem, runLeaveOneOutAnalysis = F)


# Positive control genes
forestPlot(sleFem, "MX1") # up in all datasets
forestPlot(sleFem, "IFI44") # Up in all datasets


# Check out the three genese
forestPlot(sleFem, "CD99") # up in gse49454, down in 3 others
forestPlot(sleFem, 'ZBED1') # up in 2, mess in 3
forestPlot(sleFem, 'CYBB') # Fairly cnonsistently up in SLE


# Make my filter
filter3gene= sexMetaObj$filterResults$xy
filter3gene$posGeneNames = c('POOPSNOTAGENE')
filter3gene$negGeneNames = c('CD99', 'ZBED1', 'CYBB')

# What about score? 
violinPlot(filter3gene, sleDatasets$gse50635, 'group') # p = 0.24
violinPlot(filter3gene, sleDatasets$gse39088, 'group') # p = 0.11 (weak up)
violinPlot(filter3gene, sleDatasets$gse22098_adult, 'group') # p = 0.073 (down, right direction)
violinPlot(filter3gene, sleDatasets$gse49454, 'group') # p < 5E-5 (down, right direction)
violinPlot(filter3gene, sleDatasets$gse88884, 'group') # p < 1E-6, (down, right direction)


# What about age?
age = as.character(sleDatasets_males$gse88884$pheno$age_at_baseline)
age = as.numeric(age)
hist(age)
sleDatasets_males$gse88884$pheno$age = age

# Females lower than males until ~40 yo
quickScatterPheno(sleDatasets_males$gse88884, 'CD99', 'age', 'sex') + 
  geom_smooth() +scale_color_manual(values = c(femColor, malColor)) +
  xlim(c(18, 70))

# Difference seems to increase with age
quickScatterPheno(sleDatasets_males$gse88884, 'CYBB', 'age', 'sex') + 
  geom_smooth() +scale_color_manual(values = c(femColor, malColor)) +
  xlim(c(18, 70))

# larger difference, seems to increase with age
quickScatterPheno(sleDatasets_males$gse88884, 'ZBED1', 'age', 'sex') + 
  geom_smooth() +scale_color_manual(values = c(femColor, malColor)) +
  xlim(c(18, 70))

# Too few controls to properly compare HC vs SLE with age
quickScatterPheno(sleDatasets$gse49454, 'CD99', 'age', 'group') +
  geom_smooth()

quickScatterPheno(sleDatasets_males$gse49454, 'CD99', 'age', 'sex') +
  geom_smooth()


myData = sleDatasets_males$gse88884$pheno
myData$CD99 = unlist(getSampleLevelGeneData(sleDatasets_males$gse88884, 'CD99'))
myData$CYBB = unlist(getSampleLevelGeneData(sleDatasets_males$gse88884, 'CYBB'))


#.................................#
##### SLE risk and cell types #####
#.................................#
# Purpose: 
#    - Last ditch effort to connect sex diff to SLE risk
#    - Is there a difference in cell type proportions? 

# Load dataset
load("blood/gse24706/gse24706.RData")

quickViolin(gse24706_fem_clean, 'CD99', 'group') # messy up in SLE
quickViolin(gse24706_fem_clean, 'CYBB', 'group') # more clearly up in SLE

gse24706_cells = MetaIntegrator::immunoStatesDecov(list(originalData = list(gse24706_fem_clean = gse24706_fem_clean)))

# Check quality
hist(gse24706_cells$immunoStates$gse24706_fem_clean$Correlation) # some values below 0.7

myData = gse24706_cells$originalData$gse24706_fem_clean$pheno
myData$correlation = gse24706_cells$immunoStates$gse24706_fem_clean$Correlation
myData$monocytes = gse24706_cells$immunoStates$gse24706_fem_clean$monocyte
myData$CD4 = gse24706_cells$immunoStates$gse24706_fem_clean$CD4_positive_alpha_beta_T_cell
myData$granulocyte = gse24706_cells$immunoStates$gse24706_fem_clean$granulocyte

# reorder group 
x =  factor(myData$group, levels = unique(myData$group)[c(3,1,4,2,5)])
myData$group = x
all(as.character(x) == as.character(gse24706_fem_clean$pheno$group))
gse24706_fem_clean$pheno$group = x


# Decide what threshold to cut off correlation at
hist(myData$correlation)


# SLE has a large # of samples below 0.7 correlation
ggplot(myData, aes(x=group, y= correlation)) + geom_violin(trim=F, fill='grey')+
  geom_jitter(width=0.1) + ggtitle("gse24706 - SLE and Healthy Risk Groups") + xlab('')
  geom_abline(slope=0, intercept = 0.4)

#### Monocytes
# Weakly higher in higher risk people
ggplot(myData, aes(x=group, y= monocytes)) + geom_violin(trim=F, fill='grey')+
  geom_jitter(width=0.1, aes(col=correlation < 0.5), size = 3) +
  scale_color_manual(values = c('dodgerblue', 'red'))

### CD4+ T cells
# First degree relatives seem to have reduced CD4 levels
ggplot(myData, aes(x=group, y= CD4)) + geom_violin(trim=F, fill='grey')+
  geom_jitter(width=0.1, aes(col=correlation < 0.5), size = 3) +
  scale_color_manual(values = c('dodgerblue', 'red'))

### Checking granulocytes, in case they help explain low correlation SLE people
ggplot(myData, aes(x=group, y= granulocyte)) + geom_violin(trim=F, fill='grey')+
  geom_jitter(width=0.1, aes(col=correlation < 0.5), size = 3) +
  scale_color_manual(values = c('dodgerblue', 'red'))

# CD4/Monocyte ratio is pretty gosh darn interesting, seems lower in 
# all SLE risk groups
ggplot(myData, aes(x=group, y= CD4/monocytes)) + geom_violin(trim=F, fill='grey')+
  geom_jitter(width=0.1, aes(col=correlation < 0.5), size = 3) +
  scale_color_manual(values = c('dodgerblue', 'red')) + geom_abline(slope = 0, intercept = 1)


###
myData$ratio = myData$CD4/myData$monocytes
t.test(myData$ratio[myData$group == 'HC ANA-'], myData$ratio[myData$group == 'FDR ANA+']) # 0.01
t.test(myData$ratio[myData$group == 'HC ANA-'], myData$ratio[myData$group == 'FDR ANA-']) # 0.077
t.test(myData$ratio[myData$group == 'HC ANA-'], myData$ratio[myData$group == 'HC ANA+']) # 0.08
t.test(myData$ratio[myData$group == 'HC ANA-'], myData$ratio[myData$group == 'SLE ANA+']) # 0.6
t.test(myData$ratio[myData$group == 'HC ANA-'], myData$ratio[myData$group == 'SLE ANA+' & myData$correlation > 0.5]) # 0.015



### gse49454 has directly measured stuff!
# I can't find monocytes, so it's stupid!
myData4 = sleDatasets$gse49454$pheno

ggplot(myData4, aes(x= skin, y=`cd4 tcells`/leukocytes)) + geom_violin()

for(myDataset in sleDatasets){
  print(myDataset$formattedName)
  print(colnames(myDataset$pheno))
}



#...............................#
##### Plots for Lab Meeting #####
#...............................#



### Winn's score as proof of concept it's worth looking at these gorups
sleFilter$posGeneNames
pdf('plots/gse24706_winnSig.pdf', width=8, height = 5)
violinPlot(sleFilter, gse24706_fem_clean, 'group') + xlab('') +
  ylab('SLE Score (93 genes)') + theme_cowplot(font_size = 18)
dev.off()

sum(length(sleFilter$posGeneNames) + length(sleFilter$negGeneNames)) # 93 genes

### Autosomal-iSEXS doesn't work
pdf('plots/gse24706_autoiSEXS.pdf', width=8, height=5)
violinPlot(sexMetaObj$filterResults$autosomeOnly, gse24706_fem_clean, 'group') + xlab('') +
  ylab('Autosomal-iSEXS') + theme_cowplot(font_size = 18)
dev.off()

### XY-iSEXS works!
pdf('plots/gse24706_xyiSEXS.pdf', width=8, height=5)
violinPlot(sexMetaObj$filterResults$xy, gse24706_fem_clean, 'group') + xlab('') +
  ylab('XY-iSEXS') + theme_cowplot(font_size = 18)
dev.off()

### Dissect XY-iSEXS
filterXup = sexMetaObj$filterResults$xy
filterXup$negGeneNames = c("POOPSNOTAGENE")

filterY = sexMetaObj$filterResults$xy
filterY$posGeneNames = c('POOPSNOTAGENE')
filterY$negGeneNames = filterY$negGeneNames[!filterY$negGeneNames %in% filter3gene$negGeneNames]
filterY$negGeneNames

pdf('plots/gse24706_Xfem.pdf', width=8, height=5)
violinPlot(filterXup, gse24706_fem_clean, 'group') + xlab('') +
  ylab('X genes higher in females') + theme_cowplot(font_size = 18)
dev.off()

pdf('plots/gse24706_Y.pdf', width=8, height=5)
violinPlot(filterY, gse24706_fem_clean, 'group') + xlab('') +
  ylab('Y Chromosome Genes') + theme_cowplot(font_size = 18)
dev.off()

pdf('plots/gse24706_Xmal.pdf', width=8, height=5)
violinPlot(filter3gene, gse24706_fem_clean, 'group') + xlab('') +
  ylab('X-Chr genes higher in males') + theme_cowplot(font_size = 18)
dev.off()
filter3gene


filterX = filterXup
filterX$negGeneNames = filter3gene$negGeneNames
violinPlot(filterX, gse24706_fem_clean, 'group') + xlab('') +
  ylab('X Chromosome genes') + theme_cowplot(font_size = 18)

#..................................#
##### FilterX in PBMC cohorts: #####
#..................................#
# PBMC cohorts: 
#   - gse81622
#   - gse50772

### filterX works
violinPlot(filterX, gse81622_fem, 'group') # p = 0.064 reduced in SLE patient
violinPlot(filterX, gse50772_fem, 'group')# sig reduced

### 3 genes are garbage
violinPlot(filter3gene, gse81622_fem, 'group') +theme_cowplot(font_size = 18) +
  xlab('') + ylab('X-chr genes up in males')# garbage
violinPlot(filter3gene, gse50772_fem, 'group') +theme_cowplot(font_size = 18)+ # garbage
  xlab('') + ylab('X-chr genes up in males')


### CD99 - Ugg, I don't know
quickViolin(gse81622_fem, 'CD99', 'group') # messy
quickViolin(gse50772_fem, 'CD99', 'group') # sig up


### CYBB
quickViolin(gse81622_fem, 'CYBB', 'group') # Up in SLE
quickViolin(gse50772_fem, 'CYBB', 'group') # weak up in SLE

### filterXup works well
violinPlot(filterXup, gse81622_fem, 'group') # Down in both SLE and lupus nephritis
violinPlot(filterXup, gse50772_fem, 'group')# sig reduced


#.......................#
##### SLE and iSEXS #####
#.......................#
# Purpose: 
#   - How does XY-iSEXS and Autosomal-iSEXS change with SLE proper?
#
# Results: 
#   - XY-iSEXS
#       - Consistently lower in SLE
#       - Most of the data is WB
#   - Autosomal-iSEXS
#       - Consistently lower in SLE
#       - even more consistent than XY-iSEXS

names(sleFem$originalData)

### XY-iSEXS
violinPlot(sexMetaObj$filterResults$xy, sleFem$originalData$gse88884, 'group') # WB - sig lower, but, messy
violinPlot(sexMetaObj$filterResults$xy, sleFem$originalData$gse49454, 'group') # WB - sig lower
violinPlot(sexMetaObj$filterResults$xy, sleFem$originalData$gse81622, 'group') # PBMC - sig lower in SLE, but not LN
violinPlot(sexMetaObj$filterResults$xy, sleFem$originalData$gse50772, 'group') # PBMC - No difference
violinPlot(sexMetaObj$filterResults$xy, sleFem$originalData$gse39088, 'group') # WB - sig lower
violinPlot(sexMetaObj$filterResults$xy, sleFem$originalData$gse50635, 'group') # WB - p = 0.075, lower
violinPlot(sexMetaObj$filterResults$xy, sleFem$originalData$gse22098_adult, 'group') # WB - sig lower

### Autosomal-iSEXS
violinPlot(sexMetaObj$filterResults$autosomeOnly, sleFem$originalData$gse88884, 'group') # WB -sig lower
violinPlot(sexMetaObj$filterResults$autosomeOnly, sleFem$originalData$gse49454, 'group') # WB - sig lower
violinPlot(sexMetaObj$filterResults$autosomeOnly, sleFem$originalData$gse81622, 'group') # PBMC - sig lower in SLE and LN
violinPlot(sexMetaObj$filterResults$autosomeOnly, sleFem$originalData$gse50772, 'group') # PBMC - sig lower
violinPlot(sexMetaObj$filterResults$autosomeOnly, sleFem$originalData$gse39088, 'group') # WB - sig lower
violinPlot(sexMetaObj$filterResults$autosomeOnly, sleFem$originalData$gse50635, 'group') # WB - sig lower
violinPlot(sexMetaObj$filterResults$autosomeOnly, sleFem$originalData$gse22098_adult, 'group') # WB - sig lower


#...............................#
##### SLE and Deconvolution #####
#...............................#
# Purpose: 
#   - What do CD4 and monocyte levels look like in SLE? 
#
# Results: 
#   - gse81622 - PBMC
#        - Monocytes: Significantly More monocytes in SLE
#        - CD4: significantly reduced in SLE
#        - CD4/monocytes: significantly reduced
#   - gse49454
#        - Monocytes: Not quite sig, but higher in SLE
#        - CD4 significantly reduced
#        - CD4/monocyte: significantly reduced
#                        skin looks ineresting, but otherwise no clear pattern with SLEDAI or anything

#sleFem_cells = MetaIntegrator::immunoStatesDecov(metaObject = sleFem)

### gse49454
myData = sleFem$originalData$gse49454$pheno

ggplot(myData, aes(x=sledai, y=`cd4 tcells`)) + geom_point(size = 2, alpha=0.5)

# correlation between CD4 levels and SLEDAI
ggplot(myData, aes(x=sledai, y=`cd4 tcells`/leukocytes)) + geom_point(size = 2, alpha=0.5)
cor.test(myData$sledai, myData$`cd4 tcells`/myData$leukocytes)

gse49454_cells = MetaIntegrator::immunoStatesDecov(list(originalData = list(gse49454 = sleFem$originalData$gse49454)))
myData = gse49454_cells$originalData$gse49454$pheno
myData$correlation = gse49454_cells$immunoStates$gse49454$Correlation
myData$cd4 = gse49454_cells$immunoStates$gse49454$CD4_positive_alpha_beta_T_cell
myData$monocyte = gse49454_cells$immunoStates$gse49454$monocyte

hist(myData$correlation) # most samples above 0.7

# Monocytes - non quite sig, but higher in SLE
ggplot(myData, aes(x=group, y=monocyte)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(width=0.1, size=2, aes(col=correlation < 0.6))
t.test(myData$monocyte[myData$group == 'healthy'],
       myData$monocyte[myData$group == 'SLE']) # p = 0.098

# CD4 - Significantly reduced
ggplot(myData, aes(x=group, y=cd4)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(width=0.1, size=2, aes(col=correlation < 0.6))
t.test(myData$cd4[myData$group == 'healthy'],
       myData$cd4[myData$group == 'SLE']) # p = 0.00039

# CD4/monocyte ratio - significantly reduced
ggplot(myData, aes(x=group, y=cd4/monocyte)) + geom_violin(fill='grey', trim=F) +
  geom_jitter(width=0.1, size=2, aes(col=correlation < 0.6))
myData$ratio = myData$cd4/myData$monocyte
t.test(myData$ratio[myData$group == 'healthy'],
       myData$ratio[myData$group == 'SLE']) 

### Ratio correlations

# SLEDAI - messy and hard to interpret
ggplot(myData, aes(x=sledai, y=cd4/monocyte)) + geom_point(size=2, aes(col = correlation< 0.6)) +
  geom_abline(slope = 0, intercept = 1)

# Years since diagnosis - messy and hard to interpret
ggplot(myData, aes(x=`years since diagnosis of sle`, y=cd4/monocyte)) + geom_point(size=2, aes(col = correlation< 0.6)) +
  geom_abline(slope = 0, intercept = 1)

# A-DSD - messy and hard to interpret
ggplot(myData, aes(x=`a-dsd`, y=cd4/monocyte)) + geom_point(size=2, aes(col = correlation< 0.6)) +
  geom_abline(slope = 0, intercept = 1)

### Ratio in different subgroups
# Skin - Doesn't look all that differet
ggplot(myData, aes(x=skin, y=cd4/monocyte)) + geom_violin(fill='grey', trim=F) +
  geom_jitter(width=0.1, size=2, aes(col=correlation < 0.6))

ggplot(myData, aes(x=joints, y=cd4/monocyte)) + geom_violin(fill='grey', trim=F) +
  geom_jitter(width=0.1, size=2, aes(col=correlation < 0.6))

ggplot(myData, aes(x=renal, y=cd4/monocyte)) + geom_violin(fill='grey', trim=F) +
  geom_jitter(width=0.1, size=2, aes(col=correlation < 0.6))


### gse81622
gse81622_cells = MetaIntegrator::immunoStatesDecov(list(originalData = list(gse81622 = sleFem$originalData$gse81622)))
myData = gse81622_cells$originalData$gse81622$pheno
myData$monocyte = gse81622_cells$immunoStates$gse81622$monocyte
myData$correlation = gse81622_cells$immunoStates$gse81622$Correlation
myData$cd4 = gse81622_cells$immunoStates$gse81622$CD4_positive_alpha_beta_T_cell

# Check out correlation
hist(myData$correlation) # only one sample below 0.7, so it's all good!

# Significantly more monocytes
ggplot(myData, aes(x=group, y= monocyte)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(size=2, width=0.1, aes(col=age > 40)) + ggtitle('gse81622 - PBMC')
t.test(myData$monocyte[myData$group == 'normal control'], 
       myData$monocyte[myData$group == 'SLE patient']) # p < 0.001
t.test(myData$monocyte[myData$group == 'normal control'], 
       myData$monocyte[myData$group == 'SLE patient with lupus nephritis']) # p= 0.021

# significantly fewer CD4 T cells
ggplot(myData, aes(x=group, y= cd4)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(size=2, width=0.1) + ggtitle('gse81622 - PBMC')
t.test(myData$cd4[myData$group == 'normal control'], 
       myData$cd4[myData$group == 'SLE patient']) # p < 0.018
t.test(myData$cd4[myData$group == 'normal control'], 
       myData$cd4[myData$group == 'SLE patient with lupus nephritis']) # p= 0.048

# ratio - significantly monocyte skewed
ggplot(myData, aes(x=group, y= cd4/monocyte)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(size=2, width=0.1) + ggtitle('gse81622 - PBMC') + geom_abline(slope = 0, intercept = 1)
myData$ratio = myData$cd4/myData$monocyte
t.test(myData$ratio[myData$group == 'normal control'], 
       myData$ratio[myData$group == 'SLE patient']) # p  = 0.00045
t.test(myData$ratio[myData$group == 'normal control'], 
       myData$ratio[myData$group == 'SLE patient with lupus nephritis']) # p= 0.0017

