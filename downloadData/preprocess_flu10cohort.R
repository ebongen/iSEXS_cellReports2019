# Preprocess Sandra's Cytof flu vaccination study
# December 2017

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")



source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/4_sortedCells/cellProp/sandra/")

#...............#
##### Pheno #####
#...............#
# Read in phenotype information
subjData = read.csv("SubjectData.csv", skip = 2)

# Create twin pair ID
twinPairID = paste(subjData$Group, subjData$no, sep = "")
subjData$twinPairID = twinPairID

# Create sampleIDs
sampleIDs = as.character(subjData$SLVP.ID)
sampleIDs = str_replace_all(sampleIDs, "-", "_")
rownames(subjData)= sampleIDs

# sex column 
unique(subjData$Gender)
sex = as.character(subjData$Gender)
sex = ifelse(sex == "Female", yes = "female", no = "male")
subjData$sex = sex
table(subjData$Gender, subjData$sex)
subjData$Gender = NULL

# rename age column
subjData$age = subjData$Age
subjData$Age = NULL

#..............#
##### Expr #####
#..............#
# Read in cellprop information
cellProp = read.csv("Twin_CyTOF_Complete.csv")
rownames(cellProp) = str_replace_all(cellProp$SLVP.ID, "-", "_")
cellProp$SRI.ID =NULL
cellProp$SLVP.ID = NULL

expr= t(as.matrix(cellProp))

# rownames and column names agree
all(colnames(expr) == rownames(subjData))

#...............................................#
##### Create extra parts for Dataset object #####
#...............................................#
# Create class vector
myClass = createClassVector(subjData$sex, "female", subjData)

# create keys
myKeys = rownames(expr)
names(myKeys) = myKeys
#....................................#
##### Construct a Dataset Object #####
#....................................#
# Create Dataset Object
flu10Cohort = list(pheno = subjData,
                   expr = expr, 
                   class = myClass,
                   keys = myKeys,
                   formattedName = "flu10Cohort_HIMC")

# It's looking good!
checkDataObject(flu10Cohort, "Dataset") # True

#.......................#
##### Subset by age #####
#.......................#
range(flu10Cohort$pheno$age) # 11.5 - 59
hist(flu10Cohort$pheno$age)

# 8 females, 6 males
table(flu10Cohort$pheno$age < 18, flu10Cohort$pheno$sex)

# 17 females vs 9 males [18,40]
table(flu10Cohort$pheno$age <=40 & flu10Cohort$pheno$age >18, flu10Cohort$pheno$sex)

# 26 females, 15 males
table(flu10Cohort$pheno$age > 40, flu10Cohort$pheno$sex)


phenoYoung = subset(flu10Cohort$pheno, age <=40 & age >=18)
range(phenoYoung$age) # 19-30
t.test(phenoYoung$age[phenoYoung$sex == "female"], phenoYoung$age[phenoYoung$sex == "male"]) # p=0.2
flu10Cohort_young = subsetGEMFromPheno(flu10Cohort, phenoYoung)

#..............#
##### Save #####
#..............#
save(flu10Cohort, flu10Cohort_young, file = "flu10Cohort.RData")
