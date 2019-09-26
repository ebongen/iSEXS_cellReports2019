# 5/10/2018 Preprocess gse49454

# Purpose:
#   - gse49454 has cell proportions in SLE patients
#   - I can see if XY score changes in SLE patients with CD4+T cell levels

#....................#
##### Load Stuff #####
#....................#

# Load libraries
library(MetaIntegrator)

# Load version of cleanUpPheno that can unscramble
source("/labs/khatrilab/ebongen/sexDifferences2.0/isexs/general_GEMfunx.R")
source("~/Tools/metaintegrator_public/R/cleanUpPheno.R")

# Load graphing functions
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/5_autoImmunity/SLE/blood/gse49454/")

#................................#
##### Download and fix pheno #####
#................................#
# Purpose:

# Download
gse49454 = getGEOData("GSE49454")
gse49454 = gse49454$originalData$GSE49454 

# Clean up pheno
gse49454 = cleanUpPheno(gse49454)

# Doublecheck pheno
View(gse49454$pheno[1:50,])

# Check age column - It's good!
hist(gse49454$pheno$age)

# Create sex
sex = ifelse(gse49454$pheno$gender == 'F', yes = 'female', no = 'male')
gse49454$pheno$sex = sex
gse49454$pheno$gender = NULL

#....................#
##### Check expr #####
#....................#
# Purpose: 
#   - Any NAs? Nope!
#   - Log transformed? Yup!
#   - No NAs
#   - No obvious batch effect
#        - Not able to make successful boxplot
#        - But, median, 1st Q, and 3rd Q are nearly all the same
sum(is.na(gse49454$expr))

# Range is 0-14
range(gse49454$expr)

# pdf isnt' really working
pdf(file= 'expr.pdf', width= 14, height=7)
boxplot(gse49454$expr, main = "gse49454")
dev.off()

### Medians are all about the same
allMedians = apply(gse49454$expr, 2, median)
hist(allMedians)
length(unique(allMedians))

### Only one value for 1st quartile
all1st = apply(gse49454$expr, 2, function(x) summary(x)[2])
hist(all1st)
length(unique(all1st))

### Only one value for 3rd quartile
all3rd = apply(gse49454$expr, 2, function(x) summary(x)[5])
length(unique(all3rd))

#..........................#
##### Check sex labels #####
#..........................#
# Purpose: 
#   - Remove samples where sex is mislabeled
#
# Reulst: 
#   - No samples were mislabeled

'XIST' %in% gse49454$keys # T
'RPS4Y1' %in% gse49454$keys # T
'KDM5D' %in% gse49454$keys # F

# No mislabeled samples
quickViolin(gse49454, "XIST", "sex") # a little overlap
quickViolin(gse49454, "RPS4Y1", 'sex') # clean separation

#............................#
##### Subset to Baseline #####
#............................#
# Purpose:
#   = Chaussabel has crappy sample labeling
#   = I need to only take one visit per patient

### Create a group column that takes into account both sex and disease
group2 = paste(gse49454$pheno$sex, gse49454$pheno$group)
group2 = gsub(pattern = "Healthy control of SLE", replacement = "Healthy", x = group2)
group2 = sapply(group2, simpleCap)
gse49454$pheno$group2 = group2

### Create a visit ID
visitID = strsplitVector(gse49454$pheno$title, "V", 2)
visitID[gse49454$pheno$group != "SLE"] = "1"
gse49454$pheno$visitNumber = as.numeric(as.character(visitID))
visitID = paste("Visit", visitID)
gse49454$visitID = factor(visitID, levels = sort(unique(visitID)))

### Create a sample ID
subjID = strsplitVector(gse49454$pheno$title, "V", 1)
subjID[subjID == 'C3'] = as.character(gse49454$pheno$title)[subjID=="C3"]
gse49454$pheno$subjectID = subjID


### Get the earliest sample per subject
earliestSamples = c()

for(mySubject in unique(gse49454$pheno$subjectID)){
  myPheno = subset(gse49454$pheno, subjectID == mySubject)
  earliestSamples[mySubject] = rownames(myPheno)[myPheno$visitNumber == min(myPheno$visitNumber)]
  
}

### Create baseline gse49454
phenoBaseline = gse49454$pheno[earliestSamples,]
gse49454_baseline = subsetGEMFromPheno(gse49454, phenoBaseline)

#..............#
##### Save #####
#..............#

# Contains males, females, HC, and SLE
# I'll subset how I want during analysis
table(gse49454$pheno$group, gse49454$pheno$sex)

table(gse49454_baseline$pheno$group, gse49454_baseline$pheno$sex)

save(gse49454, gse49454_baseline, file = "gse49454.RData")
