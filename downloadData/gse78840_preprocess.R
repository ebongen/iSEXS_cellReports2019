# Preprocess gse78840
# December 2017

# Background
# eQTL study in CD4s and CD8s in Estonians
# Healthy donors from an Estonian university 

## Future directions:
# Download from raw to get XY gene expression
# Divvy by age

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

# Load packages
library(MetaIntegrator)
library(data.table)
library(GEOquery)

# Load functions
source("00_tools/general_GEMfunx.R") # code for preprocessing datasets
source("00_tools/plot_funx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

# Load sexMetaObj
load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/4_sortedCells/cd4_Tcells/gse78840/")

#..................#
##### Download #####
#..................#
# Purpose: 
# Download cohort from GEO
# Clean up pheno


# Download
gse78840 = getGEOData("GSE78840")
gse78840 = gse78840$originalData$GSE78840

gse78840$rawPheno = gse78840$pheno
gse78840$pheno = cleanUpPheno(gse78840$rawPheno, T)

# Remove Duplicate columns
gse78840$pheno$`age:ch1` = NULL
gse78840$pheno$`gender:ch1` = NULL
gse78840$pheno$`sample.section:ch1` = NULL
gse78840$pheno$`sentrix.barcode:ch1` = NULL
gse78840$pheno$`cell type:ch1` = NULL
gse78840$pheno$`tissue:ch1` = NULL

# Extract subject IDs
subjectIDs = as.character(gse78840$pheno$title)
subjectIDs = strsplitVector(subjectIDs, "SAMPLE_", 2)
subjectIDs = strsplitVector(subjectIDs, "_", 1)
subjectIDs = strsplitVector(subjectIDs, " ", 1)
subjectIDs = paste("SUBJECT_", subjectIDs, sep = "")
gse78840$pheno$subjectIDs = subjectIDs

# Extract sample IDs
sampleIDs = as.character(gse78840$pheno$title)
sampleIDs = strsplitVector(sampleIDs, " RNA", 1)
gse78840$pheno$sampleIDs = sampleIDs
#............................#
##### Non-normalized expr #####
#............................#
# Background: 
#  The automatically downloaded data only contains non-sex chromosome genes
#  On GEO the thing that's listed as "GSE78840_non_normalized.txt" is actually
#  the complete expression data

# Purpose: 
#  Create expr and keys for full platform

# Grab non normalized data
nonNorm = fread("GSE78840_non_normalized.txt")
dim(nonNorm) # 47k vs 1307

# Create new expr using sample IDs
all(sampleIDs %in% colnames(nonNorm)) # True
expr = nonNorm[,..sampleIDs]
expr = as.matrix(expr)
rownames(expr) = nonNorm$ID_REF


# Change the column names to GSMs
gsmKey = rownames(gse78840$pheno)
names(gsmKey) = as.character(gse78840$pheno$sampleIDs)
colnames(expr) = gsmKey[colnames(expr)]

# Is it in log2? Nope!
summary(expr[,1:5]) # Not log2
min(expr) # 59
expr = log2(expr)

# Are there NAs?
# Nope!
sum(is.na(expr))

### Filter based on probe p-values
### I decided not to do this becuase this removed all of the highest expressed genes
### If I knew more about pre-processing, then maybe I could figure out how to do it right
pValCol = colnames(nonNorm)[grepl(pattern = "Pval", x = colnames(nonNorm))]
expr_pval = nonNorm[,..pValCol]
expr_pval = as.matrix(expr_pval)
rownames(expr_pval) = nonNorm$ID_REF
range(expr_pval) # between 0 and 1
dim(expr_pval) # 47k genes
dim(expr_pval[rownames(gse78840$originalExpr),]) # 38k
range(expr_pval[rownames(gse78840$originalExpr),]) # 38k, also ranges between 0 and 1
# allNonsig = apply(expr_pval, 1, function(x) all(x > 0.05))
# expr = expr[allNonsig,]

# No obvious batch effect
# Severely needs some normalization
png("boxplot_allSamp.png", width = 1200)
boxplot(expr, main = "Non-normalized expr - gse78840")
dev.off()

### Create keys
gpl10558 = getGEO(as.character(gse78840$pheno$platform_id[1]))
keyDict = gpl10558@dataTable@table$Symbol
names(keyDict) = gpl10558@dataTable@table$ID

newKeys = keyDict[rownames(expr)]

### Update the Dataset object
gse78840$originalExpr = gse78840$expr
gse78840$expr = expr
gse78840$exp_comment = "Obtained from GSE789840_non_normalized.txt"
gse78840$originalKeys = gse78840$keys
gse78840$keys = newKeys
gse78840$key_comment = "Manually created to match manual expr"

checkDataObject(gse78840, "Dataset")


gse78840 = GEM_normalizer(gse78840)

png("boxplot_allSamp_normalized.png", width = 1200)
boxplot(gse78840$expr, main = "Normalized Non-normalized expr - gse78840")
dev.off()
#.................................#
##### Sex Annotation and Keys #####
#.................................#
# Purpose: 
# Make sure keys look good
# Check sex annotation

# Results
# While 4 samples were flagged by imputeSex, it was unclear whether
# they were mislabeled because expression of XIST and RPS4Y1 was
# misleading
# So, I looked at the full XY ISEXS score, and those samples
# look firmly to belong to their originally annotated sex

# Keys look good
gse78840$keys[1:50]
length(gse78840$keys) # 47k
sum(is.na(gse78840$keys)) # 0 NAs

# Check for genes
"XIST" %in% gse78840$keys # True
"RPS4Y1" %in% gse78840$keys # True
"KDM5D" %in% gse78840$keys # False

# Mostly clear separation, some overlapping samples
quickViolin(gse78840, "XIST", "sex")
quickViolin(gse78840, "RPS4Y1", "sex")
quickScatter(gse78840, "RPS4Y1", "XIST", "sex")

# Impute sex
imputedSex = imputeSex(gse78840)
gse78840$pheno$imputedSex = imputedSex
gse78840$pheno$sexLabelsAgree = gse78840$pheno$sex == imputedSex
table(imputedSex, gse78840$pheno$sex)


quickScatter(gse78840, "RPS4Y1", "XIST", "sex")
quickScatter(gse78840, "RPS4Y1", "XIST", "cell_type") # CD4's and CD8's
quickScatter(gse78840, "RPS4Y1", "XIST", "imputedSex")

# The one male adn the 3 females that imputed Sex disagrees with sex labels on are hard to call
# 1 "Male" labeled "Female-expression" = low XIST and low RPS4Y1
# 3 "Female-labeled" with "male-expression" = female-like XIST and lower end of male-normal RPS4Y1 
quickScatter(gse78840, "RPS4Y1", "XIST", "sexLabelsAgree")

## Look at entier XY score, if XIST and RPS4Y1 alone are messy
xyScore = calculateScore(sexMetaObj$filterResults$xy, gse78840)
gse78840$pheno$xyScore = xyScore
ggplot(data=gse78840$pheno, aes(x=sex, y=xyScore, col= sexLabelsAgree)) + geom_violin() + geom_jitter(width = 0.1)
ggplot() + geom_jitter(data = subset(gse78840$pheno, sexLabelsAgree== TRUE), aes(x=sex, y=xyScore)) + 
  geom_jitter(data = subset(gse78840$pheno, sexLabelsAgree!= TRUE), aes(x=sex, y=xyScore), colour = "red", size = 5) +
  ggtitle("gse78840 - XY ISEXS Score - Red Dots were flagged by Imputation of sex")


# Create class vector 
gse78840$class = createClassVector(gse78840$pheno$sex, "female", gse78840$pheno)

#............................#
##### Divvy by cell type #####
#............................#
# Results: 
#   In the version that's posted to GEO, the XY chromosome genes are missing
#   and the autosomal sig works well in CD4's (AUC = 0.8), but fails entirely
#   in CD8's and WB (AUC ~0.5)
#
#   In the version that I have here, where I grab all the probes
#   It works comparably in CD4's and CD8's (AUC 0.6) and quite well
#   in WB (0.77)

# Using the original expr, AUC = 0.8, autosomal only
# But, using normalized "non-normalized" data, I get AUC 0.6
phenoCD4 = subset(gse78840$pheno, cell_type == "CD4+ T cells")
range(phenoCD4$age) # 22-84
t.test(phenoCD4$age[phenoCD4$sex == "female"], phenoCD4$age[phenoCD4$sex == "male"]) # p=0.5
gse78840_cd4 = subsetGEMFromPheno(gse78840, phenoCD4)
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse78840_cd4) # 
hist(gse78840_cd4$pheno$age) # even coverage from 20s-85


# Autosomal sig gets an AUC 0.6 in CD8+ T cells
phenoCD8 = subset(gse78840$pheno, cell_type == "CD8+ T cells")
t.test(phenoCD8$age[phenoCD8$sex == "female"], phenoCD8$age[phenoCD8$sex == "male"]) # p=0.75
gse78840_cd8 = subsetGEMFromPheno(gse78840, phenoCD8)
gse78840_cd8$formattedName = "GSE78840 CD8+ T cell"

hist(phenoCD8$age) # even coverage from 20s to 80s
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse78840_cd8) # AUC = 0.6, autosomal only


# Autosomal sig works in WB, AUC 0.77
phenoWB = subset(gse78840$pheno, cell_type == "Peripheral blood")
gse78840_wb = subsetGEMFromPheno(gse78840, phenoWB)
gse78840_wb$formattedName = "GSE78840 WB"
hist(phenoWB$age)# ~30 young people, and ~40 people over 70
table(phenoWB$age <=40, phenoWB$sex) # 19 young females, 16 young males
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse78840_wb) # AUC = 0.77



#..............#
##### Save #####
#..............#

save(gse78840_cd4, gse78840_cd8, gse78840_wb, file = "gse78840_noXY.RData")
