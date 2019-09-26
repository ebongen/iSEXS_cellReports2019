# August 2017 Preprocess GSE13485

# Background:
# gse13485 is Bali Palendran's Yellow Fever vaccine timecourse study
# It only has gene mappings given in UniGene, so I had to convert to Gene Symbol
# PBMC in people aged 18-45
# I need to ask for ages

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
library(GEOquery)
library(RMySQL)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/1_validation/gse13485_PBMC/")

ageLabels = read.csv("yellowFeverVacc_age_sex.csv")

#..................#
##### Download #####
#..................#
# Get Dataset object
gse13485 = getGEOData("GSE13485")
gse13485 = gse13485$originalData$GSE13485

# Get keys
gse13485_getGEO = getGEO("GSE13485")
keys = gse13485_getGEO$GSE13485_series_matrix.txt.gz@featureData@data

#...............#
##### Pheno #####
#...............#
# Age and sex not supplied in GEO
gse13485$rawPheno = gse13485$pheno
gse13485$pheno = cleanUpPheno(gse13485$rawPheno, T)
View(gse13485$pheno)

# Subject ID
subjID = as.character(gse13485$pheno$title)
subjID = strsplitVector(subjID, " Day", 1)
subjID = strsplitVector(subjID, "ID ", 2)
gse13485$pheno$subjectID = subjID

# Add age
ageDict = ageLabels$Age
names(ageDict) = ageLabels$Donor
ageDict = na.omit(ageDict)
gse13485$pheno$age = ageDict[gse13485$pheno$subjectID]

# Add sex
sexDict = ageLabels$Gender
sexDict = ifelse(sexDict == "F", yes = "female", no = "male")
names(sexDict) = ageLabels$Donor
sexDict = na.omit(sexDict)

gse13485$pheno$sex = sexDict[gse13485$pheno$subjectID]

#..............#
##### Expr #####
#..............#
# Expr is in good shape!

# Good range
# No obvious batch effect
boxplot(gse13485$expr, main = "gse13485")

# Positive value, good!
min(gse13485$expr) # 2.78

#......................................#
##### Semi-Manually Deal with Keys #####
#......................................#
library(org.Hs.eg.db)


#org.Hs.egSYMBOL is an R object that provides mappings between entrez gene identifiers and gene abbreviations
#         entrez gene ID --> gene abbrev
abbrevAnnot <- toTable(org.Hs.egSYMBOL)
geneSymbolDict = as.character(abbrevAnnot$symbol)
names(geneSymbolDict) = as.character(abbrevAnnot$gene_id)

# org.Hs.egUNIGENE is an R object that provides mappings between entrez gene identifiers and #UniGene identifiers.
#         entrez gene ID --> unigene
unigeneAnnot <- toTable(org.Hs.egUNIGENE)
unigeneDict<- as.character(unigeneAnnot$gene_id)
names(unigeneDict) <- as.character(unigeneAnnot$unigene_id)

# Convert from Unigene, to Entrez ID, to Gene Symbol
keysUnigene = as.character(keys$UniGene)
keysEntrezID = unigeneDict[keysUnigene]
keysSymobl = geneSymbolDict[keysEntrezID]

# Add them to the keys matrix
keys$EntrezID = keysEntrezID
keys$GeneSymbol = keysSymobl


# Create keys vector
newKeys = keys$GeneSymbol
names(newKeys) = rownames(keys)

all(rownames(gse13485$expr) == names(newKeys)) # True

gse13485$keys = newKeys[rownames(gse13485$expr)]

checkDataObject(gse13485, "Dataset") # True

#...................................#
##### Create Validation Dataset #####
#...................................#

# Take only baseline samples from people <40 years old
phenoVal = gse13485$pheno
phenoVal = subset(phenoVal, Time == 0) # Only baseline samples
phenoVal = subset(phenoVal, age <=40)

# Create Validation dataset
gse13485_validation = subsetGEMFromPheno(gse13485, phenoVal)
gse13485_validation$formattedName = "GSE13485_validation"
checkDataObject(gse13485_validation, "Dataset") # True


# Create Class vector
gse13485_validation$class = createClassVector(gse13485_validation$pheno$sex, "female", gse13485_validation$pheno)
table(gse13485_validation$class)

#....................................#
##### Make sure sex labels match #####
#....................................#

"XIST" %in% gse13485$keys # False
"RPS4Y1" %in% gse13485$keys # False
"KDM5D" %in% gse13485$keys # Triue

# All but one of the samples have consistently High or Low 
# expression of KDM5D across samples
quickViolin(gse13485_validation, "KDM5D", columnName = "sex")


# 1 male and 1 female are probably swapped
imputedSex = imputeSex(gse13485_validation)
all(imputedSex == gse13485_validation$pheno$sex) # False

# Remove samples where imputed sample does not match labeled sex
mislabeledSamples = rownames(gse13485_validation$pheno)[which(gse13485_validation$pheno$sex != imputedSex)]
for(mySample in mislabeledSamples){
  gse13485_validation = removeOneSample(gse13485_validation, mySample)
}

# Looks good!
quickViolin(gse13485_validation, "KDM5D", "sex")


# How many swapped samples are there total?
imputedSex = imputeSex(gse13485)
table(imputedSex, gse13485$pheno$sex) # 15 mislabeled samples total

#..............#
##### Save #####
#..............#

save(gse13485, file = "gse13485_full.RData")
save(gse13485_validation, file = "gse13485_validation.RData")