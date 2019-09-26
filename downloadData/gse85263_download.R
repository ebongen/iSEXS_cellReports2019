# 9/28/2018  - Preprocess gse85263

#....................#
##### Load Stuff #####
#....................#
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Load libraries
library(MetaIntegrator)
library(cowplot)

# Load functions
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

# Load iSEXS
load('isexs/sexMetaObj.RData')

# Load Aditya's datasets
SouthAmericaDatasets = readRDS("/labs/khatrilab/adityamr/random/SouthAmericaDatasets.rds")

setwd('0_datasets/1_validation/gse85263/')
#...............#
##### Fix up #####
#...............#

# Fix the phenos
for(i in 1:length(SouthAmericaDatasets)){
  colnames(SouthAmericaDatasets[[i]]$pheno) = gsub(pattern = 'Sample_', replacement = '', x = colnames(SouthAmericaDatasets[[i]]$pheno))
  SouthAmericaDatasets[[i]]$rawPheno = SouthAmericaDatasets[[i]]$pheno
  SouthAmericaDatasets[[i]]$pheno = cleanUpPheno(SouthAmericaDatasets[[i]]$rawPheno, T)
}

gse85263 = SouthAmericaDatasets$gse85263

### Subset to baseline
phenoBase = subset(gse85263$pheno, time == 'Baseline')
dim(phenoBase) # 19 people
range(phenoBase$age) # 19-41
sum(phenoBase$age > 40) # 1 person
phenoBase = subset(phenoBase, age <= 40)
table(phenoBase$sex) # 11 female, 7 male

### Does not pass muster = sample names in expr and pheno don't match
checkDataObject(gse85263, 'Dataset')
sum(rownames(gse85263$pheno) %in% colnames(gse85263$expr)) # not all here??

sampleIDs = strsplitVector(gse85263$pheno$title, 'sample_', 2)
all(sampleIDs == colnames(gse85263$expr)) # true!
colnames(gse85263$expr) = rownames(gse85263$pheno) # We can give expr teh GSM colnames now
names(gse85263$class) = rownames(gse85263$pheno)

# Yay!
gse85263_base = subsetGEMFromPheno(gse85263, phenoBase)


# Check for batch effect
# No obvious batch effect
boxplot(gse85263_base$expr, main = 'gse85263')

# Check sex labels 
'RPS4Y1' %in% gse85263_base$keys # True
quickViolin(gse85263_base, 'RPS4Y1', 'sex') # perfect separation
quickViolin(gse85263_base, 'XIST', 'sex') # perfect separation

# Significantly higher in females!
violinPlot(autoSig, gse85263_base, 'sex')


### Add class vector
gse85263_base$class = createClassVector(gse85263_base$pheno$sex, 'female', gse85263_base$pheno)

save(gse85263_base, file = 'gse85263_colombia.RData')
