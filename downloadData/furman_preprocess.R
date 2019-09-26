# Download and preprocess Furman's Testosterone data

# Purpose: 
#   - Does iSEXS correlate with testosterone levels?
#   - I'm pretty sure it doesn't, but I can't find that 
#     original analysis
#   - Download and incorporate testosterone data

#....................#
##### Load Stuff #####
#....................#
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Load packages
library(MetaIntegrator)
library(ggplot2)
library(cowplot)

# Source functions
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

load('isexs/sexMetaObj.RData')

setwd('0_datasets/7_hormone/gse41080_furman/')


testo = read.csv('furman_testo_data.csv')
#............................#
##### Download and Pheno #####
#............................#

# Ugg!!! There is no gene expression data in the series
gse41080 = getGEOData('GSE41080')
gse41080 = gse41080$originalData$GSE41080
  
stevenData = readRDS(file = '/labs/khatrilab/ebongen/friends/steven/fluCombo.rdat')

'GSE41080' %in% names(stevenData)

#...............................#
##### Do I have it already? #####
#...............................#

# Prep expr
expr = read.delim('/labs/khatrilab/ebongen/GSE_Datasets/GSE41080/GSE41080_expr_normalized.txt')
dim(expr) 
length(expr$ProbeID) == length(unique(expr$ProbeID)) # all are unique
rownames(expr) = as.character(expr$ProbeID)
expr$ProbeID = NULL
expr = as.matrix(expr)
range(expr) # 6-16, no NAs

# NO obvvious batch effect
boxplot(expr, main = 'gse41080')

#...............#
##### Pheno #####
#...............#


pheno = read.delim('/labs/khatrilab/ebongen/GSE_Datasets/GSE41080/GSE41080_pheno.txt')
rownames(pheno) = pheno$GSM
phenoClean = cleanUpPheno(pheno, T)
phenoClean$sex = tolower(phenoClean$sex)
#View(phenoClean)

# They overlap!
sum(as.character(testo$xID) %in% as.character(phenoClean$id)) # 78 overlap
testoDict = testo$T_nmol_l
names(testoDict) = as.character(testo$xID)

# Add testosterone values to phenoClean
phenoClean$T_nmol_l = testoDict[as.character(phenoClean$id)]

# Testosterone measured for both males and females
table(is.na(phenoClean$T_nmol_l), phenoClean$sex)

#..............#
##### Keys #####
#..............#
# GEt the probe IDs
gpl6747 = read.csv('/labs/khatrilab/ebongen/GSE_Datasets/GSE41080/GPL6947.csv')

# To get up-to-date probes, download a dataset that uses the same gpl
gse15530 = getGEOData('GSE15530')
myKeys = gse15530$originalData$GSE15530$keys[rownames(expr)]
length(myKeys) # 48k
sum(is.na(myKeys)) # 24k

#.................................#
##### Assemble Dataset Object #####
#.................................#


# Create class
myClass = createClassVector(phenoClean$sex, 'female', phenoClean)


# Assemble into Dataset object!
gse41080 = list(pheno = phenoClean, 
                expr = expr,
                class = myClass,
                keys = myKeys,
                formattedName = 'GSE41080', 
                rawPheno = pheno,
                exp_comment = 'Series matrix file saved from 2014?',
                key_comment = 'Taken from GSE15530'
                )

# It passes!
checkDataObject(gse41080, 'Dataset')

#..........................#
##### Check Sex Labels #####
#..........................#
# Purpose: 
#   - Make sure given sex labels match transcription
# 
# Results:
#  - one male and one female swapped
#
# Conclusions:
#  - Pheno annotation is trustworthy

'XIST' %in% myKeys   # T
'RPS4Y1' %in% myKeys # T
'KDM5D' %in% myKeys # F

# One male and one female might be swapped
quickViolin(gse41080, 'XIST', 'sex')
quickViolin(gse41080, 'RPS4Y1', 'sex')

# Impute sex
imputedSex =imputeSex(gse41080)
badSamples = rownames(gse41080$pheno)[imputedSex != gse41080$pheno$sex]
length(badSamples) # 2 bad samples

# Remove bad samples
for(mySample in badSamples){
  gse41080 = removeOneSample(gse41080, mySample)
}

#..............................................#
###### Check Testosterone Responsive Genes #####
#..............................................#
# Purpose: 
#   - Try to use published genes that are responsive to testosterone
#     to sanity check the testosoterone values Furman gave me
#   - Used a publication that highlighted genes that were up or down regulated
#     by testosterone based on a mouse model and based on a mouse cell line system
#
# Results: 
#   - AKAP1
#        - From 2012 publication PMID: 22002438 
#        - Human Protein Atlas
#            - Protein level: Maxed out in nearly all tissues, including females
#            - mRNA level: Maxed out in testes, midling in others, low in female tissues
#        - Not much of a sex difference
#        - Messy positive correlation in males with T levels (p=0.06)
#
#   - PHKG2
#        - From 2012 publication PMID: 22002438 
#        - Human Protein Atlas
#            - Protein level: Not measured
#            - mRNA level: Maxed out in testes, tiny in all others, very low in female tissues 
#        - Not much of a sex difference in expression
#        - Way non-significant trend where having higher T --> higher PHKG2 (maybe cuz not expressed strongly in blood)
#   - MRPS6
#        - From 2012 publication PMID: 22002438 
#        - Publications says it has a negative relationship with T
#        - Human Protein Atlas shows that it's a ribosomal subunit that's expreseed everywhere
#        - Way non-significant trend for people with highest T having lower MRPS6
#    - COX6B1
#        - From 2012 publication PMID: 22002438 
#        - Publications says it has a negative relationship with T  
#        - Human Protein Atlas: expressed reasonably in blood and female tissues, strong expression in male tissues
#        - Negative correlation between T and COX6B1: p=0.068
#
# Conclusions:
#    - The gene list I'm working with is from mouse and non-blood tissues
#    - So, the fact that the trends tend to be in the right direction is heartening!



# It the gene present?
'COX6B1' %in% myKeys # present! 

# Not that different between males and females, maybe higher in males
quickViolin(gse41080, 'PHKG2', 'sex')

quickScatterPheno(myDataset = gse41080, myGene = 'AKAP1', 'T_nmol_l', 'sex') # Assoc. with higher expressin with higher T


#### Test genes!
myGene = 'COX6B1'
myData = gse41080$pheno
myData$myGene = unlist(getSampleLevelGeneData(gse41080, myGene))

ggplot(myData, aes(x=T_nmol_l, y=myGene, col = sex)) + 
  geom_point(size = 2) + theme_cowplot() + ylab(myGene)

myData_male = subset(myData, sex == 'male')
cor.test(myData_male$T_nmol_l, myData_male$myGene) # r = 0.34, p = 0.0963

#..............#
##### Save #####
#..............#

save(gse41080, file = 'gse41080_testosterone.RData')
