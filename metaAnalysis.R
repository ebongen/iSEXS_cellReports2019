# Run meta-analysis in Discovery and Validation

#....................#
##### Load Stuff #####
#....................#
# Set working directory
myDir = "/labs/khatrilab/ebongen/sexDifferences2.0/isexs/"
setwd(myDir)

# Load libraries
library(MetaIntegrator)
library(ggplot2)
library(cowplot)

# Source code
source("general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

### Initialize constants
# The index number for a figure, so that they'll be listed in /plots in the order I make them
# Allows me to insert new figures without having to change plot names manually
# Does not correspond directly to "Figure 1" of the actual paper
figIndex = 0 
nReps = 750
#..............................#
##### 1 Meta-Analysis #####
#..............................#
# Purpose: 
#   - Meta-analyze Discovery cohorts
#   - Choose thresholds for isexs
#
# In the future:
#   - Adjust filter to avoid garbage genes
#         - smaller ES filter, LOO=T

# Load Discovery datasets
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/0_discovery/")
load("gse17065/gse17156_morocco.RData")
load("gse19151/gse19151_thromboEmbolism.RData")
load("gse21862/gse21862_china.RData")
load("gse47353/gse47353.RData")
load("gse53195/gse53195.RData")
load("gse60491/gse60491_UCLA.RData")

setwd(myDir)

# Create list of datasets for Old ladies vs young ladies
agingData =list(gse17065 = gse17065,
                gse53195 = gse53195)

# Remove extra datasets
rm(gse17065, gse19151, gse21862, gse47353, gse53195, gse60491)

discDatasets = list(gse17065 = gse17065_discovery, 
                    gse19151 = gse19151_discovery,
                    gse21862 = gse21862_discovery, 
                    gse47353 = gse47353_discovery, 
                    gse53195 = gse53195_discovery, 
                    gse60491 = gse60491_discovery)

# Fix formatted names
discDatasets$gse17065$formattedName = "GSE17065"
discDatasets$gse19151$formattedName = "GSE19151"
discDatasets$gse21862$formattedName = "GSE21862"
discDatasets$gse47353$formattedName = "GSE47353"
discDatasets$gse53195$formattedName = "GSE53195"
discDatasets$gse60491$formattedName = "GSE60491"

# Run meta-analysis with LOO
sexMetaObj = MetaIntegrator::runMetaAnalysis(list(originalData = discDatasets), runLeaveOneOutAnalysis = T)

# Filter genes 
test = filterGenes(sexMetaObj, isLeaveOneOut = TRUE, effectSizeThresh = 0.4, FDRThresh = 0.05, numberStudiesThresh = 1)
summarizeFilterResults(test, getMostRecentFilter(test))
length(test$filterResults[[1]]$posGeneNames)
length(test$filterResults[[1]]$negGeneNames)

# Save the distribution plot for some reason
distPlot = recordPlot()
pdf('plots/discovery_ESdistribution.pdf', width=7, height = 5)
distPlot
dev.off()

### LOO = T thresholds
# ES 0.4, FDR 5% --> 21 up, 15 down
# ES 0.4, FDR 10% --> 27 up, 20 down
# ES 0.4, FDR 15% --> 33 UP, 24 down
# ES 0.3, FDR 5% --> 31 up, 18 down
# ES 0.3, FDR 10% --> 48 up, 26 down

### LOO = F thresholds, nStudies = 2
# ES 0.4, FDR 5% --> 94 up, 50 down
# ES 0.5, FDR 5% --> 35 up, 28 down
# ES 0.5, FDR 10% --> 41 up, 35 down


# Create isexs
test = filterGenes(sexMetaObj, isLeaveOneOut = FALSE, effectSizeThresh = 0.4, FDRThresh = 0.05, numberStudiesThresh = 2)
sexMetaObj$filterResults$isexs = test$filterResults[[1]]
checkDataObject(sexMetaObj, "Meta", "Post-Filter") # True
isexs = c(sexMetaObj$filterResults$isexs$posGeneNames, sexMetaObj$filterResults$isexs$negGeneNames)
posGenes = sexMetaObj$filterResults$isexs$posGeneNames
negGenes = sexMetaObj$filterResults$isexs$negGeneNames
length(isexs)

#..............................#
##### 2 Validation Cohorts #####
#..............................#
# Purpose:
#   - Load Validation dataset objects
#   - Meta-analyize them for validation of isexs

# Load datasets
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/1_validation/")
load("gse13485_PBMC/gse13485_validation.RData")
load("gse18323_PBMC/gse18323_validation.RData")
load("gse19442/gse19442_latentTB.RData")
load("gse21311/gse21311_email.RData")
load("gse30453_PBMC/gse30453_validation.RData")
load("gse37069/gse37069.RData")
load("gse38484/gse38484_validation.RData")
load("gse61821/gse61821.RData")
load("gse65219_PBMC/gse65219_validation.RData")
load("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/2_age/postMenopause/gse58137/gse58137_young.RData")
load('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/2_age/postMenopause/gse58137/gse58137_timecourse.RData')
load('gse85263/gse85263_colombia.RData')
load('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/2_age/postMenopause/sdy212/sdy212.RData')
# go back to isexs directory
setwd(myDir)

# Save datasets good for aging
agingData$gse21311 = gse21311_email
agingData$gse38484 = gse38484_healthy
agingData$gse58137 = gse58137gpl6
agingData$gse65219 = gse65219
agingData$sdy212 = sdy212

# Remove extra datasets
rm(gse19442) # remove extra one
rm(gse21311_email, gse21311_email_over50)
rm(gse37069)
rm(gse38484, gse38484_healthy)
rm(gse58137gpl6)
rm(gse61821_eden, gse61821_SEAICRN, gse61821_SEAICRN_val)
rm(gse65219, gse65219_90s)



# Make list of validation obj
valDatasets = list(gse13485 = gse13485_validation, 
                   gse18323 = gse18323_validation, 
                   gse19442 = gse19442_val, 
                   gse21311 = gse21311_email_validation, 
                   gse30453 = gse30453_validation, 
                   gse37069= gse37069_validation, 
                   gse38484 = gse38484_validation, 
                   gse58137 = gse58137gpl6_young,
                   gse61821 = gse61821_val,
                   gse65219 = gse65219_validation, 
                   gse85263 = gse85263_base)

# Fix formatted Names
valDatasets$gse13485$formattedName = "GSE13485"
valDatasets$gse18323$formattedName = "GSE18323"
valDatasets$gse19442$formattedName = "GSE19442"
valDatasets$gse21311$formattedName = "GSE21311"
valDatasets$gse30453$formattedName = "GSE30453"
valDatasets$gse37069$formattedName = "GSE37069"
valDatasets$gse38484$formattedName = "GSE38484"
valDatasets$gse58137$formattedName = "GSE58137"
valDatasets$gse61821$formattedName = 'GSE61821'
valDatasets$gse65219$formattedName = "GSE65219"
valDatasets$gse85263$formattedName = "GSE85263"

valMetaObj = list(originalData = valDatasets)

# Run meta-analysis of validation cohorts
# Don't do LOO becuase we won't filter genes
valMetaObj = runMetaAnalysis(valMetaObj, runLeaveOneOutAnalysis = FALSE)

# Save distribution
distPlot = recordPlot()
pdf('plots/validation_ESdistribution.pdf', width=7, height = 5)
distPlot
dev.off()

#..................#
##### nSamples #####
#..................#

nValidation = 0
for(myDataset in valMetaObj$originalData){
  nValidation = nValidation + length(myDataset$class) 
}


#...............................#
##### 3 Chromosome Location #####
#...............................#
# Purpose: 
#   - Label each isexs gene by chromosome location
#        - Including PAR and X-escape
#   - Separate isexs into XY-sig and autosomal-sig

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/1_metaAnaly/chrLocation/")
# Load PAR genes
par1 = read.delim("PAR1_hugo.txt", header = T, sep = "\t")
par2 = read.delim("PAR2_hugo.txt", header = T, sep = "\t")

par1Genes = c(as.character(par1$Approved.Symbol))
par2Genes = c(as.character(par2$Approved.Symbol))

# Load X-escape genes
# From Tukiainen et al 2017, 
# Aviv Regev's group GTEX of X-inactivation across tissues
escapeGenes =read.csv("xist_escape_genes.csv")

# Create a dictionary of escape genes
escapeDict = as.character(escapeGenes$Combined.XCI.status)
escapeDict = paste("X_", escapeDict, sep = "")
names(escapeDict) = as.character(escapeGenes$Gene.name)
escapeDict = escapeDict[!escapeDict == "X_"]


setwd(myDir)

# Load dictionary of all gene locations
load('reference/allGeneLocations_dict.RData')

# Obtain ISEXS locations
ISEXSLoc = allGeneLocations_dict[isexs]
names(ISEXSLoc) = isexs

if(any(is.na(ISEXSLoc))){
  warning("isexs contains a gene of unknown chromosomal location")
}

# If BAGE5 is present, manually add it
if("BAGE5" %in% isexs){
  ISEXSLoc["BAGE5"] = "13" # Manually add in BAGE5's location
}


# Add 0's in front of single didget chromosomes (ensures correct order in histogram)
ISEXSLoc[which(ISEXSLoc %in% as.character(1:9))] = paste("0", ISEXSLoc[which(ISEXSLoc %in% as.character(1:9))], sep="")
ISEXSLoc

# Save vector of chromosome location
sexMetaObj$filterResults$isexs$chrLocation = ISEXSLoc

# Create detailed version that includes PAR and X-escape
ISEXSLoc_detailed = ISEXSLoc

# Add PAR1 annotation, but only to negGenes
if(any(negGenes %in% par1Genes)){
  # Find genes higher expressed in males that are on PAR1
  myPARgenes = negGenes[negGenes %in% par1Genes]
  
  # Give those genes the PAR1 label
  ISEXSLoc_detailed[myPARgenes] = "PAR1"
} 

# Add PAR2 annotation, but only to negGenes
if(any(negGenes %in% par2Genes)){
  # Identify negGenes that are on PAR2
  myPAR2genes = negGenes[negGenes %in% par2Genes]
  
  # Give those genes PAR2 annotation
  ISEXSLoc_detailed[myPAR2genes] = "PAR2"
}

# Add X escape annotation
if(any(posGenes %in% names(escapeDict))){
  myEscapeGenes = posGenes[posGenes %in% names(escapeDict)]
  ISEXSLoc_detailed[myEscapeGenes] = escapeDict[myEscapeGenes]
}


table(ISEXSLoc_detailed)
ISEXSLoc_detailed[ISEXSLoc_detailed == "X"] # CYBB, I can't find evidence pseudoautisomal region
escapeDict["CYBB"] # X inactive

# Create ISEXSLoc for graphing
# Convert all the "X_inactive" to just X
ISEXSLoc = ISEXSLoc_detailed
ISEXSLoc[ISEXSLoc == "X_inactive"] = "X" # Classify X_inactive as "X"
ISEXSLoc[ISEXSLoc == "X_variable"] = "X_escape"
sort(unique(ISEXSLoc))


# Save within sexMetaObj
sexMetaObj$filterResults$isexs$chrLocation = ISEXSLoc
sexMetaObj$filterResults$isexs$chrLocation_detailed = ISEXSLoc_detailed

# Create autosomal filter
autoSig = sexMetaObj$filterResults$isexs
autoSig$posGeneNames = autoSig$posGeneNames[!ISEXSLoc[autoSig$posGeneNames] %in% c('PAR1', 'X', 'X_escape', 'Y')]
autoSig$negGeneNames = autoSig$negGeneNames[!ISEXSLoc[autoSig$negGeneNames] %in% c('PAR1', 'X', 'X_escape', 'Y')]
autoSig$chrLocation = "Autosomes Only"
autoSig$chrLocation_detailed = NULL
sexMetaObj$filterResults$autoSig = autoSig

# Create XY filter
xySig = sexMetaObj$filterResults$isexs
xySig$posGeneNames = xySig$posGeneNames[ISEXSLoc[xySig$posGeneNames] %in% c('PAR1', 'X', 'X_escape', 'Y')]
xySig$negGeneNames = xySig$negGeneNames[ISEXSLoc[xySig$negGeneNames] %in% c('PAR1', 'X', 'X_escape', 'Y')]
xySig$chrLocation = "Sex Chromosomes Only"
xySig$chrLocation_detailed = NULL
sexMetaObj$filterResults$xySig = xySig

checkDataObject(sexMetaObj, "Meta", "Post-Filter") # True

#..............................................#
##### Discovery vs Validation Effect Sizes #####
#..............................................#
# Create a list of iSEXS genes
isexs = c(sexMetaObj$filterResults$isexs$posGeneNames, sexMetaObj$filterResults$isexs$negGeneNames)


# Try for LOO=T threshold
# ES = 0.3, FDR < 5% --> 49 genes
# ES = 0.3, FDR < 10% --> 74
# ES = 0.4, FDR < 10% --> 47 genes
# ES = 0.3, FDR < 100% --> 392 genes
# ES = 0.3, FDR < 15% --> 107 genes
test = filterGenes(sexMetaObj, isLeaveOneOut = T, effectSizeThresh = 0.3, FDRThresh = 0.10)
length(c(test$filterResults$FDR0.1_es0.3_nStudies1_looaTRUE_hetero0$posGeneNames, 
         test$filterResults$FDR0.1_es0.3_nStudies1_looaTRUE_hetero0$negGeneNames))

isexs2 = c(test$filterResults$FDR0.1_es0.3_nStudies1_looaTRUE_hetero0$posGeneNames,
           test$filterResults$FDR0.1_es0.3_nStudies1_looaTRUE_hetero0$negGeneNames)
sum(isexs2 %in% isexs) # 62
sum(!isexs2 %in% isexs) # 12 not in original isexs!


# Create data frame of iSEXS effect sizes
myData = cbind(sexMetaObj$metaAnalysis$pooledResults[isexs,'effectSize'], 
               valMetaObj$metaAnalysis$pooledResults[isexs,'effectSize'])
colnames(myData) = c('discovery', 'validation')
rownames(myData) = isexs
myData = as.data.frame(myData)

myLim = 0.8


a = ggplot(myData, aes(x=discovery, y=validation)) + 
  geom_point(size=2, alpha = 0.6)+
  theme_cowplot() + 
  xlab('Discovery - Summary Effect Sizes') +
  ylab('Validation - Summary Effect Sizes')

b = ggplot(myData, aes(x=discovery, y=validation)) + 
  geom_point(size=2, alpha = 0.6) +
  theme_cowplot() + xlim(c(-myLim,myLim)) + ylim(c(-myLim,myLim))+
  xlab('Discovery - Summary Effect Sizes') +
  ylab('Validation - Summary Effect Sizes') +
  geom_hline(yintercept = 0.4, col='black') + 
  geom_hline(yintercept = -0.4, col='black')+
  geom_hline(yintercept = 0.1, col='red') + 
  geom_hline(yintercept = -0.1, col='red')
 
plot_valVSdisc = plot_grid(a,b, nrow=1, ncol=2)  

save_plot('plots/disc_vs_val.pdf', plot_valVSdisc,ncol=2, nrow=1)


### Subset into positive and negative
myData_pos = subset(myData, discovery > 0)
myData_neg = subset(myData, discovery < 0)

badGenes = rownames(myData_pos)[myData_pos$validation < 0.1]
badGenes = c(badGenes, rownames(myData_neg)[myData_neg$validation > -0.1])
badGenes = na.omit(badGenes)
length(badGenes) #28
length(badGenes)/length(isexs2) # 19.4%

# 6 on X chromosome
# 22 on autosomes
table(ISEXSLoc[badGenes])

22/length(c(autoSig$posGeneNames, autoSig$negGeneNames)) # 20%
6/length(c(xySig$posGeneNames, xySig$negGeneNames)) # 17%

#...........................#
##### iSEXS Stats Table #####
#...........................#
# Purpose: 
#   - Make a table that gives: 
#       1) Gene Symbol
#       2) Chromosomal location
#       3) Discovery:
#            - Effect Size
#            - p-value
#            - FDR
#       4) Validation
#            - Effect Size
#            - p-value

# Extract Discovery info
discData = sexMetaObj$metaAnalysis$pooledResults[isexs,c('effectSize', 'effectSizePval', 'effectSizeFDR')]
colnames(discData) = paste("discovery_", colnames(discData), sep='')

# Extract Validation info
valData = valMetaObj$metaAnalysis$pooledResults[isexs,c('effectSize', 'effectSizePval')]
colnames(valData) = paste('validation_', colnames(valData), sep='')

# Combine them
sum(rownames(discData) == rownames(valData))
rownames(valData)[rownames(valData) != rownames(discData)] # Yup, there's an NA
rownames(valData) = isexs
chromosome = ISEXSLoc[isexs]
comboData = cbind(chromosome, discData, valData)

write.csv(comboData, 'reference/isexs_genes.csv', quote = F)


#...........................#
##### Tissue Annotation #####
#...........................#
# Purpose: 
#   - I need to separate PBMC and WB datasets in heatmap
#   - So, I need to standardize the tissue annotation


### Discovery
# Whole blood
sexMetaObj$originalData$gse17065$pheno$tissue = rep('Whole Blood', nrow(sexMetaObj$originalData$gse17065$pheno))
sexMetaObj$originalData$gse19151$pheno$tissue = rep('Whole Blood', nrow(sexMetaObj$originalData$gse19151$pheno))
sexMetaObj$originalData$gse53195$pheno$tissue = rep('Whole Blood', nrow(sexMetaObj$originalData$gse53195$pheno))

# PBMC
sexMetaObj$originalData$gse21862$pheno$tissue = rep('PBMC', nrow(sexMetaObj$originalData$gse21862$pheno))
sexMetaObj$originalData$gse47353$pheno$tissue = rep('PBMC', nrow(sexMetaObj$originalData$gse47353$pheno))
sexMetaObj$originalData$gse60491$pheno$tissue = rep('PBMC', nrow(sexMetaObj$originalData$gse60491$pheno))


### Validation
# Whole Blood
valMetaObj$originalData$gse19442$pheno$tissue = rep('Whole Blood', nrow(valMetaObj$originalData$gse19442$pheno))
valMetaObj$originalData$gse21311$pheno$tissue = rep('Whole Blood', nrow(valMetaObj$originalData$gse21311$pheno))
valMetaObj$originalData$gse37069$pheno$tissue = rep('Whole Blood', nrow(valMetaObj$originalData$gse37069$pheno))
valMetaObj$originalData$gse38484$pheno$tissue = rep('Whole Blood', nrow(valMetaObj$originalData$gse38484$pheno))
valMetaObj$originalData$gse58137$pheno$tissue = rep('Whole Blood', nrow(valMetaObj$originalData$gse58137$pheno))
valMetaObj$originalData$gse61821$pheno$tissue = rep('Whole Blood', nrow(valMetaObj$originalData$gse61821$pheno))

# PBMC
valMetaObj$originalData$gse13485$pheno$tissue = rep('PBMC', nrow(valMetaObj$originalData$gse13485$pheno))
valMetaObj$originalData$gse18323$pheno$tissue = rep('PBMC', nrow(valMetaObj$originalData$gse18323$pheno))
valMetaObj$originalData$gse30453$pheno$tissue = rep('PBMC', nrow(valMetaObj$originalData$gse30453$pheno))
valMetaObj$originalData$gse65219$pheno$tissue = rep('PBMC', nrow(valMetaObj$originalData$gse65219$pheno))
valMetaObj$originalData$gse85263$pheno$tissue = rep('PBMC', nrow(valMetaObj$originalData$gse85263$pheno))

#............................#
##### Country Annotation #####
#............................#
# Purpose: 
#   - Purvesh wants to add country annotation to heatmap
#   - So, I need to store that information in the Dataset object

### Discovery
sexMetaObj$originalData$gse17065$country = 'GSE17065 - Morocco'
sexMetaObj$originalData$gse19151$country = 'GSE19151 - USA'
sexMetaObj$originalData$gse21862$country = 'GSE21862 - China'
sexMetaObj$originalData$gse47353$country = 'GSE47353 - USA'
sexMetaObj$originalData$gse53195$country = 'GSE53195 - Australia'
sexMetaObj$originalData$gse60491$country = 'GSE60491 - UK'

### Validation
valMetaObj$originalData$gse13485$country = 'GSE13485 - USA'
valMetaObj$originalData$gse18323$country = 'GSE18323 - USA'
valMetaObj$originalData$gse19442$country = 'GSE19442 - South Africa'
valMetaObj$originalData$gse21311$country = 'GSE21311 - Australia'
valMetaObj$originalData$gse30453$country = 'GSE30453 - USA'
valMetaObj$originalData$gse37069$country = 'GSE37069 - USA'
valMetaObj$originalData$gse38484$country = 'GSE38484 - Netherlands/Denmark'
valMetaObj$originalData$gse58137$country = 'GSE58137 - USA'
valMetaObj$originalData$gse61821$country = 'GSE61821 - Singapore'
valMetaObj$originalData$gse65219$country = 'GSE65219 - Finland'
valMetaObj$originalData$gse85263$country = 'GSE85263 - Colombia'

#..........................................#
##### Younger Females vs Older Females #####
#..........................................#
# Purpose: 
#  - Run a meta-analysis of younger females vs older females
#  - Identify iSEXS genes that change with age in females

# Function that
#   - removes males
#   - Removes females 41-49 years old
#   - Remove females younger than 18
#   - Sets old females as controls
oldLadySetup = function(myDataset){
  # Subset pheno
  oldPheno = subset(myDataset$pheno, sex == 'female') # remove males
  oldPheno = subset(oldPheno, age <=40 | age >=50) # remove 41-49 yo
  oldPheno = subset(oldPheno, age >=18) # remove people younger than 18
  
  # Create new Dataset object for new pheno
  myDataset = subsetGEMFromPheno(myDataset, oldPheno)
  
  # Create new class vector
  myClass = ifelse(myDataset$pheno$age <50, yes = 1, no = 0)
  names(myClass) = rownames(myDataset$pheno)
  myDataset$class = myClass
  
  # Subset down to isexs genes
  myDataset$keys = myDataset$keys[myDataset$keys %in% isexs]
  myDataset$expr = myDataset$expr[names(myDataset$keys),]
  
  if(!checkDataObject(myDataset, 'Dataset')){
    warning('Subsetted Dataset object failed QC')
  }
  
  return(myDataset)
}

# Make formatted names better
agingData$gse17065$formattedName = 'GSE17065'
agingData$gse21311$formattedName = 'GSE21311'
agingData$gse58137$formattedName = 'GSE58137'

# Remove samples that won't be used in analysis
agingLadies = lapply(agingData, oldLadySetup)

# Test class vector and length of keys
for(myDataset in agingLadies){
  print(myDataset$formattedName)
  print(table(myDataset$class))
  print(length(myDataset$keys))
}


# Create meta-object and run meta-analysis
agingMetaObj = list(originalData = agingLadies)
agingMetaObj = runMetaAnalysis(agingMetaObj, runLeaveOneOutAnalysis = F)


# 6 genes pass FDR < 5%
# Higher in younger ladies:
#   - Confusing: RPAP2
#   - NK/T: ZNF827
# Higher in older ladies:
#    Granulocyte: CTSG, 
#    Mono/macro: ALDH2, LGALS1, ABCA13
agingMetaObj = filterGenes(agingMetaObj, isLeaveOneOut = F, FDRThresh = 0.05, numberStudiesThresh = 2)
summarizeFilterResults(agingMetaObj, getMostRecentFilter(agingMetaObj))

#..............#
##### Save #####
#..............#

save(sexMetaObj, valMetaObj,agingMetaObj,
     ISEXSLoc, isexs, xySig, autoSig, file = 'sexMetaObj.RData')



