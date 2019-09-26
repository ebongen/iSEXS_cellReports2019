# 6/20/2018 - HIV resistant ladies

# Purpose: 
#  - 

#....................#
##### Load Stuff #####
#....................#
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Load packages
library(MetaIntegrator)
library(ggplot2)
library(cowplot)

# Load some of my code
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

setwd("0_datasets/6_infection/HIV/gse33580_hivResistance/")

#............................#
##### Download and Pheno #####
#............................#
# Download
gse33580 = getGEOData('GSE33580')
gse33580 = gse33580$originalData$GSE33580

# Clean up pheno
gse33580$rawPheno = gse33580$pheno
gse33580$pheno = cleanUpPheno(gse33580$rawPheno, T)

# Add subject ID 
subject_id = strsplitVector(gse33580$pheno$title, '_', 2)
head(subject_id)
length(unique(subject_id)) == length(subject_id) # True!
gse33580$pheno$subject_id = subject_id

# Add sex column 
gse33580$pheno$sex = rep('female', nrow(gse33580$pheno))

#..............#
##### Expr #####
#..............#
# Purpose: 
#   - Make sure expr is up to snuff
# Results:
#   - No NAs
#   - Contains negative values

# No NAs
range(gse33580$expr) # -7 to 16

# Some values are below zero
# some mesiness, need to normalize quantiles
boxplot(gse33580$expr, main ='GSE33580')

# Raise above zero
gse33580$expr = gse33580$expr + abs(min(gse33580$expr)) + 1.01

# Quantile normalize
test = normalize.quantiles(gse33580$expr)
dim(test)
rownames(test) = rownames(gse33580$expr)
colnames(test) = colnames(gse33580$expr)
boxplot(test, main='Quantile Normalized')

gse33580$rawExpr = gse33580$expr
gse33580$expr = test

checkDataObject(gse33580, 'Dataset') # True

### Check sex labels
quickScatter(gse33580, 'XIST', 'RPS4Y1', 'hiv_status')

quickViolin(gse33580, 'XIST', 'hiv_status')
xist = unlist(getSampleLevelGeneData(gse33580, 'XIST'))
t.test(xist[gse33580$pheno$hiv_status == 'HIV negative'],
       xist[gse33580$pheno$hiv_status == 'HIV resistant'])

#..............#
##### Save #####
#..............#
save(gse33580, file='gse33580.RData')
