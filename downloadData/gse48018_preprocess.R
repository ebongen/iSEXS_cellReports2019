# 9/24/2018 - Preprocess gse48018


# Background:
#  - 119 males given flu vaccine

# Results: 
#   - I don't know where the vaccine reponse metrics are

#....................#
##### Load Stuff #####
#....................#
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Source functions
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

# Load packages
library(MetaIntegrator)
library(cowplot)

setwd("0_datasets/6_infection/vaccine/gse48018_males/")
#..................#
##### Download #####
#..................#
gse48018 = getGEOData('GSE48018')
gse48018 = gse48018$originalData$GSE48018


# Clean Pheno
gse48018$rawPheno = gse48018$pheno
gse48018$pheno = cleanUpPheno(gse48018$rawPheno, T)

# Sex
unique(gse48018$pheno$sex) # all male

# Time
unique(gse48018$pheno$time)
timepoint = factor(as.character(gse48018$pheno$time), levels = c('Day0', 'Day1','Day3', 'Day14'))
gse48018$pheno$timepoint = timepoint
