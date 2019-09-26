# 6/21/2018 - gse48023 Females give the flu vaccine


# Purpose:
#   - I want to see how auto-iSEXS changes with vaccination
#   - I also want to know if auto-iSEXS can predict vaccine response

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

# Libraries
library(MetaIntegrator)
library(ggplot2)
library(cowplot)

# Code
source("00_tools/general_GEMfunx.R")
source('~/Tools/Graphing Scripts/quickPlots.R')

# iSEXS
load('1_metaAnaly/sexMetaObj.RData')

setwd("0_datasets/6_infection/vaccine/gse48023_females/")

#............................#
##### Download and Pheno #####
#............................#

gse48023 = getGEOData('GSE48023')
gse48023 = gse48023$originalData$GSE48023

# Tidy up pheno
gse48023$rawPheno = gse48023$pheno
gse48023$pheno = cleanUpPheno(gse48023$pheno, T)

# Make sure timepoint is nice
# Make time into a factor
str(gse48023$pheno$time)
unique(gse48023$pheno$time)
time = factor(gse48023$pheno$time, levels = unique(gse48023$pheno$time)[c(1,2,4,3)])
table(gse48023$pheno$time, time) # All match!
gse48023$pheno$time = time

#..............#
##### Expr #####
#..............#
# Purpose: 
#   - Make sure expr is preprocessed right
#
# Results:
#   - No NAs
#   - All positive
#   - No obvious batch effect

# 5 to 16
range(gse48023$expr)

# No obvious batch effect
pdf('expr.pdf', width=14, height = 7)
boxplot(gse48023$expr, main = 'gse48023')
dev.off()

#............................#
##### Confirm sex labels #####
#............................#
# Purpose: 
#   - Make sure they really are all female
# Results:
#   - No major obvious outliers
#   - Cluster where you expect females
#   - RPS4Y1 expression not as low as I'd like
#
# Conclusions:
#   - All female

'XIST' %in% gse48023$keys #T
'RPS4Y1' %in% gse48023$keys # T

quickScatter(gse48023, 'XIST', 'RPS4Y1', 'time')

# 1Q: 6.8
# Median: 6.9
# 3Q: 7.3
summary(gse48023$expr[,1])

#...........................#
##### Add in HAI Titers #####
#...........................#
# Purpose: 
#   - Add HAI titers to pheno
#   - Add in change in titers (Day 28-Day 0) for each Ag tested
#   - Consider taking median/mean/max of across the antigens
