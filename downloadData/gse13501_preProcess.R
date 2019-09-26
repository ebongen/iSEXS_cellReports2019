# August 2017 PreProcess JIA dataset gse13501



#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")
source("00_tools/plot_funx.R")

load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/2_age/prePuberty/gse13501_ohio_JIA/")

#..................#
##### Download #####
#..................#
gse13501 = getGEOData("GSE13501")
gse13501 = gse13501$originalData$GSE13501

gse13501$rawPheno = gse13501$pheno
gse13501$pheno = cleanUpPheno(gse13501$rawPheno, F)

# 59 controls
table(gse13501$pheno$source_name_ch1)
table(gse13501$pheno$source_name_ch1, )
##### 