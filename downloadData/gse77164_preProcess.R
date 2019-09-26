# 8/16/2017 Preprocess gse77164

# Results: 
# Youngest person is 15, so doesn't cover puberty
# well enough to answer my question
# Excellent coverage of 18-22yo
#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")
source("00_tools/plot_funx.R")

setwd("0_datasets/2_age/prePuberty/gse77164_nepal/")

#..................#
##### Download #####
#..................#

gse77164 = getGEOData("GSE77164")
gse77164 = gse77164$originalData$GSE77164

gse77164$rawPheno = gse77164$pheno
gse77164$pheno = cleanUpPheno(gse77164$rawPheno, T)

View(gse77164$pheno)
hist(gse77164$pheno$age) # Most of the samples are between 18-22

# Check expr
boxplot(gse77164$expr[,1:5])
min(gse77164$expr) # positive value

# Check keys
# They're fine
gse77164$keys[1:10]

#..........................#
##### Quick Timecourse #####
#..........................#
load("/labs/khatrilab/ebongen/sexDifferences2.0/1_metaAnaly/sexMetaObj.RData")

pdf("gse77164_timecourse_blueFemales.pdf")
timecourse_iSEXS(gse77164, sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0, "age", "female", main = "GSE77164 - full iSEXS")
timecourse_iSEXS(gse77164, sexMetaObj$filterResults$autosomeOnly, "age", "female", main = "GSE77164 - Autosomal iSEXS")
dev.off()

