# Download and Pre-process GSE23350

# Conclusions: 
# Only 1 female and 4 males >20 years old
# Only 14 females and 5 males <20 years old
# Interesting trends with CD40LG and MPO
# Only useful if looking at puberty, I think 


#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/4_sortedCells/cd4_Tcells/gse23350/")

#..................#
##### Download #####
#..................#

gse23350 = getGEOData("GSE23350")
gse23350 = gse23350$originalData$GSE23350

gse23350$rawPheno = gse23350$pheno
gse23350$pheno = cleanUpPheno(gse23350$rawPheno, T)

# Fix age
age = as.character(gse23350$pheno$age_at_sampling)
age = as.numeric(strsplitVector(age, "y", 1))
gse23350$pheno$age = age

table(gse23350_cont$pheno$sex, gse23350_cont$pheno$age < 20)

#..........................#
##### Try out Controls #####
#..........................#
phenoCont = subset(gse23350$pheno, source_name_ch1 == "healthy control")
gse23350_cont = subsetGEMFromPheno(gse23350, phenoCont)


quickTimecourse(gse23350_cont, myGene = "CD40LG", timeColumn = "age", groupColumn = "sex", timeUnits = "years")
