# 8/16/2017 Download and Pre-process gse59654
# HIPC Young vs Old from Yale

# Note:
# I had to impute sex
# Hopefully I'll get sex labels later from Francesco

#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
library(ggplot2)

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")


setwd("0_datasets/2_age/postMenopause/gse59654/")

#..................#
##### Download #####
#..................#

gse59654 = getGEOData("GSE59654")
gse59654 = gse59654$originalData$GSE59654

# Tidy up Pheno
gse59654$rawPheno = gse59654$pheno
gse59654$pheno = cleanUpPheno(gse59654$rawPheno, T)

View(gse59654$pheno)

# Check expr
# No obvious batch effect
# range between 6 and 14
# Squished boxplots more than I'd like, but probably okay
gse59654$pheno$platform_id[1] # Illumina gpl10558
boxplot(gse59654$expr, main = "gse59654")

# Check keys
# Looks good!
print(gse59654$keys[1:20])

#....................#
##### Impute Sex #####
#....................#


"XIST" %in% gse59654$keys # False
"RPS4Y1" %in% gse59654$keys # True
"KDM5D" %in% gse59654$keys # False


imputedSex = imputeSex(gse59654)
table(imputedSex) # 88 female 68 male

# Create sex column in pheno, and mark that it's imputed
gse59654$pheno$sex = imputedSex
gse59654$pheno$imputedSex = rep(TRUE, nrow(gse59654$pheno))

table(gse59654$pheno$sex, gse59654$pheno$age_group)


# Create the class vector based on imputed sex
gse59654$class = createClassVector(gse59654$pheno$sex, "female", gse59654$pheno)



#............................................#
###### Divvy into old and young datasets #####
#............................................#

# Create baseline
phenoBase = subset(gse59654$pheno, blood_draw_date == 0)
gse59654_baseline = subsetGEMFromPheno(gse59654, phenoBase)

# Create young baseline
phenoYoung = subset(phenoBase, age_group == "Young")
gse59654_young = subsetGEMFromPheno(gse59654, phenoYoung)

# Create allOld baseline
phenoAllOld = subset(phenoBase, age_group !="Young")
unique(phenoAllOld$age_group) # Frail and Old
gse59654_frailOld = subsetGEMFromPheno(gse59654, phenoAllOld)

# Create Old only baseline
phenoOld = subset(phenoBase, age_group == "Old")
gse59654_old = subsetGEMFromPheno(gse59654, phenoOld)

# Creaet Frail only baseline
phenoFrail = subset(phenoBase, age_group == "Frail")
gse59654_frail = subsetGEMFromPheno(gse59654, phenoFrail)

#..............#
##### Save #####
#..............#
save(gse59654, gse59654_baseline, 
     gse59654_frail, gse59654_old, gse59654_frailOld, 
     gse59654_young, file = "gse59654_aging.RData")
