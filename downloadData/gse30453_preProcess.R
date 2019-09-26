# August 2017 Preprocessing gse30423

# Dataset background: 
# gse30423 uses an Affy exon array to figure out how
# SNPs regulate expression in different tissues
# They look at live healthy blood and autopsy brains
# I'm only using the GSE with the "gene-level" expression


#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/1_validation/gse30453_PBMC/")

#..................#
##### Download #####
#..................#

gse30453 = getGEOData("GSE30453")
gse30453 = gse30453$originalData$GSE30453

#...............#
##### Pheno #####
#...............#

# Clean up Pheno
gse30453$rawPheno = gse30453$pheno
gse30453$pheno = cleanUpPheno(gse30453$rawPheno, T)


# Check sex, looks good!
head(gse30453$pheno$sex)

# Chekc age, looks good!
head(gse30453$pheno$age)

# Create class vector
gse30453$class = createClassVector(gse30453$pheno$sex, "female", gse30453$pheno)

#....................#
##### Check Expr #####
#....................#
# Purpose: Check for batch effects or poor normalization

# Terrible batch effect, but they said they normalized brain and blood separate
# Not log2
boxplot(gse30453$expr, main = "gse30453 - Pre-log2")

# Minimum is 0, needs to be 1 to avoid probs with log2
min(gse30453$expr)
gse30453$expr = gse30453$expr - min(gse30453$expr) + 1


# Log2 normalize
gse30453$expr = log2(gse30453$expr)

# Looks much better, minimal match effect
# Some negative values
boxplot(gse30453$expr, main = "gse30453 - Post-log2")

# Minimum is 0, which is okay for downstream analyses
# Sometimes negative values messses stuff up
min(gse30453$expr)

#................................#
##### Subset into Validation #####
#................................#

# Create Dataset object of young blood
phenoVal = subset(gse30453$pheno, source_name_ch1 == "peripheral blood mononucleated cells collected from healthy controls")
phenoVal = subset(phenoVal, age <= 40)
range(phenoVal$age) # 18 to 40
table(phenoVal$sex) # 12 females, 29 males
gse30453_validation = subsetGEMFromPheno(gse30453, phenoVal)
gse30453_validation$formattedName = "GSE30453_validation"

checkDataObject(gse30453_validation, "Dataset") # True


# Create dataset object of old blood
phenoOld = subset(gse30453$pheno, source_name_ch1 == "peripheral blood mononucleated cells collected from healthy controls")
phenoOld = subset(phenoOld, age > 40)
range(phenoOld$age) # 41 to 67
table(phenoOld$age) # Doesn't get that far from my cutoff
table(phenoOld$sex) # 6 females 9 males
gse30453_oldBlood = subsetGEMFromPheno(gse30453, phenoOld)
checkDataObject(gse30453_oldBlood, "Dataset") # True


# Blood only dataset
phenoBlood = subset(gse30453$pheno, source_name_ch1 == "peripheral blood mononucleated cells collected from healthy controls")
gse30453_PBMC = subsetGEMFromPheno(gse30453, phenoBlood)
checkDataObject(gse30453_PBMC, "Dataset") # True

# Brains?
phenoBrain = subset(gse30453$pheno, source_name_ch1 != "peripheral blood mononucleated cells collected from healthy controls")
gse30453_brain = subsetGEMFromPheno(gse30453, phenoBrain)
checkDataObject(gse30453_brain, "Dataset") # True

#........................................#
##### Check Expr, again just in case #####
#........................................#
# Purpose: Check for batch effects or poor normalization


# check Validation dataset
# Looking good!
boxplot(gse30453_validation$expr, main = "gse30453 Validation")

# Old blood looks good!
boxplot(gse30453_oldBlood$expr, main = "old blood")

# Brain looks a bit messier, but good!
boxplot(gse30453_brain$expr, main = "brain")

#.......................................#
##### Check sex labels - Validation #####
#.......................................#
# All samples match their labels based on RPS4Y1 expression in validation

"XIST" %in% gse30453_validation$keys # False
"RPS4Y1" %in% gse30453_validation$keys # TRUE
"KDM5D" %in% gse30453_validation$keys # False

# Separates clearly by sex
# Doesn't look like there's any mislabeling
quickViolin(gse30453_validation, "RPS4Y1", "sex")

# It all matches, so no need to remove samples!
imputedSex = imputeSex(gse30453_validation)
all(imputedSex == gse30453_validation$pheno$sex) # True


#...................................#
##### Check sex labels - Others #####
#...................................#
# All of the subsets separate cleanly with no mislabels

quickViolin(gse30453_brain, "RPS4Y1", "sex") # Brain looks good
quickViolin(gse30453_oldBlood, "RPS4Y1", "sex") # Old blood looks good
quickViolin(gse30453_PBMC, "RPS4Y1", "sex") # PBMC looks good
quickViolin(gse30453, "RPS4Y1", "sex") # fulld dataset looks good

#..............#
##### Save #####
#..............#

save(gse30453_validation, file = "gse30453_validation.RData")
save(gse30453_validation, gse30453_oldBlood, gse30453_PBMC, file = "gse30453_agingBlood.RData")
save(gse30453, gse30453_brain, gse30453_oldBlood, gse30453_PBMC, gse30453_validation, file = "gse30453_full.RData")
