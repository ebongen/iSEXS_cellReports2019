# Download and pre-process GSE21862

## Background
# Dataset on benezene exposure in factory workers in China
# Benezene exposure can lead to cancer
# Their model included chip run and hybridization in order 
# to take batch effect into account. So, they included that info in pheno

## Dataset attributes
# Benezene exposure (control, very low, low, high, very high)
# Smoking status yes/no

## Methods
# expr looks good, left as is
# removed technical replicates. Only 19 samples had them

#...................#
##### Load Stuff ####
#...................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/0_discovery/gse21862/")
#..................#
##### Download #####
#..................#
gse21862 = getGEOData("GSE21862")
gse21862 = gse21862$originalData$GSE21862

#...............#
##### Pheno #####
#...............#
# Clean up pheno
gse21862$rawPheno = gse21862$pheno
gse21862$pheno = cleanUpPheno(gse21862$pheno, T)

# Age
gse21862$pheno$age # looks good

# Sex
gse21862$pheno$sex # looks good

# Technical replicates and SubjectID
theTitle = as.character(gse21862$pheno$title)

subjectID = gsub(pattern = "_technical rep 1", replacement = "", x = theTitle)
subjectID = gsub(pattern = "_technical rep 2", replacement = "", x = subjectID)

techRep = substring(text = theTitle, first = 14)
techRep[which(techRep == "")] = "technical rep 1"

subjectHasReplicat = subjectID[which(techRep=="technical rep 2")]
hasReplicate = subjectID %in% subjectHasReplicat

# Make sure my subjectID looks right compared to the original title
test = cbind(as.character(gse21862$pheno$title), subjectID, techRep, hasReplicate)

gse21862$pheno$subjectID = as.factor(subjectID)
gse21862$pheno$technicalReplicate = as.factor(techRep)



#..............#
##### Expr #####
#..............#
# Batch effect? Nope!
# Log2? Yup!
# Expected range? Yup!
# Above zero? Yup!
boxplot(gse21862$expr, main = "gse21862")

range(gse21862$expr) # [3.68, 16.82]

#......................#
##### Class + Keys #####
#......................#
# Create class vector
gse21862$class = createClassVector(gse21862$pheno$sex, "female", gse21862$pheno)

# Sanity check keys, looking fine
length(gse21862$keys) # 22177 probes total
length(na.omit(gse21862$keys)) # 17999 non-NA probes
length(unique(na.omit(gse21862$keys))) # 14512 unique genes

# Formatted name!
gse21862$formattedName = "GSE21862_China_PBMC"

#.................................#
##### Divvy Discovery Dataset #####
#.................................#
# Remove individuals older than 40
phenoYoung = subset(gse21862$pheno, age <= 40)

# Remove individuals with high benzene exposure
# I'm only taking people with benzene <1ppm or less
# <1ppm is the US occupational hazard standard
phenoYoung = subset(phenoYoung, benzene != ">10ppm")
phenoYoung = subset(phenoYoung, benzene != "5-10ppm")

# Remove the second technical replicates
phenoYoung = subset(phenoYoung, technicalReplicate =="technical rep 1")

# Sanity check
table(phenoYoung$age) # between 18 and 39, mostly in 20s
table(phenoYoung$sex) # 43 female, 44 male
table(phenoYoung$benzene) # only control, <<1ppm, and <1ppm

gse21862_discovery = subsetGEMFromPheno(gse21862, phenoYoung)

#..........................#
##### Check Sex Labels #####
#..........................#
# Purpose: Make sure the sex labels match expresion of known genes
# Results: All discovery samples have correct sex labels, according to RPS4Y1 expression


"XIST" %in% gse21862_discovery$keys # False
"RPS4Y1" %in% gse21862_discovery$keys # TRUE
"KDM5D" %in% gse21862_discovery$keys # FALSE

# Separates clearly by sex
quickViolin(gse21862_discovery, "RPS4Y1", "sex")


imputedSex = imputeSex(gse21862_discovery)
all(imputedSex == gse21862_discovery$pheno$sex) # Sex labels completely match!


#..............#
##### Save #####
#..............#
checkDataObject(gse21862, "Dataset") # True
checkDataObject(gse21862_discovery, "Dataset") # True

save(gse21862, gse21862_discovery, file = "gse21862_china.RData")
