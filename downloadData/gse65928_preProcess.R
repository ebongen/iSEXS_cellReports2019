# December 2017 - Download and preprocess gse65928

# Background:
# I want to know if a subset of ISEXS is diff expressed in B cells
# gse65928 has B cell subsets from healthy Malian adults
# Malaria is endemic there, so they have atypical B cells, naive B cells, and classic B cells
# I"ll create a Dataset object for each

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

# Necessary packages
library(MetaIntegrator)
library(biomaRt) # converting RefSeq IDs to gene symbols

# Useful scripts
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

# Load ISEXS
load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/4_sortedCells/b_cells/gse65928/")

#..................#
##### Download #####
#..................#
# Malaria paper about abnormal B cell type
gse65928 = getGEOData("GSE65928")
gse65928 = gse65928$originalData$GSE65928

gse65928$rawPheno = gse65928$pheno
gse65928$pheno = cleanUpPheno(gse65928$rawPheno, T)
View(gse65928$pheno)

# add cell type annotation
cellType = gse65928$pheno$source_name_ch1
cellType = strsplitVector(cellType, "Human ", 2)
gse65928$pheno$cellType = cellType

#..................................#
##### Manually Annotate Probes #####
#..................................#
# Purpose: 
# All the keys are labeled as NA
# Automatic key label has failed
# So, I'll have to do it manually by identifying the label, 
# and converting it to gene symbols

### None of them are in here!!!
"XIST" %in% gse65928$keys
"RSP4Y1" %in% gse65928$keys
"KDM5D" %in% gse65928$keys


### Keys - All the keys are NA!!!!!
length(gse65928$keys) # 53617
gse65928$keys[1:100]
all(is.na(gse65928$keys)) # True

# Get the gpl ID of this platform
platformID = as.character(gse65928$pheno$platform_id[1])
myGpl = GEOquery::getGEO(platformID)

# Get list of all GenBank IDs referenced by this platform
myGenbank = as.character(myGpl@dataTable@table$GB_ACC)
myGenbank = myGenbank[myGenbank != ""] # remove blank entries

# Load ensembl database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# This gpl has labels from RefSeq for non-coding RNA (ncrna) and mRNA (mrna)
dat_NR = getBM(attributes = c("hgnc_symbol", "refseq_ncrna"), values = "*", mart = ensembl)
dat_NM = getBM(attributes = c("hgnc_symbol","refseq_mrna"), values = "*", mart = ensembl)
dat_XM = getBM(attributes = c("hgnc_symbol", "refseq_mrna_predicted"), values = "*", mart = ensembl)
dat_XR = getBM(attributes = c("hgnc_symbol", "refseq_ncrna_predicted"), values = "*", mart = ensembl)

# Combine all four annotations together
masterCode = rbind(as.matrix(dat_NM), as.matrix(dat_XM), as.matrix(dat_NR), as.matrix(dat_XR))
colnames(masterCode) = c("geneSymbol", "RefSeq")
masterCode = masterCode[masterCode[,2] != "",]
masterCode = masterCode[masterCode[,2]%in% myGenbank,]

# Create a keys system that converts RefSeq IDs to Gene Symbols
myKeys = masterCode[,1]
names(myKeys) = masterCode[,2]
myKeys = myKeys

# Add gene symbol to gpl annotation
myGpl@dataTable@table$geneSymbol = myKeys[as.character(myGpl@dataTable@table$GB_ACC)]

# Create a keys system that converts probe IDs to Gene Symbols
myKeys2 = myGpl@dataTable@table$geneSymbol
names(myKeys2) = myGpl@dataTable@table$ID
gse65928$keys = myKeys2[rownames(gse65928$expr)]

checkDataObject(gse65928, "Dataset")
# True!

#....................#
##### Impute Sex #####
#....................#
# Purpose: 
# Impute sex using expression of XIST, RPS4Y1, and KDM5D

### Impute sex
imputedSex = imputeSex(gse65928)
gse65928$pheno$sex = imputedSex
gse65928$pheno$imputedSex = rep(TRUE, nrow(gse65928$pheno))

# Very clear separation
quickViolin(gse65928, "XIST", "sex")
quickViolin(gse65928, "KDM5D", "sex")

gse65928$class = createClassVector(gse65928$pheno$sex, "female", gse65928$pheno)

#.............................#
##### Subset by cell type #####
#.............................#
# Purpose:
# Create separate Dataset objects for each B cell subtype in this cohort

# Naive B cells
phenoNaive = subset(gse65928$pheno, cellType == "naive B-cell")
gse65928_naive = subsetGEMFromPheno(gse65928, phenoNaive)

# Classical B cells
phenoClassic = subset(gse65928$pheno, cellType == "classical B-cell")
gse65928_classic = subsetGEMFromPheno(gse65928, phenoClassic)

# Atypical B cells associated with chronic malaria infection 
phenoAtyp = subset(gse65928$pheno, cellType == "atypical B-cell")
gse65928_atyp = subsetGEMFromPheno(gse65928, phenoAtyp)

#..............#
##### Save #####
#..............#

save(gse65928, gse65928_atyp, gse65928_classic, gse65928_naive, file = "gse65928.RData")
