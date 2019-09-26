# 3/19/2018 - Download GSE49454

# Purpose: 
#   - Chaussabel dataset with really cool annotations
#   - SLE, ANA titers, complement, age, etc

#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/5_autoImmunity/SLE/blood/gse49454/")

#..................#
##### Download #####
#..................#

# For some reason this didn't work
#gse49454 = MetaIntegrator::getGEOData("GSE49454")

#...............#
##### Pheno #####
#...............#
# Purpose: 
#   - preporcess the pheno from the series matrix file

# Grab pheno from series matrix file version
pheno = read.delim("gse49454_pheno.txt")
myRowNames = gsub(x = as.character(pheno$Sample_title), pattern = "!", replacement = "")
pheno$Sample_title = myRowNames

# Make rownames unique 
# aka make different names for characteristics_ch1
test = strsplitVector(pheno$C2PL28V2[pheno$Sample_title == "Sample_characteristics_ch1"], ": ", 1)
test = gsub(pattern = " ", replacement = "_", x = test)
test = gsub(pattern = "-", replacement = "_", x = test)
myRowNames[myRowNames == "Sample_characteristics_ch1"] = test
length(myRowNames) == length(unique(myRowNames)) # true
pheno$Sample_title = myRowNames

# Transpose so that samples are rows
pheno = t(pheno)
pheno = data.frame(pheno, stringsAsFactors = F)
pheno$sampleID = rownames(pheno)
colnames(pheno) = as.character(pheno[1,])
pheno = pheno[-1,]

# make GSMs the rownames
rownames(pheno) = pheno$Sample_geo_accession


# Remove useless columns
colnames(pheno)
pheno = pheno[, -c(2:6, 43:62)]

# Make characteristics columns useful
hasColon = apply(pheno, 2, function(x) any(grepl(pattern = ": ", x = x)))
table(hasColon) #34 columns with colons

# Identify scrambled columns
isScrambled = c()
for(i in 1:ncol(pheno)){
  # Get values of this column
  thisColumn = pheno[,i]
  
  # Remove empty cells
  thisColumn = thisColumn[thisColumn != ""]
  
  # If there are colons in this column
  if(any(grepl(pattern = ":", x = thisColumn))){
    # Check for scrambled labels (on left side of colons)
    if(length(unique(strsplitVector(thisColumn, ": ", 1))) > 1){
      isScrambled = c(isScrambled, i)
    }
  }
  
}

length(isScrambled) # 12 
colnames(pheno)[isScrambled]

