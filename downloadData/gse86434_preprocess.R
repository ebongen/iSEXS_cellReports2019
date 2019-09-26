# December 2017 - Preprocess gse86434

# Background:
# gse86434 has sorted cells from UC, CD, and HC

# Conclusion: 
# While gse86434 has very valueble naive CD4 T cell expression
# The actual numbers of HC is actually very small
# Only 5 CD4 young males and 2 CD4 young females
# So, it's not worth the hassle that it takes to download it

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
library(GEOquery)
source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/4_sortedCells/cd4_Tcells/gse86434/")

#..................#
##### Download #####
#..................#
# gse86434 is formatted poorly, so it has to be manually downloaded and reformatted



#...............#
##### Pheno #####
#...............#
pheno = read.delim(file = "gse86434_pheno.txt")
pheno = t(pheno)

### Make unique column names
myCol = pheno[1,]
for (i in 1:length(myCol)){
  myCol[i] = str_replace_all(myCol[i], "!", "")
  if(grepl(pattern = "characteristics_ch1", x = myCol[i])){
    myCol[i] = strsplit(x = pheno[5,i], split = ": ")[[1]][[1]]
  }
}

# Manually change remaining duplicate names
myCol[27] ="sample_number"
myCol[28] = "sample_id"
myCol[18] = "samplesection_letter"
myCol[19] = "samplesection_number"

colnames(pheno) = myCol
pheno = pheno[-1,]


### Remove useless columns
myCol = myCol[!grepl(pattern = "contact_", x = myCol)]

uselessColumns = c("Sample_status", "Sample_submission_date", "Sample_last_update_date", 
                   "Sample_type", "Sample_channel_count", "Sample_molecule_ch1", 
                   "Sample_extract_protocol_ch1", "Sample_label_ch1", "Sample_label_protocol_ch1", 
                   "Sample_taxid_ch1", "Sample_hyb_protocol", "Sample_scan_protocol", 
                   "Sample_data_processing", "Sample_supplementary_file")

all(uselessColumns %in% myCol)

myCol = myCol[!myCol %in% uselessColumns]
pheno = pheno[,myCol]
pheno = as.data.frame(pheno)

### Make columns pretty
sex = strsplitVector(pheno$Sex, ": ", 2)
sex = ifelse(sex == "F", yes = "female", no = "male")
pheno$Sex = NULL
pheno$sex = sex

age = as.numeric(strsplitVector(pheno$age, ": ", 2))
pheno$age = age

diagnosisFull = strsplitVector(pheno$`disease state (diagnosisfull)`, ": ", 2)
diagnosisFull = strsplitVector(diagnosisFull, "\\(", 2)
diagnosisFull = strsplitVector(diagnosisFull, "\\)", 1)
pheno$diagnosisFull = diagnosisFull
pheno$`disease state (diagnosisfull)` = NULL

diagnosisSimple = strsplitVector(pheno$`disease state (simplediagnosis)`, ": ", 2)
pheno$diagnosisSimple = diagnosisSimple
pheno$`disease state (simplediagnosis)` = NULL

smokingStatus = strsplitVector(pheno$`smoking status`, ": ", 2)
unique(smokingStatus)
pheno$smokingStatus = smokingStatus
pheno$`smoking status` = NULL

cellType = strsplitVector(pheno$`cell type`, ": ", 2)
pheno$cellType = cellType
pheno$`cell type` = NULL


### Make GSMs the rownames 
all(rownames(pheno) == pheno$sample_id) # Okay, so their sample IDs are the rownames right now
rownames(pheno) = pheno$Sample_geo_accession

#.............#
##### Expr ####
#.............#

expr = read.delim(file ="gse86434_expr.txt")
dim(expr)
rownames(expr) = expr$ID_REF
length(unique(rownames(expr))) == nrow(expr)
expr$ID_REF = NULL
expr = as.matrix(expr)


#..................................#
##### Explore specific subsets #####
#..................................#
phenoCD4 = subset(pheno, cellType == "CD4")
table(phenoCD4$diagnosisSimple, phenoCD4$sex)
table(phenoCD4$diagnosisFull, phenoCD4$sex)
range(phenoCD4$age) # 18-63

phenoCD4_healthy = subset(phenoCD4, diagnosisFull %in% c("Healthy control ", "Healthy control"))
table(phenoCD4_healthy$age < 41, phenoCD4_healthy$sex)
cbind(phenoCD4_healthy$age, phenoCD4_healthy$sex)

