# Preprocess Tsang cohort
# December 2017
# Data obtained from ImmPort (SDY80)
# Expression data available as gse47353


# Cohort purpose
# Studied vaccine responses in 63 people

### So far: 
# I got sex and cell proportions from ImmPort (SDY80)
# I got age annotations from GEO
# I mapped Immport IDs to GSMs

### Future: 
# Create a keys system 
#   - R friendly nomenclature
#   - HIMC nomenclature
#   - Consistent names across cohorts

#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(MetaIntegrator)
library(dplyr)
library(reshape2)
library(data.table)

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/4_sortedCells/cellProp/SDY80/")

#..................#
##### gse47353 #####
#..................#
gse47353 = getGEOData("GSE47353")
gse47353 = gse47353$originalData$GSE47353

# Clean up pheno
gse47353$rawPheno = gse47353$pheno
gse47353$pheno = cleanUpPheno(gse47353$rawPheno, T)

# Remove duplicate columns
gse47353$pheno = gse47353$pheno[,!grepl(pattern = ":ch1", x = colnames(gse47353$pheno))]

# Examine age
# Large spike a college age
hist(gse47353$pheno$age)

# Examine sex
# 170 females, 122 males
table(gse47353$pheno$sex)

# Fix timepoint
timepoint= as.character(gse47353$pheno$sample_collection_time)
timepoint = strsplitVector(timepoint, '\\(day', 2)
timepoint = strsplitVector(timepoint, "\\)", 1)
gse47353$pheno$timepoint = as.numeric(timepoint)

### Check keys
gse47353$pheno$platform_id[1]
sum(is.na(gse47353$keys))/length(gse47353$keys) # 27% are NAs, that's normal
head(gse47353$keys[!is.na(gse47353$keys)]) # Looks good

### Check expr 

# No obvious batch effect
boxplot(gse47353$expr, main = "gse47353")

# No negative values
range(gse47353$expr, na.rm = T) # 1.6 - 13.7


### Check sex labels

"XIST" %in% gse47353$keys # FALSE
"RPS4Y1" %in% gse47353$keys # TRUE
"KDM5D" %in% gse47353$keys # TRUE

# Looks like 5 samples have incorrect sex labels
quickViolin(gse47353, "RPS4Y1", "sex") # 3 females, 2 males look wrong
quickViolin(gse47353, "KDM5D", "sex") # 3 females, 2 males look wrong

# Find GSMs of mislabeled samples
imputedSex = imputeSex(gse47353)
sum(imputedSex != gse47353$pheno$sex) # 5
wrongSex = rownames(gse47353$pheno)[imputedSex != gse47353$pheno$sex]

# Does baseline pheno have mislabeled samples?
phenoBaseline = subset(gse47353$pheno, sample_collection_time == "immediately before vaccination (day0)")
dim(phenoBaseline) # 57 rows
sum(wrongSex %in% rownames(phenoBaseline)) # 2 of the mislabeled samples are in baseline

# Since we need all the GSMs for mapping, I won't remove samples yet

checkDataObject(gse47353, "Dataset") # True

#.............................#
##### Mapping ImmPort IDs #####
#.............................#
# Purpose: 
#   GEO has age, but ImmPort doesn't, the assholes
#   So, Map age from gse47353 onto subjectInfo 

# Subjects labeled by SUBJECT_ACCESSION
subjectInfo = read.delim("subject.txt", stringsAsFactors = F)

# SUBJECT_ACCESSION
# BIOSAMPLE_ACCESSION
biosample = read.delim("biosample.txt", stringsAsFactors = F)

# BIOSAMPLE_ACCESSION
# EXPSAMPLE_ACCESSION
exp_2_bio = read.delim("expsample_2_biosample.txt", stringsAsFactors = F)

# EXPSAMPLE_ACCESSION
# REPOSITORY_ACCESSION (GEO)
exp_2_geo = read.delim("expsample_public_repository.txt", stringsAsFactors = F)


### 1) Add a biosample ID column ot exp_2_geo
all(exp_2_geo$EXPSAMPLE_ACCESSION %in% exp_2_bio$EXPSAMPLE_ACCESSION) # True
bioDICTexp = exp_2_bio$BIOSAMPLE_ACCESSION
names(bioDICTexp) = exp_2_bio$EXPSAMPLE_ACCESSION
exp_2_geo$biosample = bioDICTexp[exp_2_geo$EXPSAMPLE_ACCESSION]


### 2) Add a Subject Accesion to exp_2_geo
all(exp_2_geo$biosample %in% biosample$BIOSAMPLE_ACCESSION) # True
subjDICTbio = biosample$SUBJECT_ACCESSION
names(subjDICTbio) = biosample$BIOSAMPLE_ACCESSION
exp_2_geo$subjectID = subjDICTbio[exp_2_geo$biosample]


### 3) Add sex from ImmPort's subject info
all(exp_2_geo$subjectID %in% subjectInfo$SUBJECT_ACCESSION) # True
sexDICTsubj = subjectInfo$GENDER
names(sexDICTsubj) = subjectInfo$SUBJECT_ACCESSION
exp_2_geo$immPortSex = tolower(sexDICTsubj[exp_2_geo$subjectID])

### 4) Add sex from GEO
all(exp_2_geo$REPOSITORY_ACCESSION %in% rownames(gse47353$pheno)) # True
sexDICTgsm = as.character(gse47353$pheno$sex)
names(sexDICTgsm) = rownames(gse47353$pheno)
exp_2_geo$geoSex = sexDICTgsm[exp_2_geo$REPOSITORY_ACCESSION]
table(exp_2_geo$immPortSex, exp_2_geo$geoSex) # 14 samples disagree

### 5) Add age from GEO
ageDICTgsm =gse47353$pheno$age
names(ageDICTgsm) = rownames(gse47353$pheno)
exp_2_geo$age = ageDICTgsm[exp_2_geo$REPOSITORY_ACCESSION]


### 6) Explore samples that don't add up
missmatch = exp_2_geo[exp_2_geo$immPortSex != exp_2_geo$geoSex,]
table(missmatch$subjectID, missmatch$age)

### 7) create final sex column
sex = exp_2_geo$geoSex
sex[sex != exp_2_geo$immPortSex] = NA
exp_2_geo$sex = sex

### 8) Add timepoint from GEO
timepoint = gse47353$pheno$timepoint
names(timepoint) = rownames(gse47353$pheno)
all(exp_2_geo$REPOSITORY_ACCESSION %in% names(timepoint)) # true
exp_2_geo$timepoint = timepoint[as.character(exp_2_geo$REPOSITORY_ACCESSION)]

#...................#
##### Cell Prop #####
#...................#

# Load cell proportions data
cellProp = read.delim("fcs_analyzed_result.txt")

# Remove columns that are entirely NA
allNA = apply(cellProp, 2, function(x) all(is.na(x)))
cellProp = cellProp[,!allNA]

# Remove columns that are all the same
allSame = apply(cellProp,2 ,function(x) length(unique(x))== 1)
cellProp= cellProp[,!allSame]

# Subset down to timepoint 0
#cellProp = subset(cellProp, STUDY_TIME_COLLECTED == 0)
sexDict =  exp_2_geo$sex

# Create the expr matrix
# expr = cellProp %>% reshape2::dcast(POPULATION_NAME_REPORTED ~ SUBJECT_ACCESSION, value.var = "POPULATION_STATISTIC_REPORTED")
# rownames(expr) = expr$POPULATION_NAME_REPORTED
# expr$POPULATION_NAME_REPORTED = NULL
# expr = as.matrix(expr)

#...............#
##### Pheno #####
#...............#
# I no longer want to make a Dataset object with cell prop
# So, this isn't really important

# Create phenotype matrix that matches expr
# pheno = exp_2_geo
# pheno = subset(pheno, subjectID %in% colnames(expr))
# pheno = pheno[,c("subjectID", "immPortSex", "geoSex", "sex", "age")]
# pheno = unique(pheno)
# 
# length(pheno$subjectID) == length(unique(pheno$subjectID)) # True!
# rownames(pheno) = as.character(pheno$subjectID)
# 
# # True, they're all in each other!
# all(rownames(pheno) %in% colnames(expr))
# all(colnames(expr) %in% rownames(pheno))
# pheno = pheno[colnames(expr),]
#................................................#
##### Create What cell Subsets we care about #####
#................................................#
# Purpose:
#   - Extract cell types that correspond to gse65133
#
# gse65133 cell types:
#   CD4+ T cells 
#      - resting, active, memory, and naive
#   CD8+ T cells 
#   B cells
#      - memory, naive
#   Monocytes
#   gd T cells
#   NK cells


### Identify corresponding cell types
cellTypes = unique(as.character(cellProp$POPULATION_NAME_REPORTED)) # all 113 available cell types


# Extract major cell types
majorCellTypes = cellTypes[grep(pattern = '\\)$', x = cellTypes)]


# Create R-friendly names for the major cell types
majorCellTypes_names = majorCellTypes
majorCellTypes_names[5] =  "ID34, CD45RA+ of CD4+ T cells (Naive CD4+ T)"
majorCellTypes_names = strsplitVector(majorCellTypes_names, "\\(", 2)
majorCellTypes_names = strsplitVector(majorCellTypes_names, "\\)", 1)
majorCellTypes_names = gsub(pattern = ' ', replacement = '_', x = majorCellTypes_names)
majorCellTypes_names = gsub(pattern = "\\+_", replacement = "pos_", x = majorCellTypes_names)
majorCellTypes_names = gsub(pattern = "\\+", replacement = "pos_", x = majorCellTypes_names)
majorCellTypes_names = gsub(pattern = "\\-_", replacement = "neg_", x = majorCellTypes_names)
majorCellTypes_names = gsub(pattern = "\\-", replacement = "neg_", x = majorCellTypes_names)

# Create a dictionary to turn ugly names into R-friendly names
majorCellTypes_dict = majorCellTypes_names
names(majorCellTypes_dict) = majorCellTypes

#............................................................#
##### Connect the cell types we're interested in to GSMs #####
#............................................................#
# Purpose: 
#   - Subset cellProp to major cell types
#   - Give better names
#   - Connect Subject IDs and timepoints to GSMs

# All cell types presnt in pop name reported!
all(names(majorCellTypes_dict) %in% as.character(cellProp$POPULATION_NAME_REPORTED))

cellProp = cellProp[,c("POPULATION_NAME_REPORTED",
                       "POPULATION_STATISTIC_REPORTED",
                       "STUDY_TIME_COLLECTED",
                       "SUBJECT_ACCESSION")]

# Add column with major cell type name
cellProp = subset(cellProp, POPULATION_NAME_REPORTED %in% names(majorCellTypes_dict))
cellProp$cellType = majorCellTypes_dict[as.character(cellProp$POPULATION_NAME_REPORTED)]

# Remove ugly cell type name
cellProp$POPULATION_NAME_REPORTED = NULL

# Add subjectKey
cellProp$subjectKey = paste(cellProp$SUBJECT_ACCESSION, cellProp$STUDY_TIME_COLLECTED, sep="_")
exp_2_geo$subjectKey = paste(exp_2_geo$subjectID, exp_2_geo$timepoint, sep='_')

subjectKey_dict = exp_2_geo$REPOSITORY_ACCESSION
names(subjectKey_dict) = exp_2_geo$subjectKey

cellProp$GSM = subjectKey_dict[as.character(cellProp$subjectKey)]
View(cellProp)

# Get Dataframe where rownames are GSMs and columns are cell types
cellProp = subset(cellProp, !is.na(GSM))
cellProp2 = cellProp %>% reshape2::dcast(GSM ~ cellType, value.var = "POPULATION_STATISTIC_REPORTED")
rownames(cellProp2) = cellProp2$GSM
cellProp2$GSM = NULL

cellProp2 = cellProp2[rownames(gse47353$pheno),]
rownames(cellProp2) = rownames(gse47353$pheno)

# Add cell types to pheno
all(rownames(cellProp2) == rownames(gse47353$pheno)) # true
gse47353$pheno = cbind(gse47353$pheno, cellProp2)
checkDataObject(gse47353, "Dataset") # T

#...................................................#
##### Sanity check gene expression vs celltypes #####
#...................................................#
# Purpose: 
#   - Correlate well known cell type gene combinations
#   - Make sure they generally correlate

### CD19 and B cells --- Messy, but significant correlation
"CD19" %in% gse47353$keys
bCells = gse47353$pheno$Total_B_cells
cd19 = unlist(getSampleLevelGeneData(gse47353, "CD19"))
length(bCells)
length(cd19)
cor.test(bCells, cd19)
plot(bCells, cd19)

### CD14 and monocytes --- Messy, but signififcant correlation
'CD14' %in% gse47353$keys
mono = gse47353$pheno$Total_Monocytes
cd14 = unlist(getSampleLevelGeneData(gse47353, "CD14"))
cor.test(mono, cd14) # r = 0.53
plot(mono, cd14)

### CD3D and CD4+T cells --- Much messier, but significant correlation
'CD3D' %in% gse47353$keys
cd4 = gse47353$pheno$Total_T_cells
cd3d = unlist(getSampleLevelGeneData(gse47353, 'CD3D'))
cor.test(cd4, cd3d)
plot(cd4, cd3d)

#............................................#
##### Remove Samples with mislabeled Sex #####
#............................................#
# Purpose: 
#   - Remove samples whose labeled sex conflicts with
#     gene expression


"XIST" %in% gse47353$keys #F
"RPS4Y1" %in% gse47353$keys # T
"KDM5D" %in% gse47353$keys # T

# 2 males and 3 females mislabeled
quickViolin(gse47353, "RPS4Y1", "sex")
quickViolin(gse47353, "KDM5D", "sex")
quickScatter(gse47353, "KDM5D", 'RPS4Y1', 'sex')

# Impute sex
imputedSex = imputeSex(gse47353)

# Go thorugh each mislabeled sample and remove it
mislabeledSex = rownames(gse47353$pheno)[gse47353$pheno$sex != imputedSex]
for(mySample in mislabeledSex){
  gse47353 = removeOneSample(gse47353, mySample)
}

# No more weirdos!
quickViolin(gse47353, "RPS4Y1", "sex")

#............................#
##### Subset to baseline #####
#............................#
# Purpose: 
#   - Subset to baseline

# Use -7 because it has more samples
table(gse47353$pheno$timepoint, gse47353$pheno$sex)

# Create baseline pheno
pheno7 = subset(gse47353$pheno, timepoint == '-7')
dim(pheno7)
unique(pheno7$timepoint) # only -7

# Create baseline dataset object
gse47353_baseline = subsetGEMFromPheno(gse47353, pheno7)

#..............#
##### Save #####
#..............#
# Save final dataset objects

save(gse47353, gse47353_baseline, file ="gse47353.RData")
