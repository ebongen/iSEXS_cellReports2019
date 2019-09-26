# October 2018 - Adapt BMI 212 code for sex differences project
# Author: Erika Bongen
# Contact: ebongen@stanford.edu
# Date: Started 5/19/2016, Adapted October 2018

# Purpose: 
# Preprocess and analyze differences in cell proportion from 
# Gabi's Cytof data

# Background: 
# Cell types
# CD45+CD66- = mononuclear cells
# 

## Gabi's cell explanation
# CD11b and CD14 were purely co-expressed, so any CD11b+ monocyte is also CD14+
# Unclear cell types: 
# [1] "CD7+/HLA-DR-" NK Cells (these are CD3-CD19-but CD7+, so these are NK cells and contain both of the CD56hi and lo populations)
# [1] "CD16+/CD11b-" CD16+ monocytes
# [1] "CD11b+/CD16-" Classical monocytes (CD11b and CD14 were purely co-expressed in the human data, so this is equivalent to CD14+CD16-)
# [1] "CD11b-/CD16-"
# [1] "CD16+/CD11b+" CD14+CD16+ monocytes


#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

# Packages
library(reshape2)
library(ggplot2)
library(cowplot)

# Source useful functions
source('00_tools/general_GEMfunx.R')


setwd('0_datasets/4_sortedCells/cellProp/gabi/')

# Phenotypic information for each donor
pheno = read.csv(file = "gabi_cytofdemographics.csv", header = T)
colnames(pheno) = tolower(colnames(pheno))


# Raw cell counts from each sample
rawCounts = read.csv(file="gabi_raw_cell_counts.csv", header = T)



#.......................#
##### Explore Pheno #####
#.......................#
# Get informationa bout cohort
summary(na.omit(pheno[,"age"])) # range 19-63
length(unique(pheno$donor)) # 93 people
summary(na.omit(pheno[, "gender"])) # 42 females, 49 males, but one female has NA for age
table(pheno$age <=40, pheno$gender)

#...........................#
##### Explore rawCounts #####
#...........................#
length(unique(rawCounts$File)) # 86 different files


#............................#
##### Re-name cell types #####
#............................#

cellTypes = as.character(unique(rawCounts$Population))
names(cellTypes) = gsub(pattern = " ", replacement = "_", x = cellTypes) # remove spaces
names(cellTypes) = gsub(pattern = "\\+", replacement = "", x = names(cellTypes)) # remove pluses

names(cellTypes)[which(cellTypes == "CD45+ CD66-")] = "mononuclearCells"
names(cellTypes)[which(cellTypes == "CD7+/HLA-DR-")] = "NK_cells"
names(cellTypes)[which(cellTypes == "CD4+/CD8+")] = "CD4_CD8_Tcells"
names(cellTypes)[which(cellTypes == "CD16+/CD11b-")] = "CD16_monocytes"
names(cellTypes)[which(cellTypes == "CD11b+/CD16-")] = "CD14_monocytes" # Classical monocyte
names(cellTypes)[which(cellTypes == "CD11b-/CD16-")] = "CD14neg_CD16neg_cells" # Gabi didn't say what kind of cell this includes
names(cellTypes)[which(cellTypes == "CD16+/CD11b+")] = "CD14_CD16_monocytes"

# Oops, i need the reverse
processedNames = names(cellTypes)
names(processedNames) = cellTypes

# Change Population to be our accepted formatted names
newPop = as.character(rawCounts$Population)
newPop = processedNames[newPop]
rawCounts$oldPopulation = rawCounts$Population
rawCounts$Population = as.factor(newPop)
#............................#
##### Pre-process Counts #####
#............................#
## Extract sampleID
sampleID = strsplitVector(as.character(rawCounts$File), mySplit = "_", whichPart = 1)
rawCounts$sampleID = as.factor(sampleID)

## Extract Donor
donor = strsplitVector(as.character(rawCounts$File), mySplit = "_", whichPart = 2)
rawCounts$donor = as.factor(donor)

## Are all sampleIDs from rawCounts also in pheno? Yes!
all(as.character(rawCounts$sampleID) %in% as.character(pheno$sampleid)) # True!

## Are all donors from rawCounts also in pheno? Yes!
all(as.character(rawCounts$donor) %in% as.character(pheno$donor))# True!

## Is there the same number of samples as donors? 
# aka is there one sample per donor? 
length(unique(rawCounts$donor)) # 84 donors
length(unique(rawCounts$sampleID)) # 86 samples

## Who has duplicates? 
# 7605 and 7962 have 28 samples where all others have 14 samples
# They just seem like normal people, nothing special about them
# I can take the average counts for them
table(rawCounts$donor)[which(table(rawCounts$donor) != 14)]

#................................................#
##### Create new matrix in the style of expr #####
#................................................#
# columns = sample ID
# rows = cell type

cellCounts = acast(rawCounts,Population ~ donor,value.var = "Event_Count", fun.aggregate= mean)


#.......................................#
##### Convert counts to proportions #####
#.......................................#
cellProp = cellCounts

granulocyte = c("Neutrophils", "basophils", "singlets")

for(cell in rownames(cellCounts)) {
  if(cell %in% granulocyte){
    # If it's a granulocyte, do as a percentage of singlets
    cellProp[cell,] = cellCounts[cell,] / cellCounts["singlets",]
  } else {
    # if it's not a granulocyte, do as a percentage of mononuclear cells
    cellProp[cell,] = cellCounts[cell,]/ cellCounts["mononuclearCells",]
  }
}

#.....................................#
##### Create a total monocyte row #####
#.....................................#
# We were given a three monocyte subsets: CD14+, CD16+ and CD14+CD16+
# I need to add them together into one large "monocyte" group 

monocyte = colSums(cellProp[c("CD14_monocytes", "CD16_monocytes", "CD14_CD16_monocytes"),])
cellProp = rbind(cellProp, monocyte)

#...............................#
##### Create data for Manoj #####
#...............................#
phenoMatrix = unique(pheno[,c("donor", "gender", "age")])
phenoMatrix = as.matrix(phenoMatrix)
rownames(phenoMatrix) = as.character(phenoMatrix[,"donor"])


# Create sex as 1=female, 0=male
sex = as.character(phenoMatrix[,"gender"])
sex[which(sex == "")] = NA
sex[which(sex == "M")] = "Male"
sex[which(sex== "F")] = "Female"

phenoMatrix = cbind(phenoMatrix, sex)
table(phenoMatrix[,'gender'], phenoMatrix[,'sex']) # they all match

# Add age and sex into cytofResults data frame
cytofResults = as.data.frame(t(cellProp))
cytofResults$sex = phenoMatrix[rownames(cytofResults), "sex"]
cytofResults$age = phenoMatrix[rownames(cytofResults), "age"]

# Remove samples that have NA's
# One sample has NA for age and sex
cytofResults = na.omit(cytofResults)

#...............................#
##### Quick look at results #####
#...............................#
# Purpose:
#   - See if sex diff in monocyte proportion is lost with age

# Most folk over 40 are between 40-50 years old, so may not work so well
cytofResults$age = as.numeric(as.character(cytofResults$age))
hist(cytofResults$age)


# Create group column 
group = ifelse(cytofResults$age <=40, yes = 'Younger', no = 'Older')
group = paste(group, cytofResults$sex)
cytofResults$group = factor(group, levels = unique(group)[c(2,3,1,4)])


pValYoungOld = function(myCell){
  youngPval = t.test(cytofResults[[myCell]][cytofResults$group == 'Younger Female'],
                     cytofResults[[myCell]][cytofResults$group == 'Younger Male'])$p.value
  youngPval = signif(youngPval, 2)
  
  # Older people
  oldPval = t.test(cytofResults[[myCell]][cytofResults$group == 'Older Female'],
                     cytofResults[[myCell]][cytofResults$group == 'Older Male'])$p.value
  oldPval = signif(oldPval, 2)
  
  myResult = c(youngPval, oldPval)
  names(myResult) = c('young', 'old')
  return(myResult)
}


### Check cell type significance
pValYoungOld('monocyte')# young sig, Old close to but not
pValYoungOld("CD14_monocytes") # young sig, old close, but not
pValYoungOld("CD16_monocytes") # neither significant
pValYoungOld("CD14_CD16_monocytes") # neither significant
pValYoungOld("CD4T_cells") # Young close to, old not at all


### Bulk monocytes are sig in young people
ggplot(cytofResults, aes(x=group, y=monocyte)) +
  geom_violin(fill='grey', trim=F)+
  geom_jitter(width=0.1, size=2, aes(col=sex)) + 
  scale_color_manual(values = c('red', 'darkblue'))

## CD14+ monocytes are sign in young people
ggplot(cytofResults, aes(x=group, y=CD14_monocytes)) +
  geom_violin(fill='grey', trim=F)+
  geom_jitter(width=0.1, size=2, aes(col=sex)) + 
  scale_color_manual(values = c('red', 'darkblue')) +
  theme_cowplot() + theme(legend.position = 'none') +
  xlab('') + ylab('CD14+ Monocyte (% of Mononuclear Cells)')

#..............#
##### Save #####
#..............#

write.table(x = cytofResults, file = "cytof_cellProp.txt", sep = "\t" )
