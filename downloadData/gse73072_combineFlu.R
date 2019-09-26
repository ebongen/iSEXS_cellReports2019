# July 3rd 2018 - Combine influenza challenge studies from gse73072

# Background: 
#   - By combining challenge studies and only looking at
#     symptomatic-shedders, we can examine sex differences in how auto-iSEXS
#     changes during infection
# Purpose:
#   - Create a mega-dataset object of Challenges A, B, and C from gse73072

#....................#
##### Load Stuff #####
#....................#

# Load libraries
library(sva)
library(MetaIntegrator)
library(cowplot)

# source my code
source('/labs/khatrilab/ebongen/viralChallenge_clean/00_tools/general_GEMfunx.R')

# Load quick graphing code
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/6_infection/infection/gse73072/')

#............................#
##### Prep fullChallenge #####
#............................#
# Purpose: 
#   - copy code from viral challenge code v6.0
#   - Create a list of timecourse dataset objects
#   - Only include DEE2, DEE3, and DEE5
#   - Only use influenza with roughly equal numbers of symptomatic-shedder male and females
#   - Exclude DEE4 because all symptomatic-shedders are female

# Load full challenges
load('/labs/khatrilab/ebongen/viralChallenge_clean/0_data/gse73072_2018/gse73072_fullCohorts.RData')

# Remove unnecessary challenges
# Datasets excluded if not flu
# DEE4 excluded because only females are symptomatic shedders
# I can't tell apart batch effect from sex effect
rm(gse73072_dee1, gse73072_dee4, gse73072_duke, gse73072_uva)

# fix formatted names
gse73072_dee2$formattedName = "Challenge B - H3N2"
gse73072_dee3$formattedName = "Challenge A - H1N1"
gse73072_dee5$formattedName = "Challenge C - H3N2"


# Create list of datasets where we want to compare 
# male and female sick people
fullChallenge = list(gse73072_dee3 = gse73072_dee3,
                     gse73072_dee2 = gse73072_dee2,
                     gse73072_dee5 = gse73072_dee5)

# Remove weird samples
for(i in 1:length(fullChallenge)){
  phenoClean = subset(fullChallenge[[i]]$pheno, group %in% c("symptShed", "asympNonShed"))
  fullChallenge[[i]] = subsetGEMFromPheno(fullChallenge[[i]], phenoClean)
  print(fullChallenge[[i]]$formattedName)
  print(table(fullChallenge[[i]]$pheno$group))
}

#........................#
##### Prep comboTime #####
#........................#
# Purpose: 
#   - copy code from viral challenge project
#   - create a combined timecourse Dataset object
#     from Challenges A, B, and C from gse73072

# Are the genes the same?
genes2 = unique(fullChallenge$gse73072_dee2$keys)
genes3 = unique(fullChallenge$gse73072_dee3$keys)
genes5 = unique(fullChallenge$gse73072_dee5$keys)
all(genes2 == genes3) # true
all(genes2 == genes5) # true

# Get expression of all genes
expr2 = getSampleLevelGeneData(fullChallenge$gse73072_dee2, genes2)
expr3 = getSampleLevelGeneData(fullChallenge$gse73072_dee3, genes2)
expr5 = getSampleLevelGeneData(fullChallenge$gse73072_dee5, genes2)

all(rownames(expr2) == rownames(expr3)) # True
all(rownames(expr2) == rownames(expr5)) # True

# Combine into one expr
expr = cbind(expr2, expr3, expr5)
dim(expr) # 933 samples!

# Vector separating samples by batch
batch = c(rep(1, ncol(expr2)), rep(2, ncol(expr3)), rep(3, ncol(expr5)))

# Run ComBat
expr_combat = ComBat(expr, batch)


# Make pretty plot!
all(colnames(fullChallenge$gse73072_dee2$pheno) == colnames(fullChallenge$gse73072_dee3$pheno)) # True
all(colnames(fullChallenge$gse73072_dee2$pheno) == colnames(fullChallenge$gse73072_dee5$pheno)) # False!!!
myCol = c("subject", "virus", "time_point", "sex", "group")

discData = rbind(fullChallenge$gse73072_dee2$pheno[,myCol],
                 fullChallenge$gse73072_dee3$pheno[,myCol],
                 fullChallenge$gse73072_dee5$pheno[,myCol])

# discData gonna be our pheno
colnames(discData)

# Create Dataset object
myClass = rep(0, nrow(discData))
names(myClass) = rownames(discData)

myKeys = rownames(expr_combat)
names(myKeys) = rownames(expr_combat)

comboTime = list(pheno = discData,
                 expr = as.matrix(expr_combat), 
                 class = myClass,
                 keys = myKeys,
                 formattedName = "Combined Flu Challenges")

checkDataObject(comboTime, 'Dataset') # True!


### Create an object with only symptomatic shedders
table(comboTime$pheno$group)

phenoSympt= subset(comboTime$pheno, group == 'symptShed')
comboTime_sympt  = subsetGEMFromPheno(comboTime, phenoSympt)
checkDataObject(comboTime_sympt, 'Dataset') # works!

# Sanity check 
quickTimecourse(comboTime_sympt, myGene = 'IFI44L', timeColumn = 'time_point', groupColumn = 'sex') + 
  xlim(c(-24, 200))

#.........................................#
##### Why am I looking at SLE filters #####
#.........................................#
# Results: 
#   - All subsets of SLE-signatures increase in females faster
#     and higher than males

load('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/5_autoImmunity/SLE/winn/filterObjects.RData')

myData = comboTime_sympt$pheno
myData$sleFilt = calculateScore(sleFilter, comboTime_sympt)
myData$sleIFN = calculateScore(sleIfnFilter, comboTime_sympt)
myData$sleMystery = calculateScore(sleMysteryFilter, comboTime_sympt)
myData$sleNeut = calculateScore(sleNeutroFilter, comboTime_sympt)

myData = subset(myData, time_point < 200)

# IFN sig increase faster and higher in females
ggplot(myData, aes(x=time_point, y=sleIFN, col = sex))+ geom_smooth()+
  geom_point(size=2.5) + scale_color_manual(values = c(femColor, malColor))+
  theme_cowplot() + xlab('Time (hours)') + ggtitle("GSE73072 - Symptomatic-Shedders") +
  ylab('IFN Genes')

# Full Filter: Rises faster and higher in females
ggplot(myData, aes(x=time_point, y=sleFilt, col = sex))+ geom_smooth()+
  geom_point(size=2.5) + scale_color_manual(values = c(femColor, malColor))+
  theme_cowplot()+ xlab('Time (hours)') + ggtitle("GSE73072 - Symptomatic-Shedders") +
  ylab('SLE Meta-Signature')

# Faster and higher in females
ggplot(myData, aes(x=time_point, y=sleMystery, col = sex))+ geom_smooth()+
  geom_point(size=2.5) + scale_color_manual(values = c(femColor, malColor))+
  theme_cowplot() + xlab('Time (hours)') + ggtitle("GSE73072 - Symptomatic-Shedders")+
  ylab("SLE Mystery Genes")

# Faster and higher in females
ggplot(myData, aes(x=time_point, y=sleNeut, col = sex))+ geom_smooth()+
  geom_point(size=2.5) + scale_color_manual(values = c(femColor, malColor))+
  theme_cowplot() + xlab('Time (hours)') + ggtitle("GSE73072 - Influenza Challenge")

#..............#
##### Save #####
#..............#

gse73072_flu_symptShed = comboTime_sympt
gse73072_flu = comboTime

save(gse73072_flu_symptShed, file = 'gse73072_flu_ComBat.RData')
