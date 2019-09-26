# July 20, 2018 - Download Piasecka's gene expression data

# Background: 
#   - Patin et al. published flow on Milieu cohort
#   - Piasecka et al. published Nanostring on 560 immune genes

# Purpose: 
#   - Is Piasecka's data usable? 
#        - What is the overlap with iSEXS? 

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/4_sortedCells/cellProp/patin/")

source('/labs/khatrilab/ebongen/sexDifferences2.0/00_tools/general_GEMfunx.R')

library(ggplot2)
library(cowplot)
library(MetaIntegrator)

load('/labs/khatrilab/ebongen/sexDifferences2.0/isexs/sexMetaObj.RData')

#....................#
##### Grab demographics #####
#....................#

demographics = read.delim('demographics.txt')
dim(demographics)
View(demographics)


#...........................#
##### Download Piasecka #####
#...........................#
# Purpose: 
#  - prep nanostring data

# Load data frame
piasecka = read.delim('piasecka_nanostring.txt')
dim(piasecka) # 562 columns, 5.6k rows
colnames(piasecka) # gene names

# Subset down to just non-stimulated
piasecka = subset(piasecka, stimulus == 'NS')
dim(piasecka) # 562 columns, 804 rows


# get sex dictionary 
sex = demographics$Sex
names(sex) = as.character(demographics$SUBJID)
piasecka$sex = sex[as.character(piasecka$id)]

### Sanity check sex
# Paper highlighted CD14 higher in males - we see it!
ggplot(piasecka, aes(x=sex, y=CD14)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(width=0.1, size=2, aes(col=sex)) + theme_cowplot()
t.test(piasecka$CD14[piasecka$sex == 'Male'],
       piasecka$CD14[piasecka$sex == 'Female']) # way significant

# Is ICOS higher in females? Yup!
ggplot(piasecka, aes(x=sex, y=ICOS)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(width=0.1, size=2, aes(col=sex)) + theme_cowplot()
t.test(piasecka$ICOS[piasecka$sex == 'Male'],
       piasecka$ICOS[piasecka$sex == 'Female']) # way significant

# get age dictionary
age = demographics$Age
names(age) = as.character(demographics$SUBJID)
piasecka$age = age[as.character(piasecka$id)]

### Sanity check age
# Only works for flu, not N.S.
# So, not worth digging into
test = subset(piasecka, age <40)
test$group = ifelse(test$age <30, yes = 'twenties', no = 'thirties')

#  IFNA2 - weird, not significant
ggplot(test, aes(x=group, y=IFNA2)) + geom_boxplot()
t.test(test$IFNA2[test$group == 'twenties'],
       test$IFNA2[test$group == 'thirties']) # crazy significant

# FCGRT
ggplot(test, aes(x=group, y=FCGRT)) + geom_boxplot()
t.test(test$FCGRT[test$group == 'twenties'],
       test$FCGRT[test$group == 'thirties']) # crazy significant

#.........................#
##### Dataset Object #####
#.........................#
# Purpose: 
#   - Turn Piasecka into a dataset object 


# Expr
expr = piasecka[,colnames(piasecka)[!colnames(piasecka) %in% c('age', 'sex', 'stimulus', 'id')]]
rownames(expr) = paste('sample', as.character(piasecka$id), sep = "_")
expr = t(expr)
is.matrix(expr) # true

### Pheno
pheno = piasecka[,c('id', 'age', 'sex', 'stimulus')]
rownames(pheno) = paste('sample', as.character(pheno$id), sep = '_')
all(rownames(pheno) == colnames(expr)) # True
is.data.frame(pheno) # true

### class
myClass = ifelse(pheno$sex == 'Female', yes = 1, no = 0)
table(myClass, pheno$sex)
names(myClass) = rownames(pheno)

### Keys
myKeys = rownames(expr)
names(myKeys) = rownames(myKeys)
tail(myKeys)

### Dataset obj
milieuGEM = list(expr = expr,
                 pheno = pheno, 
                 class = myClass, 
                 keys = myKeys, 
                 formattedName = 'Milieu Interieur')

# Remove samples with NA sex
phenoNAs = subset(milieuGEM$pheno, !is.na(sex))
dim(phenoNAs)
milieuGEM = subsetGEMFromPheno(milieuGEM, phenoNAs)

# Looks good!
checkDataObject(milieuGEM, 'Dataset')
#..............................#
##### Analyze Young people #####
#..............................#
# Purpose: 
#   - Only keep people 18-40 years old
#   - Get effect sizes
#   - Compare effect sizes to sexMetaObj and valMetaObj 
#     effect sizes

# 1 # Remove oldies
phenoYoung = subset(milieuGEM$pheno, age <=40)
hist(phenoYoung$age)
range(phenoYoung$age) # even range 20-40

# Age distributions aren't so much the same
# but, mostly comparable
hist(phenoYoung$age[phenoYoung$sex == 'Male'])
hist(phenoYoung$age[phenoYoung$sex == 'Female'])

milieuGEM_young  = subsetGEMFromPheno(milieuGEM, phenoYoung)



# 2 # Run meta-analysis for effect sizes
milieuMeta = list(originalData = list(milieuGEM_young = milieuGEM_young, extra = milieuGEM))
milieuMeta = runMetaAnalysis(milieuMeta, runLeaveOneOutAnalysis = F)

# 3 # Get effect sizes

# Get isexs genes that are shared
isexs = c(sexMetaObj$filterResults$isexs$posGeneNames, sexMetaObj$filterResults$isexs$negGeneNames)
isexs = isexs[isexs %in% myKeys]
length(isexs) # 13 genes are shared

# Create matrix of effect sizes
isexs_shared = cbind(milieuMeta$metaAnalysis$datasetEffectSizes[isexs,"milieuGEM_young"],
             sexMetaObj$metaAnalysis$pooledResults[isexs,'effectSize'],
             valMetaObj$metaAnalysis$pooledResults[isexs,'effectSize'])

colnames(isexs_shared) = c('milieu', 'discovery', 'validation')
isexs_shared = as.data.frame(isexs_shared)

# 4 # Plot the comparison

a = ggplot(isexs_shared, aes(x=discovery, y=milieu)) + 
  geom_hline(yintercept = 0, col='darkgrey', size=1.5, linetype='dashed') + 
  geom_vline(xintercept = 0, col='darkgrey', size=1.5, linetype='dashed')+
  geom_point(size=3)+
  theme_cowplot() + ggtitle('Milieu Cohort vs Discovery Effect Sizes') +
  xlab('Discovery Effect Sizes') + ylab('Milieu Interieur Effect Sizes')


b = ggplot(isexs_shared, aes(x=validation, y=milieu)) + 
  geom_hline(yintercept = 0, col='darkgrey', size=1.5, linetype='dashed') + 
  geom_vline(xintercept = 0, col='darkgrey', size=1.5, linetype='dashed')+
  geom_point(size=3)+
  theme_cowplot()+ ggtitle('Milieu Cohort vs Validation Effect Sizes') +
  xlab('Validation Effect Sizes') + ylab('Milieu Interieur Effect Sizes')

milieuPlot = plot_grid(a,b, labels = c('a','b'), nrow = 1, ncol = 2)

save_plot(filename = 'milieuPlot.pdf',milieuPlot,ncol = 2, nrow=1)

#..............#
##### Save #####
#..............#
save(piasecka, file='piasecka_nanostring.RData')

save(milieuMeta, file = 'milieuMeta.RData')
