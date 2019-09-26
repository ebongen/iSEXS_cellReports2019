# Explore GSE56033 - Naieve CD4+ T cells

# Results: 
# CCR7 may be higher expressed in female T cells
# Other genes look disappointing
# Autosomal Signature diverges 18-40, lost after 50
# XY Signature is friggin weird
#   - 3 tiers of female expression
#   - Most likely an artifact
#   - BUt, there's an off chance that it has to do with cycle

#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

library(GEOquery)
library(MetaIntegrator)
source("00_tools/general_GEMfunx.R")
source("00_tools/plot_funx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/4_sortedCells/cd4_Tcells/gse56033/")

#..........................#
##### Download Dataset #####
#..........................#
gse56033 = getGEOData("GSE56033")
gse56033 = gse56033$originalData$GSE56033

gse56033$rawPheno = gse56033$pheno
gse56033$pheno = cleanUpPheno(gse56033$rawPheno, T)

# Remove redundant columns
redundantCols = colnames(gse56033$pheno)[grepl(pattern = ":ch1", x = colnames(gse56033$pheno))]
gse56033$pheno = gse56033$pheno[,!colnames(gse56033$pheno) %in% redundantCols]

# Create subject ID 
subjectID = strsplitVector(gse56033$pheno$title, "\\.", 1)
length(unique(subjectID)) #454, 479 should be unique samples
dim(unique(cbind(subjectID, gse56033$pheno$age))) # 461, doesn't get us all the way up to 479
gse56033$pheno$subjectID = subjectID

# Clean up age
gse56033$pheno$age = as.numeric(gse56033$pheno$`age_(yrs)`)
gse56033$pheno$`age_(yrs)` = NULL
hist(gse56033$pheno$age) # Mostly 18-60
min(gse56033$pheno$age)

#..............#
##### Expr #####
#..............#
# Purpose:
#   Look for batch effect
#   Examine amount of NAs
#   Move al values above zero

# No obvious batch effect
boxplot(gse56033$expr, main = "gse56033")

# Normal amount of NAs
sum(is.na(gse56033$expr)) # 868 NAs

# Move expr above zero
min(gse56033$expr, na.rm = T) # min is -6.27
gse56033$expr = gse56033$expr - min(gse56033$expr, na.rm = T) + 1

#....................#
##### Sex Labels #####
#....................#

# Check sex labels
"XIST" %in% gse56033$keys # false
"RPS4Y1" %in% gse56033$keys # true
"KDM5D" %in% gse56033$keys # true


# Grab expression for Y-chr genes
rps4y1 = getSampleLevelGeneData(gse56033, "RPS4Y1")
kdm5d = getSampleLevelGeneData(gse56033, "KDM5D")

# RPS4Y1 has NAs, so hard to impute sex with it
sum(is.na(rps4y1)) # 50 samples have NA for RSP4Y1
sum(is.na(kdm5d))

# Clear separation via KDM5D alone
quickViolin(gse56033, "KDM5D", "sex") # clear separation by KDM5D alone
imputedSex = imputeSex(gse56033, malGenes = "KDM5D")
any(imputedSex != gse56033$pheno$sex) # nope! No mislabled samples!

# IGTB187 says its age is 9 yo
# But, the paper says everyone was 18-50yo
# I looked in gse56034 (monocyte cohort), and IGTB187 is a 49yo Male
# gse56033 says he's a 9yo male, so I'm guessing they lost the 4
# So, I'll put his age as 49
kids = subjectID[gse56033$pheno$age < 18]
gse56033$pheno$age[gse56033$pheno$subjectID == kids] = 49


#......................#
##### Class Vector #####
#......................#

gse56033$class = createClassVector(gse56033$pheno$sex, "female", gse56033$pheno)


#..............................#
##### Full ISEXS signature #####
#..............................#

pdf("gse56033_sig_timecourses.pdf", width = 10)
timecourse_iSEXS(gse56033, sexMetaObj$filterResults$xy, "age", "sex", "years")
timecourse_iSEXS(gse56033, sexMetaObj$filterResults$autosomeOnly, "age", "sex", "years")
dev.off()

#......................................#
##### Exploring 4 XY-gene clusters #####
#......................................#
# Purpose: 
# Classify each sample as being either:
# Female hi, female mid, female low,or male



xyScore = calculateScore(datasetObject = gse56033, filterObject = sexMetaObj$filterResults$xy)

set.seed(0)
xyGroups = kmeans(xyScore, 4)
xyGroups$centers
groupDict = c("male", "fem1lo", "fem2mid", "fem3hi")
names(groupDict) = as.character(order(xyGroups$centers))

xyGroupLabels = groupDict[as.character(xyGroups$cluster)]
gse56033$pheno$xyGroups = xyGroupLabels

# XY-score is different, but autosomal score not so much
pdf("gse56033_xyScore_violins.pdf", width = 8)
violinPlot(sexMetaObj$filterResults$xy, gse56033, "xyGroups")+ ggtitle("GSE65033 - XY ISEXS")
violinPlot(sexMetaObj$filterResults$autosomeOnly, gse56033, "xyGroups") + ggtitle("GSE65033 - Autosomal ISEXS")
dev.off()

### Do the different groups separate by PCA? 
# Could this be some sort of batch effect
# Answer: I don't think so, no clear separation
pcaData = gse56033$expr
pcaData[is.na(pcaData)] = 0
pcaTotal = prcomp(t(pcaData), scale. = T, center = T)

# No obvious clustering of samples by XY expression if you look at all samples
pdf("gse56033_xyScore_pca.pdf", width = 10)
ggplot() + geom_point(aes(x=pcaTotal$x[,1], y=pcaTotal$x[,2], col = xyGroupLabels)) +
  ggtitle("GSE56033 - PCA of All Genes") + ylab("PC2") + xlab("PC1")
dev.off()

#......................#
##### Divvy by Age #####
#......................#
### Young cohort, 18-40 yo
phenoYoung = subset(gse56033$pheno, age <=40)
range(phenoYoung$age) # 18-40
table(phenoYoung$sex) # 255 females, 143 males
t.test(phenoYoung$age[phenoYoung$sex == "female"], phenoYoung$age[phenoYoung$sex == "male"]) # p = 0.13
gse56033_young = subsetGEMFromPheno(gse56033, phenoYoung)


### 18-25yo
# Significant difference in age, but it's such a tight range, I don't think it'll ruin anything
pheno_25 = subset(gse56033$pheno, age <= 25)
range(pheno_25$age) # 18-25
table(pheno_25$sex) # 145 female, 74 male
t.test(pheno_25$age[pheno_25$sex == "female"], pheno_25$age[pheno_25$sex == "male"]) # p = 0.009
ggplot(data = pheno_25, aes(age, col = sex)) + geom_density(size = 2) + ggtitle("pheno_25") + scale_color_manual(values = c(femColor, malColor))
gse56033_25yo = subsetGEMFromPheno(gse56033, pheno_25)
gse56033_25yo$formattedName = "GSE56033 18-25yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56033_25yo) # AUC = 0.52


### 26-30yo
pheno_30 = subset(gse56033$pheno, age <= 30 & age > 25)
range(pheno_30$age) #26-30
table(pheno_30$sex) # 72 females, 32 males
t.test(pheno_30$age[pheno_30$sex == "female"], pheno_30$age[pheno_30$sex == "male"]) # p = 0.8
gse56033_30yo = subsetGEMFromPheno(gse56033, pheno_30)
gse56033_30yo$formattedName = "GSE56033 26-30yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56033_30yo) # AUC = 0.63

### 30-40yo
pheno_40 = subset(gse56033$pheno, age <= 40 & age > 30)
range(pheno_40$age) # 31-40
table(pheno_40$sex) # 38 females, 37 males
t.test(pheno_40$age[pheno_40$sex == "female"], pheno_40$age[pheno_40$sex == "male"]) # p = 0.03
ggplot(data = pheno_40, aes(age, col = sex)) + geom_density(size = 2) + ggtitle("pheno_40") + scale_color_manual(values = c(femColor, malColor))
ggplot(data = pheno_40, aes(age, fill = sex)) + geom_histogram() + ggtitle("pheno_40") + scale_fill_manual(values = c(femColor, malColor))
gse56033_40yo = subsetGEMFromPheno(gse56033, pheno_40)
gse56033_40yo$formattedName = "GSE56033 31-40yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56033_40yo) # 0.76

### 30-35yo
pheno_35 = subset(gse56033$pheno, age <= 35 & age > 30)
range(pheno_35$age)
table(pheno_35$sex) # 26 females 18 males
t.test(pheno_35$age[pheno_35$sex == "female"], pheno_35$age[pheno_35$sex == "male"]) # p = 0.02
gse56033_35yo = subsetGEMFromPheno(gse56033, pheno_35)
gse56033_35yo$formattedName = "GSE56033 30-35yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56033_35yo) # 0.81

### 36-40yo
pheno_3640 =  subset(gse56033$pheno, age <= 40 & age > 35)
range(pheno_3640$age) # 36-40
table(pheno_3640$sex) # 12 females, 19 males

gse56033_3640yo = subsetGEMFromPheno(gse56033, pheno_3640)
gse56033_3640yo$formattedName = "GSE56033 36-40yo"
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse56033_3640yo) # 67

### 41-50
pheno_50 = subset(gse56033$pheno, age <=50 & age > 40)
table(pheno_50$sex) # 31 females, 54 males
t.test(pheno_50$age[pheno_50$sex == "female"], pheno_50$age[pheno_50$sex == "male"]) # p = 0.5
gse56033_50yo = subsetGEMFromPheno(gse56033, pheno_50)
gse56033_50yo$formattedName = "GSE56033 41-50yo"

### Over  50 
pheno_54 = subset(gse56033$pheno, age > 50)
table(pheno_54$sex) # 8 males, 8 females
t.test(pheno_54$age[pheno_54$sex == "female"], pheno_54$age[pheno_54$sex == "male"]) # p = 0.17
gse56033_54yo = subsetGEMFromPheno(gse56033, pheno_54)
gse56033_54yo$formattedName = "GSE56033 51-54yo"


#..............#
##### Save #####
#..............#
save(gse56033, gse56033_young,file = "gse56033_full.RData")
save(gse56033_25yo, gse56033_30yo, 
     gse56033_35yo, gse56033_3640yo, gse56033_40yo, 
     gse56033_50yo, gse56033_54yo, file = "gse56033_ageBins.RData")
