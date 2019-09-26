# 8/16/2017 Download and preprocess gse69832

# Purpose:
# Aging timecourse of iSEXS expression

#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/2_age/postMenopause/gse69832/")

sexLabels = read.delim("gse69832_ageSex.txt", header = T)
#............................#
##### Download -gse69832 #####
#............................#

gse69832 = getGEOData("GSE69832")
gse69832 = gse69832$originalData$GSE69832

# Tidy up Pheno
gse69832$rawPheno = gse69832$pheno
gse69832$pheno = cleanUpPheno(gse69832$rawPheno, T)

# Only 43 entries, but spread from 14-94
View(gse69832$pheno)
hist(gse69832$pheno$age)
range(gse69832$pheno$age)

# Create subjectID
subjID = gse69832$pheno$title
subjID = strsplitVector(subjID, "control ", 2)
gse69832$pheno$subjectID = subjID

# Add in sex labels
sex = tolower(sexLabels$Gender)
names(sex) = sexLabels$Sample_name
gse69832$pheno$sex = sex[as.character(gse69832$pheno$subjectID)]

# Check keys
# look good!
gse69832$keys[1:20]

# Set class
gse69832$class = createClassVector(gse69832$pheno$sex, "female", gse69832$pheno)

# Check expr
# it needs to be log2 and there's negative valuees
# it's pretty messy with weird low outliers once I'm done
boxplot(gse69832$expr,main="gse69832")

# Bring above zero and take the log2
gse69832$expr = gse69832$expr - min(gse69832$expr) + 1
gse69832$expr = log2(gse69832$expr)

boxplot(gse69832$expr,main="gse69832")


#..........................#
##### Check Sex labels #####
#..........................#
# Sex labels are perfectly fine!
# Nothing I have to do to change things

"XIST" %in% gse69832$keys # FALSE
"RPS4Y1" %in% gse69832$keys #TRUE
"KDM5D" %in% gse69832$keys # TRUE

quickViolin(gse69832, "RPS4Y1", "sex")
quickViolin(gse69832, "KDM5D", "sex")
quickScatter(gse69832, "RPS4Y1", "KDM5D", "sex")

imputedSex = imputeSex(gse69832)
any(imputedSex != gse69832$pheno$sex) # False
all(imputedSex == gse69832$pheno$sex) # True

#..............#
##### Save #####
#..............#
save(gse69832, file = "gse69832_spain_timecourse.RData")
