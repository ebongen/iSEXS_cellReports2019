# December 2017 - Download and preprocess gse70044

# Background: 
# I want to see what isexs is like in sorted cells
# gse70044 has neutrophils from college aged Singaporeans 

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("~/Tools/Graphing Scripts/quickPlots.R")
source("00_tools/general_GEMfunx.R")
load("1_metaAnaly/sexMetaObj.RData")

setwd("0_datasets/4_sortedCells/neutrophils/gse70044/")

# I had to download gse70044 from local machine
# So, load the Meta-object from there
load("gse70044_raw.RData")

#...............#
##### Pheno #####
#...............#

# Clean up pheno
gse70044$rawPheno = gse70044$pheno
gse70044$pheno =cleanUpPheno(gse70044$rawPheno, T)


# Some of the entries are misssing age and sex annotation
# So, their columns got a little shifted
gse70044$pheno$raceEthnicity = rep("chinese", nrow(gse70044$pheno))
gse70044$pheno$disease_status = rep("healthy", nrow(gse70044$pheno))
gse70044$pheno$ageUnits = rep("yr", nrow(gse70044$pheno))
gse70044$pheno$sex[!gse70044$pheno$sex %in% c("male", "female")] = NA
gse70044$pheno$sex = as.factor(as.character(gse70044$pheno$sex))

# Examine age
# One person's age is listed as 0
# The youngest person in the paper is 19 years old
# So, we'll set the 0 aged person to NA
sort(gse70044$pheno$age)
gse70044$pheno$age[gse70044$pheno$age == 0] = NA
sum(is.na(gse70044$pheno$age)) # Now 10 people have NAs


#..............#
##### Expr #####
#..............#
# Ranges 6-15
# No obvious batch effect
boxplot(gse70044$expr, main = "gse70044")
range(gse70044$expr, na.rm = T)
sum(is.na(gse70044$expr)) # Only 1 NA


#....................#
##### Impute sex #####
#....................#
# A few samples are missing sex annotation
# Since this cohort is 18-29 years old, we can impute sex on the unknown samples
# Since, we know what range they're in, even if we don't know their actual ages

# Look for sex labeling probes
"XIST" %in% gse70044$keys #t
"RPS4Y1" %in% gse70044$keys #T
"KDM5D" %in% gse70044$keys# F

# Impute sex
imputedSex = imputeSex(gse70044)

# All known labels are correct
table(imputedSex, as.character(gse70044$pheno$sex)) # They match perfect! No need to relabel

# Add labels to samples given as NA
gse70044$pheno$sex[is.na(gse70044$pheno$sex)] = imputedSex[is.na(gse70044$pheno$sex)]
table(gse70044$pheno$sex) # 46 females and 68 males

# They separate well by Y-chromosome genes!
quickViolin(gse70044, "XIST", "sex") # Doesn't separate well
quickViolin(gse70044, "RPS4Y1", "sex") # Separates beautifully!

# Create class vector
gse70044$class = createClassVector(gse70044$pheno$sex, "female", gse70044$pheno)


#............................#
##### Difference in ages #####
#............................#
# Problem: statistically significant difference in distribution 
#          male and female ages
# Observations:
#   Females 18-22
#   Males: 18-24/29
#
# Solution:
#   Sure, they have different distributions, but these ages
#   aren't all that different biologically, so probably no 
#   big skewing effect

# There is s t.test difference between ages of males and females
t.test(gse70044$pheno$age[gse70044$pheno$sex == "female"], gse70044$pheno$age[gse70044$pheno$sex == "male"]) # p = 0.008

# Yeah, their distributions are different
hist(gse70044$pheno$age[gse70044$pheno$sex == "female"], main = "female ages")
hist(gse70044$pheno$age[gse70044$pheno$sex == "male"], main = "male ages")

# Even if you remove the outlier, it's still significant
femAge = na.omit(gse70044$pheno$age[gse70044$pheno$sex == "female"])
malAge = na.omit(gse70044$pheno$age[gse70044$pheno$sex == "male"])
malAge = malAge[malAge != max(malAge)]

hist(femAge)
hist(malAge)
t.test(malAge, femAge) # Still significant

#...........................#
##### Quick ISEXS Check #####
#...........................#
# Works surprsingly well in autosomal genes!
rocPlot(sexMetaObj$filterResults$xy, gse70044) # AUC = 1.0
rocPlot(sexMetaObj$filterResults$autosomeOnly, gse70044) # AUC = 0.72

#..............#
##### Save #####
#..............#
save(gse70044, file = "gse70044_processed.RData")
