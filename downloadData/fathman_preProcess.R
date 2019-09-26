# 8/16/2017 gse74752 Timecourse from 40 to 100 yo

# Fathman's datasets: 
# GSE74752
# GSE73089

# Conclusions:
# XIST and RPS4Y1 do not spearate samples into two lear clusters
# Sex labels are a total mess
# So, I can't use these datasets


#....................#
##### Load Stuff #####
#....................#

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")

source("00_tools/general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd("0_datasets/2_age/postMenopause/gse74752/")

#............................#
##### Download -gse73089 #####
#............................#

gse73089 = getGEOData("GSE73089")
gse73089 = gse73089$originalData$GSE73089

# Clean up pheno
gse73089$rawPheno = gse73089$pheno
gse73089$pheno = cleanUpPheno(gse73089$rawPheno, T)

View(gse73089$pheno)

# Create class
gse73089$class = createClassVector(gse73089$pheno$sex, "female", gse73089$pheno)

# Check Keys
# Lots of NAs, but looks okay
gse73089$keys[1:20]
length(gse73089$keys)

# Raise Expr above zero
gse73089$expr = gse73089$expr - min(gse73089$expr, na.rm = T) + 1

# Check expr
# Agilent
# Medians around zero, ugly and messy, but no obvious batch effect
boxplot(gse73089$expr[,1:10], main = "first 10 samples gse73089")
pdf("gse73089_pbmc_boxplot.pdf")
boxplot(gse73089$expr, main = "gse73089 PBMC")
dev.off()



# Set formatted name
gse73089$formattedName = "GSE73089_Fathman_PBMC"

#...........................................#
##### Remove Control Samples - GSE73089 #####
#...........................................#
# There's 10 control samples that seem to come from teh same individual
# I think they're technical controls, so I'm removing all of them

phenoNoCont = subset(gse73089$pheno, individual == "patient")
gse73089 = subsetGEMFromPheno(gse73089, phenoNoCont)


#............................#
##### Download -GSE74752 #####
#............................#
# Download the WB version of the dataset
# There are duplicate patients, but I dunno to what extent

gse74752 = getGEOData("GSE74752")
gse74752 = gse74752$originalData$GSE74752

# Clean up pheno
gse74752$rawPheno = gse74752$pheno
gse74752$pheno = cleanUpPheno(gse74752$rawPheno, T)

View(gse74752$pheno)

# Create Class
gse74752$class = createClassVector(gse74752$pheno$sex, "female", gse74752$pheno)

# Check keys
gse74752$keys[1:20]
length(gse74752$keys) # 45016, 

# Formatted name
gse74752$formattedName = "GSE74752_WB"

# Check expr
# Pretty messed up, here's how to fix it:
# 1) Remove last probe (entirely NA)
# 2) Raise so min is 1
# 3) Take the log2
# 4) Check boxplot

# Remove last probe
gse74752$expr = gse74752$expr[1:nrow(gse74752$expr)-1,]
gse74752$keys = gse74752$keys[rownames(gse74752$expr)]

# Raise so min is 1
gse74752$expr = gse74752$expr - min(gse74752$expr, na.rm = T) +1

# take log2
gse74752$expr = log2(gse74752$expr)


pdf("gse74752_WB_boxplot.pdf")
boxplot(gse74752$expr, main = "gse74752 - Whole Blood")
dev.off()

#..........................#
##### Check Sex Labels #####
#..........................#
# XIST and RPS4Y1 do not spearate into two clusters
# I don't think this dataset is usable


# TRUE
"XIST" %in% gse73089$keys
"XIST" %in% gse74752$keys

# TRUE
"RPS4Y1" %in% gse73089$keys
"RPS4Y1" %in% gse74752$keys

# FALSE
"KDM5D" %in% gse73089$keys
"KDM5D" %in% gse74752$keys


quickScatter(gse73089, "RPS4Y1", "XIST", "sex")
quickScatter(gse74752, "RPS4Y1", "XIST", "sex")

#..............#
##### Save #####
#..............#

save(gse73089, gse74752)