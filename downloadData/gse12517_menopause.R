# October 2018 - Preprocess gse21517 Menopause

# Purpose: 
#   - I previously found that there's no diff in auto-iSEXS with menopause
#   - I abandoned that finding at first, but now it's useful!
#   - Negative data means something!

#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/2_age/postMenopause/gse12517_menopause/")


#............................#
##### Download and Pheno #####
#............................#
gse12517 = getGEOData("GSE12517")
gse12517 = gse12517$originalData$GSE12517

# Clean up pheno
gse12517$rawPheno = gse12517$pheno
gse12517$pheno = cleanUpPheno(gse12517$rawPheno, T)

# Make nicer treatment column
treatment = as.character(gse12517$pheno$title)
treatment[which(grepl(x = treatment, pattern = "Aromatase", ignore.case = T))] = "Aromatase Inhibitor"
treatment[which(grepl(x = treatment, pattern = "Tamoxifen", ignore.case = T))] = "Tamoxifen"
treatment[which(!treatment %in% c("Aromatase Inhibitor", "Tamoxifen"))] = "None"
gse12517$pheno$treatment = treatment

# Create group, where first group is pre-menopause
group = paste(gse12517$pheno$Stage, gse12517$pheno$treatment)
group = factor(group,levels = unique(group))
gse12517$pheno$group = group

# Create group2, where first group is post-menopause
group = paste(gse12517$pheno$Stage, gse12517$pheno$treatment)
group = factor(group,levels = unique(group)[c(2,3,4,1)])
gse12517$pheno$group2 = group

#..............#
##### Expr #####
#..............#
# Results:
#   - No NAs
#   - 2-14 range
#   - No major batch effect


# Looks great!
gse12517$pheno$platform_id[1] # gpl96, yay, Affy!

range(gse12517$expr) # 2-14

# Is there batch effect?
boxplot(gse12517$expr, main = "gse12517")

# Any males? - Nope, just a giant cloud!
quickScatter(gse12517, 'RPS4Y1', 'XIST', 'treatment')

#..............#
##### iSEXS ####
#..............#

# Full iSEXS - maybe some patterns?
violinPlot(sexMetaObj$filterResults$isexs, gse12517, "group")

# XY-iSEXS - no Differences
violinPlot(sexMetaObj$filterResults$xySig, gse12517, "group")

# Auto-iSEXS - no difference between pre- and post-menopause
violinPlot(sexMetaObj$filterResults$autoSig, gse12517, "group")


# Auto-iSEXS - Significantly lower scores in aromatase inhibitor group
violinPlot(sexMetaObj$filterResults$autoSig, gse12517, "group2")

#..............#
##### Save #####
#..............#

save(gse12517, file='gse12517.RData')



