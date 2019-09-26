# May 3 2018 - Obtain chromosome labels for all genes

# Background: 
#   - The package I use to get chromosome locations fails randomly
#      and takes a long time to run
#   - It can throw a wrench if it's part of a larger script, cuz
#     then the larger script keeps failing

#....................#
##### Load Stuff #####
#....................#

library(biomaRt) # For chromosomal location of genes

setwd("/labs/khatrilab/ebongen/sexDifferences2.0/isexs/reference/")


load('sexMetaObj.RData') # Load iSEXS
#..................................#
##### Get chromosomal location #####
#..................................#

# Keep hitting this line until it goes through
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Get location of all genes from biomart
allGeneLocations <- getBM(attributes=c('hgnc_symbol','chromosome_name'), filters = 'chromosome_name', values =c(as.character(1:22), "X", "Y"), mart = ensembl)

# Create a dictionary to easily get ISEXS locations
allGeneLocations_dict = allGeneLocations$chromosome_name
names(allGeneLocations_dict) = allGeneLocations$hgnc_symbol

save(allGeneLocations_dict, file= "allGeneLocations_dict.RData")


#...........................................#
##### Are X and Y Chromosomes enriched? #####
#...........................................#
# Purpose: 
#    - Are more X and Y genes in iSEXS than due to random chance?
# Results: 
#    - Yup! Both are crazy significant


### Test Y-chromosome genes
yData = rbind(c(length(ISEXSLoc), sum(ISEXSLoc == 'Y')), c(nrow(allGeneLocations), sum(allGeneLocations$chromosome_name == 'Y')))
colnames(yData) = c('all', 'Y')
rownames(yData) = c('iSEXS', 'genome')
chisq.test(yData) # p < 3.39E-7


### Test X-chromosome genes
xData = rbind(c(length(ISEXSLoc), sum(grepl(pattern = 'X', x = ISEXSLoc))), 
              c(nrow(allGeneLocations), sum(allGeneLocations$chromosome_name == 'X')))
colnames(xData) = c('all', 'X')
rownames(xData) = c('iSEXS', 'genome')
chisq.test(xData) # p < 3.66E-11
