# August 2017: Downloading, getting key labels, and ComBat normalizing gse11235
# Based off of Take 4.1 from original Sex Differences folder

### Unique to Take 4: 
# Take 2 focused on mapping and did initial COCONUT
# Take 3 focused on: 
#    - More thouroughly validate probes
#    - Finalize COCONUT by doing mean imputation (like Winn suggested)
#    - Create separate datasets for 
#       - XY-female vs XY-male 
#       - XY-female vs XX-female
# Take 4.1 will: 
#   - Take a more conservative approach to normalization
#       - Remove probes with >50% NAs
#       - Remove probes with too high or too low variance
#       - ComBat with zero imputation
#       - Remove sample P-032, because it's cord blood and has normal T levels


### Notes: 
# Most of the pre-processing is copied from Take 3


### Background:
# gse11235 is a very valueble dataset of blood of androgen insensitivity syndrome humans
# Unfortunately, it's only annotated by "Composite Sequence Identifier"
#     - each probe is identified by a different identifier (e.g. symbol, ID number, description)
# gse11235 came from the Stanford HEEBO arrays, so I'm using their annotation table to map
# gse11235 is split into platforms (I'm guessing different batches of the same arrays), but probes are 
#    just listed in numerical order. Probe IDs are not given. So, I'll have to combine at gene level


### Purpose: 
# 1) Map the Composite Sequence Identifiers to Gene Symbols 
# 2) Combine 4 platforms into one big expr matrix 


### Outline: 
# 1)  Download gse11235 and its annotations from GEO
# 2)  Convert annotation tables to $keys strux Dataset obj
# 3)  Create a giant table of all possible identifiers from "Composite Sequence Identifier" 
# 4)  Create a function that will map Comp Seq Identifier to gene symobls using table from #3
# 5)  Actually map Comp Seq Identifiers to gene sympbols
# 6)  Reformat mapped gene symbols to better match $keys format
# 7)  Replace Comp Seq Ident in $keys with gene symbols 
# 8)  Combine sub-dataset expr into one giant expr
# 9)  Create combined gse11235 Dataset object
# 10) Identify batch effect
# 11) Pre-process for ComBat
# 12) ComBat
# 13) Create Dataset objects
#      - expr_noNA and expr_NA
#      - Separate into useful groups
#           -XY_female vs XY_male
#           -XY_female vs XX_female


#....................#
##### Load Stuff #####
#....................#
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/3_humanPerturbations/gse11235_AIS/")

# Packages for downloading Datasets and anotations
library(MetaIntegrator)
library(GEOmetadb)
library(GEOquery)
library(sva)
require(org.Hs.eg.db)



# Script for graphing gene expression
source("~/Tools/Graphing Scripts/quickPlots.R")
library(ggfortify) # package for PCA plots

# My general functions
source("/labs/khatrilab/ebongen/sexDifferences2.0/00_tools/general_GEMfunx.R")

# Load HEEBO table for gene mapping
heebo= read.csv("HEEBO_Human_Set_v1.00.csv", stringsAsFactors = F)

# Load iSEXS
load("/labs/khatrilab/ebongen/sexDifferences2.0/1_metaAnaly/sexMetaObj.RData")
iSEXS = c(sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$posGeneNames,
          sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0$negGeneNames)


#................................................................#
##### 1) Directly download gse11235 and annotations from GEO #####
#................................................................#
# Download the annotation table from GEO
gse11235_getGEO = getGEO("GSE11235")


# Download gse11235 from GEO
gse11235 = getGEOData("GSE11235")
names(gse11235$originalData) # four different platforms
gse11235 = gse11235$originalData
save(gse11235, file= "gse11235_raw.RData")


# Extract GPLs
gpl6765 = gse11235_getGEO$`GSE11235-GPL6765_series_matrix.txt.gz`@featureData@data
gpl6766 = gse11235_getGEO$`GSE11235-GPL6766_series_matrix.txt.gz`@featureData@data
gpl6767 = gse11235_getGEO$`GSE11235-GPL6767_series_matrix.txt.gz`@featureData@data
gpl6768 = gse11235_getGEO$`GSE11235-GPL6768_series_matrix.txt.gz`@featureData@data

#..................................#
##### 2) Convert gpls to $keys #####
#..................................#
# Purpose: create a function to add key vector to each Dataset object
# Create keys vector using Composite Sequence Identifiers


# addKeys creates a $keys vector for gse11235 using an annotaiton table
# it is specific to gse11235. 
# It preprocesses the Composite Sequence Identifiers by converting empty
# ones to NA, Removing weird \\" characters, and making it all uppercase
#
# Inputs: 
#    gem: a Dataset object for gse11235, obtained using MetaIntegrator
#    annotationTable: the gpl table obtained using getGEO
#
# Outputs: 
#   gem: same as input, except $keys have been annotated using Composite Sequence Identifier
addKeys <- function(gem, annotationTable){
  # Create keys
  keys = as.character(annotationTable$`CompositeSequence Identifier`)
  names(keys) = as.character(annotationTable$ID)
  
  # Replace empty strings "" with NA
  keys = replace(x = keys, list = which(keys == ""), values = rep(NA, length(which(keys==""))))
  
  # remove weird \" characters from keys
  keys = gsub(pattern = "\\\"", x =keys, replacement = "", ignore.case = T)
  
  # Make it all upper case
  keys = toupper(keys)
  
  # Put Keys into your dataset object, and check to make sure everything's okay
  gem$keys = keys
  print(checkDataObject(gem, "Dataset"))
  
  return(gem)
}

# Create individual Dataset Objects
gse11235_gpl6765 = addKeys(gse11235$GSE11235_GPL6765, gpl6765)
gse11235_gpl6766 = addKeys(gse11235$GSE11235_GPL6766, gpl6766)
gse11235_gpl6767 = addKeys(gse11235$GSE11235_GPL6767, gpl6767)
gse11235_gpl6768 = addKeys(gse11235$GSE11235_GPL6768, gpl6768)


# Check if the numerical "probe ID" always maps to the same "comp ID"
allProbeNubs = unique(c(names(gse11235_gpl6765$keys),
                        names(gse11235_gpl6766$keys),
                        names(gse11235_gpl6767$keys),
                        names(gse11235_gpl6768$keys)))

# Get a vector of all probes where its compID does not match across platforms
mismatch = c()
for(myProbe in allProbeNubs){
  mv5 = gse11235_gpl6765$keys[myProbe]
  mv6 = gse11235_gpl6766$keys[myProbe]
  mv7 = gse11235_gpl6767$keys[myProbe]
  mv8 = gse11235_gpl6768$keys[myProbe]
  if(length(unique(c(mv5, mv6, mv7, mv8))) > 1){
    mismatch = c(mismatch, myProbe)
  }
}
length(mismatch) # 42727, 99% of probes

# printProbe prints out the compID for a single probe in each gpl
printProbe <- function(probeID){
  print(paste("gpl_5", "     ", gse11235_gpl6765$keys[probeID]))
  print(paste("gpl_6", "     ", gse11235_gpl6766$keys[probeID]))
  print(paste("gpl_7", "     ", gse11235_gpl6767$keys[probeID]))
  print(paste("gpl_8", "     ", gse11235_gpl6768$keys[probeID]))
}

# It seems like the probeIDs are shared across the gpls
# But, diff gpls could have a compID for the same gene
printProbe(mismatch[50]) # each gpl refers to ATP8B3 by a diff name

# Combine all the possible keys into one compID vector
# I will figure out what gene symbols these compIDs correspond to
compID = c(gse11235_gpl6765$keys, gse11235_gpl6766$keys, gse11235_gpl6767$keys, gse11235_gpl6768$keys)
names(compID) = NULL
compID = na.omit(unique(compID))

# Clean up compID by removing "**" from it
hasStar = sapply(compID, function(x) grepl(pattern = "\\*", x = x))
compID[which(hasStar)] = gsub(pattern = "\\**", replacement = "", x = compID[which(hasStar)])

#..................................................................#
##### 3) Create a table to map composite gene identifiers from #####
#..................................................................#
# Purpose: Tweak the heebo table to optimize for mapping Composite Sequence Identifiers to Gene Symbols
# Do this by: 
#   1) Remove uninformative columns
#   2) Put everything in uppercase
#   3) Remove rows with no Gene Symbol
#   4) Add a UniGene annotation column


# Prune heebo down to the unique entries in columns that might map to compID
# Sequence_ID: all sequence IDs contain "hSQ", none of the compIDs do
usefulColumns = c("LocusLink_ID", "Accession",
                  "Updated_GeneID", "GeneID.v12",
                  "Symbol.v12", "Synonyms.v12", "description.name.v12",
                  "Full.name.from.nomenclature.authority.v12")
heebo = heebo[,usefulColumns]
heebo = unique(heebo)

# Make description all caps, to match compID
heebo$description.name.v12 = toupper(heebo$description.name.v12)
heebo$Full.name.from.nomenclature.authority.v12 = toupper(heebo$Full.name.from.nomenclature.authority.v12)

# Remove rows with gene symbols == "."
heebo = heebo[which(heebo$Symbol.v12 != "."),]
length(unique(heebo$Symbol.v12)) # 24281 

# Replace - with NA in Synonym and Full Name columns
heebo$Synonyms.v12 = replace(heebo$Synonyms.v12, list = which(heebo$Synonyms.v12== "-"), values = rep(NA, length(which(heebo$Synonyms.v12== "-"))))
heebo$Full.name.from.nomenclature.authority.v12 = replace(heebo$Full.name.from.nomenclature.authority.v12, 
                                                          list = which(heebo$Full.name.from.nomenclature.authority.v12== "-"), 
                                                          values = rep(NA, length(which(heebo$Full.name.from.nomenclature.authority.v12== "-"))))


# Put Unigene annotatations into HEEBO
# I think this gets ~7,000 of ~10,000 Unigene codes in Composite identifier
unigeneAnnot <- toTable(org.Hs.egUNIGENE)
unigene <- unigeneAnnot$unigene_id
names(unigene) <- unigeneAnnot$gene_id

heebo$unigene = toupper(unigene[heebo$GeneID.v12])

#........................................#
##### 4) Create Function for Mapping #####
#........................................#
# mapOneGene takes a composite sequence identifier and find what gene symbol
# it maps to using heebo. 
#
# To save time, it uses clever rules: 
# 1) If contains "HS." --> Only check Unigene column
# 2) if nchar > 15 --> Check gene description, synonyms, and full name
# 3) If entirely numbers --> check LocusLinc, Gene ID, updated Gene ID
# 4) If not in groups 1-3 -->  check Accession, Symbol, Synonyms, description, Full name
#
# Inputs: myID: a string, a composite sequence identifier (e.g. "3342", "GAPDH", "INTEGRIN BETA")
#
# Outputs: geneSymbol: a vector containing the gene symbol(s) a composite sequence identifier maps to
mapOneGene <- function(myID){
  geneSymbol = c()
  
  ### 1 ### If contains "HS.", then it's a Unigene code and we can 
  if(grepl(pattern = "HS\\.", x = myID, ignore.case = T)){
    if(myID %in% heebo$unigene){
      geneSymbol = unique(c(geneSymbol, heebo$Symbol.v12[which(heebo$unigene == myID)]))
      
      # If it's HS., then I'm pretty sure it's unigene, and we don't have to check other columns
      return(geneSymbol)
    } else{
      return(geneSymbol)
    }
  }
  
  
  ### 2 ### If nchar > 15, then it's in one of three columns: 
  # Synonyms, Description, Full.name
  if(nchar(myID) > 15){
    # check if in Synonyms
    if(myID %in% heebo$Synonyms.v12){
      geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$Synonyms.v12 == myID)])
    }
    # Check if Description
    if(myID %in% heebo$description.name.v12){
      geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$description.name.v12 == myID)])
    }
    # Check if Full.name
    if(myID %in% heebo$Full.name.from.nomenclature.authority.v12){
      geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$Full.name.from.nomenclature.authority.v12 == myID)])
    }
    
    # Return geneSymbol, becuase we've checked all three columns where nchar > 15
    return(unique(geneSymbol))
  }
  
  
  ### 3 ### If a number, then check the number columns
  # Locus link ID, gene ID, updated gene ID
  if(suppressWarnings(!is.na(as.numeric(myID)))){
    
    # Check if in Link Locus ID
    if(myID %in% as.character(heebo$LocusLink_ID)){
      geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$LocusLink_ID == myID)])
    }
    # Check if in Gene ID
    if(myID %in% as.character(heebo$GeneID.v12)){
      geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$GeneID.v12 == myID)])
    }
    # Check if in updated gene ID
    if(myID %in% as.character(heebo$Updated_GeneID)){
      geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$Updated_GeneID == myID)])
    }
    
    return(unique(geneSymbol))
  }
  
  ### 4) Remaining chr
  # 15 characters or less, combination of numbers and letters, not a Unigene ID
  # c("Accession", "Symbol.v12", "Synonyms.v12", "description.name.v12", "Full.name.from.nomenclature")
  
  # check if in Accession
  if(myID %in% heebo$Accession){
    geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$Accession == myID)])
  }
  
  # check if in Symbol
  if(myID %in% heebo$Symbol.v12){
    geneSymbol = c(geneSymbol, myID) # clear case, just copy in the symbol
  }
  
  # check if in Synonyms
  if(myID %in% heebo$Synonyms.v12){
    geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$Synonyms.v12 == myID)])
  }
  # Check if Description
  if(myID %in% heebo$description.name.v12){
    geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$description.name.v12 == myID)])
  }
  # Check if Full.name
  if(myID %in% heebo$Full.name.from.nomenclature.authority.v12){
    geneSymbol = c(geneSymbol, heebo$Symbol.v12[which(heebo$Full.name.from.nomenclature.authority.v12 == myID)])
  }
  
  return(unique(geneSymbol))
  
}

#..................................................#
##### 5) Map Composite Seq IDs to Gene Symbols #####
#..................................................#
# Purpose: use mapOneGene() function to map all comp seq identifiers to gene symbols
#          And, explore results
#
# Results: 
# 78% of compIDs map to one gene symbol
# 22% of compIDs map to nothing
#     - LOCs and orfs
#     - Defunct gene codes, IDs, and symbols
# < 1% of compIDs map to more than onen gene
#     - LOCs and orfs
#     - Gene families


# Conclusion: I'm satisfied with this level of mapping!

# Run code!
timestamp()
mappedCompID = lapply(compID, mapOneGene)
timestamp()


### Explore mappedcompID to see how well it worked
length(mappedCompID)



## Loop through mappedCompIDs to find the ones that failed to map 
# and the ones that mapped to more than one gene
lengthCI = c() # number of genes mapped to each compID
isNullCI = c() # compIDs that failed to map are listed as NULL, which ones are NULL
for(i in 1:length(mappedCompID)){
  lengthCI[i] = length(mappedCompID[[i]])
  isNullCI[i] = is.null(mappedCompID[[i]])
}

sum(isNullCI) # 11290
table(lengthCI) # max nGenes is  22, min is 0
round(table(lengthCI)/length(mappedCompID), 2) # 21% is 0, 79% has one gene, a handful has many

# What's up with the ones that fail to map?
# they seem to be ORFs and transcribed locuses, not well established genes
# I'm okay with losing them
# A lot of long descriptions and defunct Unigene codes
compID[isNullCI][1:100]


# What's up with the ones that map multiple times?

# Maps to 22 protocadheren family members
compID[which(lengthCI ==22)]
mappedCompID[which(lengthCI ==22)]

table(lengthCI)
n =2
compID[which(lengthCI ==n)][1:5]
mappedCompID[which(lengthCI ==n)][1:5]


### Observations:
# 79% of CompIDs map to only 1 gene
# 21% of CompIDs fail to map, I'm okay with them
# A relatively tiny fraction of compIDs map to 2-22 gene symbols
# This mostly seems to be due to putative sequences like LOC and orfs, but
# is also due to gene families and alias gene symbols
# So, I'll list all of them, like a probe mapping to mult genes

names(mappedCompID) = compID
#save(mappedCompID, file = "exploring4_plots/mappedCompID.RData")


#...............................................................#
##### 6) Reformat mappedCompID to better match $keys format #####
#...............................................................#
# Purpose: mappedCompID is a list now, but I need a vector
# compIDs now store multiple genes as a vector, but I need it as a long string
# e.g. c("GAPDH", "ACTB", "IL4") --> "GAPDH, ACTB, IL4"


# combineVector takes a character vector, and pastes
# it together into one string
#
# e.g. c("GAPDH", "ACTB", "IL4") --> "GAPDH, ACTB, IL4"
#
# Inputs: myEntry: a character vector 
# Outputs: comboEntry: a string 
combineVector = function(myEntry){
  if(is.null(myEntry)){
    return(NA)
  }
  
  if(length(myEntry) == 1){
    return(myEntry)
  }
  
  comboEntry = paste(myEntry, collapse = ",")
  return(comboEntry)
}



# Convert mappedCompID from an unwieldly list to a vector 
mappedCompID = sapply(mappedCompID, combineVector)
names(mappedCompID) = compID
save(mappedCompID, file = "mappedCompID.RData")

#.....................#
##### 7) Fix Keys #####
#.....................#
# Purpose: Convert keys for gse11235 from Composite Sequence Identifiers
# to gene Symbols


# fixKeys uses mappedCompID to convert the 
# Composite Sequence Identifiers in the keys of gse11235
# to gene symbols
#
# Inputs: gem: a Khatri lab Dataset object with gse11235 CompIDs as keys
# Outputs: gem: same as input, except keys have been converted to gene symbols
fixKeys = function(gem){
  newKeys = mappedCompID[gem$keys]
  names(newKeys) = rownames(gem$expr)
  gem$keys = newKeys
  print(checkDataObject(gem, "Dataset"))
  return(gem)
}


# Fix keys on all datasets
gse11235_gpl6765_fix = fixKeys(gse11235_gpl6765)
gse11235_gpl6766_fix = fixKeys(gse11235_gpl6766)
gse11235_gpl6767_fix = fixKeys(gse11235_gpl6767)
gse11235_gpl6768_fix = fixKeys(gse11235_gpl6768)


# All of the platforms have the same gene expr now!
mismatch2 = c()
for(myProbe in allProbeNubs){
  mv5 = gse11235_gpl6765_fix$keys[myProbe]
  mv6 = gse11235_gpl6766_fix$keys[myProbe]
  mv7 = gse11235_gpl6767_fix$keys[myProbe]
  mv8 = gse11235_gpl6768_fix$keys[myProbe]
  if(length(unique(na.omit(c(mv5, mv6, mv7, mv8)))) > 1){
    mismatch2 = c(mismatch2, myProbe)
  }
}

length(mismatch2) # 24,060 mismatches
length(mismatch) - length(mismatch2) # 18667 fewer mismatches

# printProbe prints out the compID for a single probe in each gpl
printProbe_fix <- function(probeID){
  print(paste("gpl_5", "     ", gse11235_gpl6765_fix$keys[probeID]))
  print(paste("gpl_6", "     ", gse11235_gpl6766_fix$keys[probeID]))
  print(paste("gpl_7", "     ", gse11235_gpl6767_fix$keys[probeID]))
  print(paste("gpl_8", "     ", gse11235_gpl6768_fix$keys[probeID]))
}

# This one is fixed!
printProbe(mismatch[50])
printProbe_fix(mismatch[50])

# gpl_7 maps to SAA4, but others map to SAA2
# it has a Unigene code and that's what the code maps to now. 
printProbe(mismatch2[10])
printProbe_fix(mismatch2[10])

# I trust my probe mapping

#........................................#
##### Sanity Check for Probe Mapping #####
#........................................#
# Purpose: 
# Check how well XIST will work as a sanity check gene
# Is its probes given different lables on different arrays? 

# Answer:
# all 4 batches have probes labled by gene symbol, "XIST"
# But, they vary in whether they have Gene ID, Gene Name, or Unigene code

xistIDs = names(mappedCompID[which(mappedCompID == "XIST")])

for(myName in xistIDs){
  print(paste("gpl_5:     ", paste(gse11235_gpl6765$keys[which(gse11235_gpl6765$keys %in% xistIDs)])))
}

print(paste("gpl_5:     ", cat(unique(gse11235_gpl6765$keys[which(gse11235_gpl6765$keys %in% xistIDs)]), sep = ", ")))
print("     ")
print(paste("gpl_6:     ", cat(unique(gse11235_gpl6766$keys[which(gse11235_gpl6766$keys %in% xistIDs)]), sep = ", ")))
print("     ")
print(paste("gpl_7:     ", cat(unique(gse11235_gpl6767$keys[which(gse11235_gpl6767$keys %in% xistIDs)]), sep = ", ")))
print("     ")
print(paste("gpl_8:     ", cat(unique(gse11235_gpl6768$keys[which(gse11235_gpl6768$keys %in% xistIDs)]), sep = ", ")))
print("     ")

#......................................................#
##### 8) Combine sub-datasets inton one giant Expr #####
#......................................................#

# Purpose: gse11235 does not have any probe IDs, just numbers
# The same probe number in one batch does not necessarily map 
# to the same probe number on a different batch
# Though, there's more consistency between gpl_5+gpl_8
# as wellas gpl_6+gpl_7



# prepForCombo takes a dataset object from gse11235 and converts
# it from probe level to gene expression level
#
# Inputs: 
# gem: a gse11235 Dataset object with probe level expression
#
# Outputs: 
# gem: a gse11235 Dataset object with gene level expression

prepForCombo <- function(gem){
  # Get list of all available genes
  allGenes = gem$keys
  names(allGenes) = NULL
  allGenes= unique(allGenes)
  
  # make expr's rownames probe IDs and not numbers
  rownames(gem$expr) = paste("probe", rownames(gem$expr), sep = "")
  names(gem$keys) = rownames(gem$expr)
  
  # Put expr in gene level
  exprG = as.matrix(getSampleLevelGeneData(gem, allGenes))
  
  # incorporate exprG into gem
  gem$expr_probe = gem$expr
  gem$expr_gene = exprG
  gem$expr = exprG
  
  # tweak keys
  gem$keys_probe = gem$keys
  gem$keys_gene = rownames(exprG)
  names(gem$keys_gene) = rownames(exprG)
  gem$keys = gem$keys_gene
  
  print(checkDataObject(gem, "Dataset"))
  
  return(gem)
}


# Create gene level Dataset objects
gse11235_gpl6765G = prepForCombo(gse11235_gpl6765_fix)
gse11235_gpl6766G = prepForCombo(gse11235_gpl6766_fix)
gse11235_gpl6767G = prepForCombo(gse11235_gpl6767_fix)
gse11235_gpl6768G = prepForCombo(gse11235_gpl6768_fix)


# Combine them into one giant expr
library(dplyr)

# prepExpr reformats expr to prepare for using a full_join
# to combine each expr into one matrix
# aka transposes it, and makes sampleID be its own column
#
# Inputs: gem: Dataset object
#
# Outputs: reformatted expr
prepExpr = function(gem){
  expr = as.data.frame(t(gem$expr))
  expr$sampleIDs = rownames(expr)
  
  return(expr)
}

#Reformat each expr to prepare for joining
expr5 = prepExpr(gse11235_gpl6765G)
expr6 = prepExpr(gse11235_gpl6766G)
expr7 = prepExpr(gse11235_gpl6767G)
expr8 = prepExpr(gse11235_gpl6768G)


# Create giant expr
giantExpr = expr5 %>% 
  full_join(expr6) %>%
  full_join(expr7) %>%
  full_join(expr8)


# Yes, giantExpr has the expected dimensions!
dim(giantExpr) # 34 rows, 14420 columns
length(unique(c(colnames(expr5), colnames(expr6), colnames(expr7), colnames(expr8)))) # 14420 unique genes
sum(c(nrow(expr5), nrow(expr6), nrow(expr7), nrow(expr8))) # 34 samples

# Reformat giantExpr into expr format
rownames(giantExpr) = giantExpr$sampleIDs
giantExpr$sampleIDs = NULL
giantExpr = t(as.matrix(giantExpr))

#....................................................#
##### 9) Create a proper gse11235 Dataset object #####
#....................................................#
# Purpose: Take the giant Expr and make it into a Dataset object
gse11235_raw = gse11235


# Create Pheno
giantPheno = rbind(gse11235_gpl6765G$pheno, gse11235_gpl6766G$pheno, gse11235_gpl6767G$pheno, gse11235_gpl6768G$pheno)
#giantPheno = cleanUpPheno(giantPheno, removeChar = T)


# Unscramble Disease_STate characteristics and Genotypes columns 
for(i in 1:nrow(giantPheno)){
  if(!grepl(pattern = "Disease State", x = giantPheno$characteristics_ch2[[i]])){
    # Identify the genotype and diseaseStates in wrong columns
    genotype = giantPheno$characteristics_ch2[[i]]
    diseaseState = giantPheno$characteristics_ch2.1[[i]]
    
    # Put them in their correct columns
    giantPheno$characteristics_ch2[[i]] = diseaseState
    giantPheno$characteristics_ch2.1[[i]] = genotype
  }
}

# Make giantPheno pretty 
giantPheno = cleanUpPheno(giantPheno, removeChar = T)


# Reformat age
age = strsplitVector(giantPheno$Measurement, mySplit = " ", whichPart = 2)
age = gsub(pattern = ";", replacement = ".", x = age)
age = as.numeric(age)
giantPheno$age = age
giantPheno$Measurement = NULL


# Create a more informative column of both genetalia and Chr stats
chr = giantPheno$Genotype
chr[which(chr == "46, XX")] = "XX"
chr[which(chr == "46, XY")] = "XY"
chr = factor(x = chr, levels = c("XX", "XY", "45, X0/46, XY"))
giantPheno$chromosomes = chr


# Tidy up Phenotype
pty = as.factor(giantPheno$Phenotype)
pty = factor(x = pty, levels = c("normal female", "P1", "P1-2", "P2", "P4", "normal male"))

giantPheno$Phenotype = pty

# Combine Genotype and Phenotype information
genePhen = factor(paste(chr, pty), levels = c("XX normal female", 
                                              "XY normal female", 
                                              "XX P1", 
                                              "XY P1-2",
                                              "XY P2",
                                              "XX P4", 
                                              "45, X0/46, XY P4",
                                              "XY P4",
                                              "XY normal male" ))
giantPheno$genePhen = genePhen

# Create a caseControl column, 
# to differentiate controls from samples of disorders of sexual development
caseControl = as.character(giantPheno$Disease_State)
caseControl = replace(caseControl, which(caseControl == "female control"), "Control")
caseControl = replace(caseControl, which(caseControl == "male control"), "Control")
caseControl = replace(caseControl, which(caseControl != "Control"), "DSD")
names(caseControl) = rownames(giantPheno)
giantPheno$caseControl = as.factor(caseControl)


# Create keys
giantKeys = rownames(giantExpr)
names(giantKeys) = rownames(giantExpr)

# Create class
giantClass = rep(0, nrow(giantPheno))
names(giantClass) = rownames(giantPheno)

# Create Dataset Object
gse11235 = list(expr = giantExpr,
                pheno = giantPheno,
                keys = giantKeys,
                class = giantClass, 
                formattedName = "GSE11235", 
                comment = "Composite Seq Identifiers mapped to Gene Symbols using Stanford Functional Genomics Facility HEEBO annotations")

checkDataObject(gse11235, "Dataset")

#save(gse11235, file = "gse11235_preCOCONUT.RData")

#...................................#
##### 10) Examine Batch Effects #####
#...................................#
# Question: is there batch effect across the platforms? 
# Answer: Yes, when I look at all genes in controls alone, they separate by platform

# Conclusion:
# It's worth it to ComBat normalize samples

# create data frames containing the gene expression of 
# only the genes that have no NA's (roughly 50% of genes)
myData = as.data.frame(t(na.omit(gse11235$expr)))
myDataControls = myData[which(gse11235$pheno$caseControl == "Control"),]

# Create pdf of plots looking at the relationship of the
# batch effect to different variables. 
# The variance in samples is best explained by platform (in controls)
pdf("plots/gse11235_preComBat.pdf", width = 8)
boxplot(gse11235$expr, col = gse11235$pheno$platform_id, main = "Pre-ComBat GSE11235")
autoplot(prcomp(myData, scale. = T), data = gse11235$pheno, colour = "platform_id", size = 5) + ggtitle("Pre-ComBat GSE11235 - All Samples")
autoplot(prcomp(myDataControls, scale. = T), data = gse11235$pheno[which(gse11235$pheno$caseControl == "Control"),],
         colour = "platform_id",  shape = "chromosomes",size = 5) + ggtitle("Pre-ComBat GSE11235 - Controls")
autoplot(prcomp(myData, scale. = T), data = gse11235$pheno, colour = "Phenotype", size = 5)
autoplot(prcomp(myData, scale. = T), data = gse11235$pheno, colour = "Genotype", size = 5)
autoplot(prcomp(myData, scale. = T), data = gse11235$pheno, colour = "age", size = 5)
autoplot(prcomp(myData, scale. = T), data = gse11235$pheno, colour = "Disease_State", size = 5)
dev.off()


#..........................................#
##### 11) Prepare gse11235 for ComBat  #####
#..........................................#

expr = gse11235$expr

# Remove genes missing from >50% samples
naRows <- rowSums(is.na(expr))
expr = expr[which(naRows <(0.5*ncol(expr))),]

# Remove genes with low variance
geneVariance <- apply(expr, 1, var, na.rm=T)

# Find the zero variance genes
# I think there's an NA in zeroVarianceGenes because one gene has variance of NA (only one sample has it measured)
zeroVarianceGenes <- names(geneVariance)[geneVariance < 0.1 | geneVariance > 100]

length(zeroVarianceGenes) # 13
expr = expr[which(!(rownames(expr) %in% zeroVarianceGenes)),]

# replaceIndividualNA takes a vector of values and replaces
# NA's or NaN's with zeros or with the mean value
# If you're doing ComBat, use zero imputation
# 
# Inputs: x, a vector
#         method, a string, either "zero" or "mean", defaults to "zero"
#
# Outputs: x, that same vector, but with any NA's replaced with mean(x)
replaceIndividualNA <- function(x, method = "zero"){
  # If no NA's just return x
  if(sum(is.na(x)) == 0){
    return(x)
  }
  
  # If, meethod = "mean", then do mean imputation
  # replace NAs with mean value
  if(method == "mean"){
    xMean = mean(na.omit(x))
    x[which(is.na(x))] = xMean
  }
  
  # Default is zero imputation
  # Replaces NAs with zeros
  x[which(is.na(x))] = 0 # zero imputation
  return(x)
}


# It works!
expr = t(apply(expr, 1, replaceIndividualNA))

#....................#
##### 12) Combat #####
#....................#
# Vector separating samples by batch
batch = as.numeric(gse11235$pheno$platform_id)

# Run ComBat
expr_combat = ComBat(expr, batch)

# Create a Dataset object with combat expr
gse11235_combat = gse11235
gse11235_combat$expr_preComBat = gse11235$expr
gse11235_combat$expr = expr_combat
gse11235_combat$keys = gse11235_combat$keys[rownames(gse11235_combat$expr)]
checkDataObject(gse11235_combat, "Dataset")


# Sanity check plots for combat!
myData = as.data.frame(t(na.omit(gse11235_combat$expr)))
myDataControls = myData[which(gse11235_combat$pheno$caseControl == "Control"),]

pdf("plots/gse11235_postComBat.pdf")
autoplot(prcomp(myData, scale. = T), data = gse11235_combat$pheno, colour = "platform_id", size = 5) + ggtitle("Combat GSE11235 - all genes")
autoplot(prcomp(myDataControls, scale. = T), data = gse11235_combat$pheno[which(gse11235$pheno$caseControl == "Control"),],
         colour = "platform_id",  shape = "chromosomes",size = 5) + ggtitle("COMBAT GSE11235 - Controls")
quickViolin(gse11235_combat, myGene = "RPS4Y1", "chromosomes")
quickViolin(gse11235_combat, myGene = "XIST", "chromosomes")
dev.off()

#...............................................................#
##### 13) Create Datasets with only controls and XY-females #####
#...............................................................#
# Purpose: 
# 1) Subset dataset to only have controls and XY-females
# 2) Look at PCA to see how samples separate
#      - Do XY males and XX females separate? 
#             -Poorly, but some signal
#      - Do the XY females cluster together? Is baby an outlier?
#             - Baby is not outlier, clusters with others
#             - P-019 is an outlier: 16yo gonadal dysgen. on estrogen

samplesOfInterest = c("female control", "male control", "gonadal dysgenesis", "P450scc deficiency")
samples_XXf_XYf_XYm = gse11235_combat$pheno$Disease_State %in% samplesOfInterest
pheno_XXf_XYf_XYm = gse11235_combat$pheno[samples_XXf_XYf_XYm,]
gse11235_XXf_XYf_XYm = subsetGEMFromPheno(gse11235_combat, pheno_XXf_XYf_XYm)

checkDataObject(gse11235_XXf_XYf_XYm, "Dataset") # True


# Put NAs back in
gse11235_XXf_XYf_XYm$expr_preComBat = gse11235_XXf_XYf_XYm$expr_preComBat[rownames(gse11235_XXf_XYf_XYm$expr), colnames(gse11235_XXf_XYf_XYm$expr)]
gse11235_XXf_XYf_XYm$expr_postCB = gse11235_XXf_XYf_XYm$expr
gse11235_XXf_XYf_XYm$expr_postCB_NA = gse11235_XXf_XYf_XYm$expr_postCB
gse11235_XXf_XYf_XYm$expr_postCB_NA[which(is.na(gse11235_XXf_XYf_XYm$expr_preComBat))] = NA

all(is.na(gse11235_XXf_XYf_XYm$expr_preComBat) == is.na(gse11235_XXf_XYf_XYm$expr_postCB_NA)) # True

# removeNAGenes
# Removes genes whose NA count exceeds my inclusion criteria
#
# Inclusion criteria: 
#   1) Measured in at least 3 XY-females
#
# Inputs: 
#   - Dataset object
#datasetObj = gse11235_XXf_XYf_XYm
#myGroupCol = "genePhen"
#myGroup = "XY normal female"
#nNA = 2
removeNAGenes = function(datasetObj,myGroupCol, myGroup, nNA){
  groupSamples = (datasetObj$pheno[[myGroupCol]] == myGroup)
  
  sampExpr = datasetObj$expr_postCB_NA[,groupSamples]
  geneNAs = apply(sampExpr, 1, function(x) sum(is.na(x)))
  geneOk = (geneNAs <= nNA)
  
  # Remove genes with too many NAs from datasetObj
  datasetObj$expr_postCB_NA = datasetObj$expr_postCB_NA[geneOk,]
  datasetObj$expr = datasetObj$expr_postCB_NA
  datasetObj$keys = datasetObj$keys[geneOk]
  
  print(checkDataObject(datasetObj, "Dataset"))
  return(datasetObj)
}


gse11235_XXf_XYf_XYm = removeNAGenes(gse11235_XXf_XYf_XYm, "genePhen", "XY normal female", 1)


### PCA plots
library(ggfortify) # package for PCA plots

sampleID = as.character(gse11235_XXf_XYf_XYm$pheno$source_name_ch2)
sampleID[which(gse11235_XXf_XYf_XYm$pheno$genePhen == "XY normal male")] = "normal male"
sampleID[which(gse11235_XXf_XYf_XYm$pheno$genePhen == "XX normal female")] = "normal female"
gse11235_XXf_XYf_XYm$pheno$groupSampleID = sampleID

myData = as.data.frame(t(na.omit(gse11235_XXf_XYf_XYm$expr)))
myDataiSEXS =as.data.frame(t(na.omit(gse11235_XXf_XYf_XYm$expr[iSEXS[which(iSEXS %in% rownames(gse11235_XXf_XYf_XYm$expr))],])))

pdf("plots/gse11235_postComBat2.pdf", width=8)
autoplot(prcomp(myData, scale. = T), data = gse11235_XXf_XYf_XYm$pheno, colour = "groupSampleID", size = 5) + 
  ggtitle("gse11235  - all genes")
autoplot(prcomp(myDataiSEXS, scale. = T), data = gse11235_XXf_XYf_XYm$pheno, colour = "groupSampleID", size = 5) + 
  ggtitle("gse11235 - iSEXS (70 genes w/o NA)")
autoplot(prcomp(myDataiSEXS, scale. = T), data = gse11235_XXf_XYf_XYm$pheno, colour = "platform_id", size = 5) + 
  ggtitle("gse11235 - iSEXS (70 genes w/o NA)")
dev.off()


### Subset to have only one type of control 
pheno_XXf_XYf = subset(gse11235_XXf_XYf_XYm$pheno, genePhen %in% c("XX normal female", "XY normal female"))
nrow(pheno_XXf_XYf)
gse11235_XXf_XYf = subsetGEMFromPheno(gse11235_XXf_XYf_XYm, pheno_XXf_XYf)
checkDataObject(gse11235_XXf_XYf, "Dataset")

pheno_XYm_XYf = subset(gse11235_XXf_XYf_XYm$pheno, genePhen %in% c("XY normal male", "XY normal female"))
nrow(pheno_XYm_XYf)
gse11235_XYm_XYf = subsetGEMFromPheno(gse11235_XXf_XYf_XYm, pheno_XYm_XYf)
checkDataObject(gse11235_XYm_XYf, "Dataset")

#.........................................#
##### 14) Number of NAs in each group #####
#.........................................#
# Question: Are there sufficient numbers of NAs in female controls and male controls to worry about?
# Answer: nope!
#   - Only 3 genes are missing from 4 samples, and only 1 gene is missing from 5
#   - I don't think that's enough to worry about


expriSEXS = gse11235_XXf_XYf_XYm$expr[iSEXS[which(iSEXS %in% rownames(gse11235_XXf_XYf_XYm$expr))],]
dim(expriSEXS) # 100 genes

naMatrix = matrix(nrow = nrow(expriSEXS), ncol = 3)
rownames(naMatrix) = rownames(expriSEXS)
colnames(naMatrix) = c("XX_female", "XY_male", "XY_female")

for(i in 1:nrow(naMatrix)){
  myGene = expriSEXS[i,]
  
  XXf = myGene[which(gse11235_XXf_XYf_XYm$pheno$genePhen == "XX normal female")]
  XYm = myGene[which(gse11235_XXf_XYf_XYm$pheno$genePhen == "XY normal male")]
  XYf = myGene[which(gse11235_XXf_XYf_XYm$pheno$genePhen == "XY normal female")]
  
  
  naMatrix[i,1] = sum(is.na(XXf))
  naMatrix[i,2] = sum(is.na(XYm))
  naMatrix[i,3] = sum(is.na(XYf))
}

table(naMatrix)


#...............#
###### Save #####
#...............#

save(gse11235_raw,gse11235, gse11235_combat, gse11235_XXf_XYf_XYm, gse11235_XXf_XYf, gse11235_XYm_XYf, file= "gse11235.RData")






