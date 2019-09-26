# General functions for Gene Expression Microarray and their analysis

# Purpose: 
# Make one function that can extract gene expression from expr
# whether or not there's probes that map to multiple genes

library(stringr)
library(MetaIntegrator)
library(preprocessCore)
require(samr)

#................................#
##### getSampleLevelGeneData #####
#................................#
# getSampleLevelGeneData first expands expr so that probes that map to multiple 
# genes are duplicated, and assigned to each gene individually. This allows
# us top pick up expression from the AMY family that is only measured on one
# probe that picks up expression of multiple family members
#
# Inputs: 
#   gem: a Khatri lab Dataset object
#   genes: a character vector of gene symbols: e.g. c("GAPDH", "SMYD3", "SETD6")

getSampleLevelGeneData <- function(gem, genes){
  gem = gem_expand.df(gem)
  geneExpr= extractDataFromGEM(gem, genes)
  return(geneExpr)
}


# gem_expand.df takes a DAtaset object (gem) and expands probes that map to multiple 
# genes by duplicating those rows within expr 
#
# duplicate probes are separated by a comma in the default settings
#
# Inputs: 
#   gem: a Dataset object
#   key.name: the name of the column within df, where the genes for each probe will be stored
#   keys.sep: the character that separates gene symbols when a probe maps to more than one
#
# Outputs: 
#   gem: a modified Dataset object where the probes that map to multiple genes are duplicated
gem_expand.df <- function (gem, key.name = "keys", keys.sep = ",") {
  # I added this part, so that you could directly input a gem
  df = as.data.frame(gem$expr)
  df$keys = gem$keys
  
  # Purvesh's code: 
  ## every row in df may have multiple identities
  ## these identities are defined in keys, with individual keys separated by keys.sep
  ## this function expands every row that has 1:n mapping to n x 1:1 rows
  
  keys <- as.character(df[ ,key.name])
  skey <- strsplit( keys, split = keys.sep)
  df   <- df[rep(1:nrow(df), sapply(skey, length)), ]
  df[, key.name] <- unlist(skey)
  
  ## Convert it back into Dataset format
  # Create new keys
  gem$keys = df$keys
  names(gem$keys) = rownames(df)
  
  # create new expr
  expr = df
  expr$keys = NULL
  gem$expr = as.matrix(expr)
  return(gem)
}

# extractDataFromGEM takes a vector of genes, and returns a matrix
# of those genes' expression in that dataset
# I think Purvesh originally wrote this function
#
# Only accepts probes that maps to a single gene
# If you have probes that map to multiple genes, use gem_expand.df first
# and then use extractDataFromGEM
#
# Inputs: 
#   gem: a Khatri lab Dataset object
#   genes: a character vector of gene symbols (e.g. c("SETD6", "EZH2", "SMYD3"))
#
# Returns: 
#   
extractDataFromGEM <- function(gem, genes) {
  tempExprs = NULL
  #junk = apply(as.matrix(genes), 1, function(x, keys) which(keys == x), keys=gem$keys)
  junk = lapply(as.matrix(genes), function(x, keys) which(keys == x), keys=gem$keys)
  for(j in 1:length(junk)) {
    if(length(junk[[j]]) == 0) {
      next
    }
    temp = gem$expr[junk[[j]],]
    if(!is.vector(temp)) {
      temp = t(as.matrix(colMeans(temp, na.rm = T) ))
    } else {
      temp = t(as.matrix(temp))
    }
    rownames(temp) = genes[j]
    tempExprs = rbind(tempExprs, temp)
  }
  tempExprs = data.frame(tempExprs)
  #tempExprs$ID = rownames(tempExprs)
  return(tempExprs)
}


#................................#
##### Preprocessing Datasets #####
#................................#

# subsetGEMFromPheno
# Takes a Dataset object (gem) and subsets to suit an already subsetted pheno file
#
# Inputs: gem: Dataset object
#         newPheno: dataframe, a subsetted version of gem$pheno
#
# Outputs: newGem: a new Dataset object that only contains the samples within newPheno
subsetGEMFromPheno = function(gem, newPheno) {
  newGem = gem
  newGem$rawPheno = gem$pheno
  
  # Reduce samples in pheno, expr and class
  newGem$pheno = newPheno
  newGem$expr = gem$expr[,rownames(newPheno)]
  newGem$class = gem$class[rownames(newPheno)]
  
  # Print that it worked!
  print(checkDataObject(newGem, "Dataset"))
  
  return(newGem)
}


# cleanUpPheno takes a pheno dataframe, as given by GEOquery and cleans it up by
# 1) Removing a hard coded list of columns that are always useless
# 2) Renaming the "characteristics" columns by the characteristic they contain
#
# Inputs: 
# pheno: dataframe of phenotype info, rows are samples, cols are different categories
# removeChar: boolean, defaults to TRUE. Whether or not you want to do 2) renaming characteritics
#
# Outputs: 
# pheno: same as input, except better looking
cleanUpPheno = function(pheno, removeChar) {
  ### Remove unwanted columns
  unwantedColumns <- c("geo_accession", "status", "submission_date", "last_pudate_date", "type", "channel_count", 
                       "molecule_ch1", "extract_protocol_ch1", "label_ch1", "label_protocol_ch1", "taxid_ch1", 
                       "hyb_protocol", "scan_protocol", "data_processing", "contact_name", "contact_email", 
                       "contact_institute", "contact_address", "contact_city", "contact_zip/postal_code", 
                       "contact_country", "supplementary_file", "data_row_count", "last_update_date", "contact_state", 
                       "contact_fax", "contact_phone", "channel_count", "contact_web_link", "contact_department", 
                       "molecule_ch2", "extract_protocol_ch2", "label_ch2", "label_protocol_ch2", "taxid_ch2",
                       "contact_laboratory")
  
  
  for (thisRow in unwantedColumns) {
    pheno[thisRow] = NULL
  }
  
  
  ### Clean up "characteristics" column
  if (removeChar) {
    for (thisColName in colnames(pheno)){
      # If it's a "characteritics_ch1.1" sort of column
      if (strsplit(thisColName, "_")[[1]][[1]] == "characteristics") {
        thisCol = as.character(pheno[[thisColName]])
        
        # If the first entry of this column doesn't have ":",
        # Then it doesn't have the "sex: female" structure, so skip this column
        if (!grepl(":", thisCol[1])){
          next
        }
        
        newCol <- c()
        newColName <- str_replace_all(strsplit(as.character(thisCol[1]), ": ")[[1]][[1]], " ", "_")
        
        # Check for scrambled columns, like when sex is in age column for some samples
        allColNames = strsplitVector(myVector = thisCol, mySplit = ": ", whichPart = 1)
        if (length(unique(allColNames[allColNames != ""])) != 1) {
          warning(paste(newColName, "Characteristics columns might be scrambled"))
        }
        
        # Extract values from each sample
        newCol = strsplitVector(myVector = thisCol, mySplit = ": ", whichPart = 2)
        
        # Remove the the current column from pheno
        pheno[thisColName] = NULL
        
        ### Create a replacement column with special rules for age, sex, race columns
        # If it's an age column, then separate the numbers from the units
        if (tolower(newColName) == "age"){
          ### if it's an age column
          
          # If something's left over when you remove all numbers, then save that 
          # in a column called "ageUnits"
          ageUnits = gsub(pattern = "[0-9]", replacement = "", x = newCol)
          if (!all(ageUnits == "")){
            pheno["ageUnits"] = ageUnits
          }
          
          newCol = gsub(pattern = "[[:punct:]]", replacement = "", x = newCol) # remove punctuation
          pheno["age"] = as.numeric(gsub(pattern="[a-z]", replacement = "", x = newCol)) # remove letters and make numeric
        } else if (tolower(newColName) == "sex" || tolower(newColName)=="gender") {
          ### If it's a sex/gender column
          newCol = tolower(newCol)
          
          # If sex is annotated with just M and F, convert to "male" and "female"
          if("m" %in% unique(tolower(newCol)) || "f" %in% unique(tolower(newCol))){
            newCol = gsub(pattern = "m", replacement = "male", x = newCol)
            newCol = gsub(pattern = "f", replacement = "female", x = newCol)
          }
          pheno["sex"] = as.factor(tolower(newCol))
        } else if (tolower(newColName) =="race" || tolower(newColName) == "ethnicity") {
          ### if it's a race/ethnicity column 
          pheno["raceEthnicity"] = as.factor(tolower(newCol))
        } else {
          ### if it's not age, sex, or raceEthnicity
          # Create a new column with the correct column name, and remove the old one
          pheno[newColName] = newCol
        }
        
        
      }
    }
  }
  
  return(pheno)
}


# strsplit_emptyString returns an empty string if you
# try to split an empty string
#
# Inputs: 
#    myString: string, the string you want to split
#    mySplit: string, the pattern you're splitting around
#    whichPart: integer, whether you want the 1st, 2nd, 3rd, etc part around split
strsplit_emptyString = function(myString, mySplit, whichPart) {
  # If empty string, return empty string
  if(myString == "") {
    return("")
  }
  # if not empty, split it!
  return(strsplit(myString, split=mySplit)[[1]][[whichPart]])
}


# strsplitVector 
# String splits all members of a vector the same way
#
# Inputs: 
#       myVector: the vector you want to split, aka column from pheno
#       mySplit: the character that you'll split each item in the vector by (e.g " ")
#       whichPart: integer, do you want the first, second, etc part of the split vector
strsplitVector = function(myVector, mySplit, whichPart) {
  myVector = as.character(myVector)
  
  myNewVector = as.vector(sapply(myVector, function(x) strsplit_emptyString(myString = x, mySplit = mySplit, whichPart = whichPart)))
  return(myNewVector)
}

# Function from Stackoverflow
# Caplitalizes the first letter of each word in a string
simpleCap <- function(x) {
  if(is.na(x)){
    return(x)
  }
  
  s <- strsplit(as.character(x), " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# createClassVector takes a column from pheno, and turns it into a vector of 1's and 0's. Aka your class vector!
# Inputs: groupColumn <- a column from pheno that denotes your two classes (aka EBV+ vs EBV-)
# casesAre <- a string that tells you the exact wording withing groupColumn that tells you whether that row is a case or control. e.g. "ebv status: positive" or "ebv+"
# Generates: a table comparing the original groupColumn to the generated classVector. Eyeball it to make sure that they match up how you want them to!
# Returns: a vector of 1's and 0's. Where, 1's are cases and 0's are controls
createClassVector <- function(groupColumn, casesAre, pheno) {
  groupColumn <- as.character(groupColumn)
  classVector <- c()
  
  for (i in 1:length(groupColumn)) {
    if (groupColumn[i] == casesAre) {
      classVector[i] <- 1
    } else {
      classVector[i] <- 0
    }
    
    
  }
  print(cbind(groupColumn, classVector))
  
  names(classVector) <- rownames(pheno)
  return(classVector)
  
}


#...........................#
# Remove One Sample #
#...........................#
# Removes a single sample from a Dataset object
#
# Inputs: 
#   gem: a Dataset object
#   sampleID: a string of the sample ID, must be rownames(pheno) and colnames(expr) etc
#      e.g. "DEE4XH1N14037"
# 
# Outputs: 
#   gem: a Dataset object with that one sample removed

removeOneSample <- function(gem, sampleID) {
  n = nrow(gem$pheno)
  
  gem$expr = gem$expr[,-which(colnames(gem$expr) == sampleID)]
  gem$pheno = gem$pheno[-which(rownames(gem$pheno) == sampleID),]
  gem$class = gem$class[-which(names(gem$class) == sampleID)]
  
  # If it's a deconvoluted dataset, remove sample from 
  if("cellProp" %in% names(gem)){
    gem$cellProp = gem$cellProp[,-which(colnames(gem$cellProp) == sampleID)]
    gem$geneExpr = gem$geneExpr[,-which(colnames(gem$geneExpr) == sampleID)]
    gem$cellPropStats = gem$cellPropStats[-which(rownames(gem$cellPropStats) == sampleID),]
  }
  
  if(all(colnames(gem$expr) != rownames(gem$pheno))){
    warning(paste(gem$formattedName, "Samples in expr do not match pheno"))
  }
  
  if(!all(colnames(gem$expr) == names(gem$class))){
    warning(paste(gem$formattedName, "Samples in expr do not match class"))
  }
  
  if(n==ncol(gem$expr)){
    warning(paste(gem$formattedName, "Number of samples has not changed"))
  }
  return(gem)
}


#GEM_normalizer
#~[by Francesco 2014/10/31]
#~this function normalizes a GEM expression matrix
GEM_normalizer <- function(GEM){
  #normalize by quantile
  norm_expr <- normalize.quantiles(as.matrix(GEM$expr))
  
  #This method is based upon the concept of a quantile-quantile 
  #plot extended to n dimensions. No special allowances are made 
  #for outliers. If you make use of quantile normalization please 
  #cite Bolstad et al, Bioinformatics (2003).
  row.names(norm_expr) <- row.names(GEM$expr)
  colnames(norm_expr)  <-  colnames(GEM$expr)
  
  #return GEM
  GEM$expr <- norm_expr
  return(GEM)
}


# combineDataFrame takes all the rows in a data frame and combines
# them into one vector/row. If all the entries in a column are the
# same (e.g. "Female"), then the resulting entry is just that ("Female")
# If the entries in a column are different, they are combined and separated
# with a semicolon (e.g. "-24; 0")
#
# Inputs:
#   df: a dataframe, with 1 or more rows
# 
# Outputs: 
#   comboDF: a data frame, with only 1 row
combineDataFrame <- function(df){
  # If there's only one row in df, just return it as a vector
  if (nrow(df) == 1){
    comboDF = df # This leaves it as a one row data frame
    rownames(comboDF) = df$SUBJECTID
    
    # If "SAMPLEID" AND "CEL" are columns (like in the Synapse data)
    # Then, make their content be characters, instead of factors
    if("SAMPLEID" %in% names(comboDF)){
      comboDF[,"SAMPLEID"] = as.character(comboDF[,"SAMPLEID"])
    }
    if("CEL" %in% names(comboDF)){
      comboDF[,"CEL"] = as.character(comboDF[,"CEL"])
    }
    return(comboDF)
  }
  
  # If there's more rows in df
  # Combine them into a one row data frame
  comboDF = df[-c(1:nrow(df)),]
  for(i in 1:ncol(df)){
    if(length(unique(df[,i])) == 1){
      # If there's only one unique entry in this column (e.g. 'GENDER')
      comboDF[1,i] = unique(df[,i])
      
    } else{
      # If there's more than one unique entry in this column, (e.g. 'TIMEPOINT')
      # Paste all the different possibilities into a string
      comboDF[,i] = as.character(comboDF[,i])
      comboDF[1,i] = paste(unique(as.character(df[,i])),collapse = "; ")
    }
  }
  
  rownames(comboDF) = comboDF$SUBJECTID
  return(comboDF)
  
}


# averageSubjectSamples takes a Khatri lab Dataset object and averages
# samples from the same subjectID together
#
# Inputs: 
#   gem: a Khatri lab Dataset object
#   subjectIDCol = string, the column name of the column in pheno that contains subject ID info
#
# Outputs: 
#   newGem: a Khatri lab Dataset object that only has one entry per subject
averageSubjectSamples <- function(gem, subjectIDCol) {
  # vector of all subjectIDs
  subjectIDs = unique(as.character(gem$pheno[,subjectIDCol]))
  
  # Create a new expr matrix where columns are subjectIDs and rows are probes
  newExpr = matrix(nrow = nrow(gem$expr), ncol = length(subjectIDs), 
                   dimnames = list(rownames=rownames(gem$expr), colnames=subjectIDs))
  
  for(subject in subjectIDs){
    sampleIDs = rownames(gem$pheno)[which(gem$pheno[[subjectIDCol]] == subject)]
    
    # If there's only one sample for that subject, copy it over to newExpr
    if(length(sampleIDs) == 1){
      newExpr[,subject] = gem$expr[,sampleIDs]
    } else{
      # If there's more than one sample for that subject, average them together
      newExpr[,subject] = rowMeans(gem$expr[,sampleIDs])
    }
  }
  
  # Collapse pheno to subject level data
  #newPheno = data.frame(matrix(nrow=length(subjectIDs), ncol=ncol(gem$pheno), 
  #                            dimnames=list(rownames=subjectIDs, colnames=colnames(gem$pheno))))
  
  # Initiate an empty new data frame with the correct colnames to be newPheno
  newPheno = gem$pheno[-c(1:nrow(gem$pheno)),]
  
  # Convert factor columns in newPheno to character, so we can add new values
  for(column in colnames(newPheno)){
    if(is.factor(newPheno[[column]])){
      newPheno[[column]] = as.character(newPheno[[column]])
      gem$pheno[[column]] = as.character(gem$pheno[[column]])
    }
  }
  
  for(subject in subjectIDs){
    subjectPheno = gem$pheno[which(gem$pheno[,subjectIDCol] == subject),]
    newPheno[subject,] = combineDataFrame(subjectPheno)
  }
  
  # Establish newGem
  newGem = gem
  newGem$pheno = newPheno
  newGem$expr = newExpr
  newGem$class = rep(0, nrow(newGem$pheno))
  names(newGem$class) = rownames(newGem$pheno)
  
  return(newGem)
}

#...............................#
##### Meta-Analysis Helpers #####
#...............................#

# sigGeneCountsSummary will an object that's the result of 
# runMetaAnalysis() and tell you how many genes are present
# at different summary effect size and FDR thresholds
#
# Inputs: 
#   metaObj: Khatri Lab Meta Object
# 
# Outputs: 
#   thresholdMatrix: a matrix where columns are FDR thresholds and rows are effect size thresholds

sigGeneCountsSummary <- function(metaObj){
  # Summarize the statistics for all the genes
  allGenes = metaObj$metaAnalysis$pooledResults
  
  # Establish the interval of thresholds
  effectSizeThresh = seq(from=0.1, to = 1, by = 0.05)
  fdrThresh = seq(from = 1, to = 0.05, by = -0.05)
  
  # Initialize thresholdMatrix
  thresholdMatrix = matrix(nrow = length(effectSizeThresh), ncol = length(fdrThresh))
  rownames(thresholdMatrix) = effectSizeThresh
  colnames(thresholdMatrix) = fdrThresh
  
  # Loop through all combinations of effectSizeThresh and fdrThresh 
  # and record the number of genes in thresholdMatrix
  for(i in effectSizeThresh){
    for(j in fdrThresh){
      threshGenes = allGenes[which(allGenes$effectSize > i | allGenes$effectSize < -i),]
      threshGenes = threshGenes[which(threshGenes$effectSizeFDR < j),]
      thresholdMatrix[as.character(i),as.character(j)] = nrow(threshGenes)
    }
  }
  
  # Add es and fdr to the rownames and colnames so you remember which is which
  rownames(thresholdMatrix) = paste("es", rownames(thresholdMatrix), sep= "")
  colnames(thresholdMatrix) = paste("fdr", colnames(thresholdMatrix), sep="")
  
  return(thresholdMatrix)
}

# sigGeneStats takes a metaObj and a filter object
# and gives you stats behind the genes in that filter object
# it's meant to help you figure out how reproducible the genes
# in your signature are 
#
# Column Names:
# nTotal: total number of genes in signature
# nMiss: number of signature genes not measured in dataset
# nSignRight: number of signature genes whose dataset effect size is in the same direction as summary effect size
# nSignRight_sig: number of signature genes whose dataset effect size - standard error is greater than zero
# nSignRight_nonsig: nSignRight - nSig, aka the genes that are non-significantly in the right direction
# nSignWrong: number of signature genes whose effect size is in the wrong direction in that dataset
#
# Inputs: 
#   metaObj: Khatri lab Meta object
#   filterObj: Khatri lab filter object
#   key.split: character that splits gene symbols when multiple genes map to the same probe
#   as.proportion: whether you want everything outputed as a proportions
#
# Outputs: 
#   finalTable: matrix where rows are datasets and columns are useful stats (listed above)
sigGeneStats = function(metaObj, filterObj, key.split = ",", as.proportions = FALSE){
  # Setup gene lists 
  posGenes = filterObj$posGeneNames
  negGenes = filterObj$negGeneNames
  allGenes = c(posGenes, negGenes)
  
  # Get the total number of genes in the signature
  # repeat it so that you can put it in a table at the end
  nTotal = rep(length(allGenes), length(metaObj$originalData))
  
  # Get the number of genes missing in each dataset
  nMiss = sapply(metaObj$originalData, function(x) sum(!(allGenes %in% unlist(strsplit(x$keys, split = key.split)))))
  
  # Extract effect sizes
  posES = metaObj$metaAnalysis$datasetEffectSizes[posGenes,]
  negES = metaObj$metaAnalysis$datasetEffectSizes[negGenes,]
  
  # Get the total number of genes whose dataset effect size is 
  # in the same direction as the summary effect size
  nSignRight_pos = apply(posES, 2, function(x) sum(x > 0, na.rm = T))
  nSignRight_neg = apply(negES, 2, function(x) sum(x < 0, na.rm = T))
  nSignRight = nSignRight_pos + nSignRight_neg
  
  
  # Extract Standard Error
  posStanEr = metaObj$metaAnalysis$datasetEffectSizeStandardErrors[posGenes,]
  negStanEr = metaObj$metaAnalysis$datasetEffectSizeStandardErrors[negGenes,]
  
  # Subtract Effect Sizes by standard error
  posSig = posES - posStanEr
  posSig = apply(posSig, 2, function(x) sum(x >0, na.rm=T)) 
  negSig = negES + negStanEr
  negSig = apply(negSig, 2, function(x) sum(x < 0, na.rm=T))
  
  # Get the total number of "significant" genes
  # aka number of genes where effect size - standard error doesn't cross zero
  nSignRight_sig = posSig + negSig
  
  # Get the number of non-significant genes in the right direction
  nSignRight_nonsig = nSignRight - nSignRight_sig
  
  # Get the number of genes with teh sign wrong
  nSignWrong = nTotal - nSignRight - nMiss
  
  # Combine all the stats into one big table
  finalTable = cbind(nTotal, nMiss, nSignRight, nSignRight_sig, nSignRight_nonsig, nSignWrong)
  
  # If as.proportions == TRUE
  # then give all the values as proportions of total number of genes
  if(as.proportions == TRUE){
    # Calculate the proportion of genes measured that fall into each category
    # Remove missing genes from total genes when calculating proportion
    propTable = cbind(nSignRight, nSignRight_sig, nSignRight_nonsig, nSignWrong)
    propTable = apply(propTable, 2, function(x) round(x/(nTotal - nMiss), digits = 2))
    
    finalTable = cbind(nTotal, nMiss, propTable)
  }
  
  return(finalTable)
}

# sigGeneStats_graph creates a stacked bar plot
# showing in the datasets within a metaObj, how
# many genes from your filter (filterObj) are:
#
# 1) Missing (nMiss)
# 2) Significantly in the right direction (standard error doesn't cross 0) (nSignRight_sig)
# 3) In the right direction, but standard error crosses 0 (nSignRight_nonsig)
# 4) In the wrong direction (nSignWrong)
#
# Inputs: 
#   metaObj: a Khatri lab meta object
#   filterObj: a Khatri lab Filter Object
#
# Generates: 
#   a stacked bar plot showing what proportons of your signature genes are missing, 
#   in the right direction, nonsignificant, or in the wrong direction for each dataset
sigGeneStats_graph <- function(metaObj, filterObj, plotTitle = NULL){
  myStats = as.data.frame(sigGeneStats(metaObj, filterObj))
  
  # Remove columns that won't be in stacked bar plot
  myStats$nTotal = NULL
  myStats$nSignRight = NULL
  
  # Melt data frame to prep for graphing
  myStats$dataset = factor(rownames(myStats), levels = rownames(myStats))
  
  
  myStats = melt(myStats, id.vars = "dataset")
  
  ggplot(data = myStats, aes(x = dataset, y = value, fill = variable)) + geom_bar(stat = "identity") + 
    coord_flip() + theme_set(theme_grey(base_size = 24)) + ggtitle(plotTitle)
}

# calcES is getES renamed
# Purvesh wrote getES, it's part of the guts of MetaIntegrator
# getES calculates a bunch of basic statistics between two groups
#
# Inputs: 
#    v: numeric vector, values of the trait you measured in your cohort
#    g: binary vector (1 or 0), marks each sample as being in group 1 or group 0
#
# Outputs: 
# A list containing many things including: 
#    g = Hedges' g, aka (mean1 - mean0)/pooledStandDeviation

calcES <- function(v, g){
  ## function to calculate basic statistics for two-class comparison for a gene
  
  stopifnot( identical( length(v), length(g) ) )
  
  # Separate the values belonging to the two different classes
  # x = values for cases, y = values for controls
  # cleanNA removes NAs and inf 
  x <- cleanNA( v[ which(g==1) ] )
  y <- cleanNA( v[ which(g==0) ] )
  
  # Get the number of samples in two classes
  n1 <- length(x); n2 <- length(y)
  
  # If either class has less than two values in it
  # then return NAs for everything
  if( n1 < 2 | n2 < 2 )
    return( c(n1=NA, m1=NA, sd1=NA,
              n2=NA, m2=NA, sd2=NA,
              diff=NA, pooled.sd=NA,
              g=NA, se.g=NA) )
  
  # Get the means for both classes and the diff between the means
  m1   <- mean(x); m2 <- mean(y)
  diff <- m1 - m2
  
  # Get the standard deviations for both classes
  sd1  <- sd(x);  sd2 <- sd(y)
  
  # Calculate the pooled standard deviation
  sp   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) )
  
  # I think this is a correction factor made for small sample sizes
  cf   <- 1 - 3/( 4*(n1 + n2) - 9 )
  
  # Hedges' g = difference in means/pooled standard deviation
  # I think cf is a correction factor
  g    <- cf * diff/sp
  
  # Standard error of Hedges' g
  se.g <- sqrt( (n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94) )
  
  return( c(n1=n1, m1=m1, sd1=sd1,
            n2=n2, m2=m2, sd2=sd2,
            diff=diff, pooled.sd=sp,
            g=g, se.g=se.g) )
  
}

# Purvesh wrote cleanNA
# It removes NA and inf values from a vector
# Inputs: 
#   x: a vector, preferably, numbers
# Outputs:
#   x: same as input vector, but without NA or inf values
cleanNA <- function(x) return( x[!is.na(x) & is.finite(x) ] )


# geomMean
# Calculates the geometric mean of a vector of values#
# 
# Inputs: 
#    x: a numeric vector
#
# Outputs: 
#   geometic mean: a numeric value
geomMean = function (x, na.rm = FALSE) {
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
    }
  
  if (na.rm) {
    x <- x[!is.na(x)]
    }
  if (any(x < 0)) {
    stop("'x' contains negative value(s)")
    }
  
  ### direct method causes overflow errors, use log method instead
  ### return(prod(x)^(1/length(x)))  
  return(exp(sum(log(x))/length(x)))
}

# signatureDiffPval
# produces the p-value from a t-test of the signature score
# comparing the signature score in the Cases (class == 1)
# and controls (class == 0)
#
# Inputs: 
#   gem: Khatri lab Dataset object
#   filterObject: a Khatri lab Filter Object
#
# Outputs:
#   p-values: float, the output of t.test()$p.value
signatureDiffPval = function(gem, filterObj){
  # Calculate score for all samples
  scoreAll = calculateScore(filterObj, gem)
  
  # Subset by class
  score0 = scoreAll[gem$class == 0]
  score1 = scoreAll[gem$class == 1]
  
  # Calculate t.test
  return(t.test(score0, score1)$p.value)
}

#....................#
##### Impute Sex #####
#....................#

# imputeSex
# imputes the sex of each sample in a Dataset object by performing K means
# clustering on the genes XIST, RPS4Y1, and KDM5D
#
# These genes were chosen because of Toker 2016's publication
# on mislabeled samples 
#
# Inputs: 
#    myDataset: a Khatri lab Dataset object
#    femGenes: vecotor of gene symbols of genes higher expressed in females
#    malGenes: vector of gene symbols of genes higher expressed in males
#
# Outputs: 
#    imputedSex: vecotor indicating whether each sample is classified as "male" or "female" 

imputeSex = function(myDataset, femGenes = "XIST", malGenes = c("RPS4Y1", "KDM5D")){
  # Make sure at least one of the sex annotation genes exists in this dataset
  if (!any(c(femGenes, malGenes) %in% myDataset$keys)){
    warning("Sex classification genes not present in dataset")
    return(NULL)
  }
  
  # Get sex genes in dataset
  femGenes = femGenes[which(femGenes %in% myDataset$keys)]
  malGenes = malGenes[which(malGenes %in% myDataset$keys)]
  sexGenes = c(femGenes, malGenes)
  sexExpr = getSampleLevelGeneData(myDataset, sexGenes)
  
  # Run kmeans
  kmeansResults = kmeans(t(sexExpr), centers = 2)
  
  
  # Figure out which cluster is female
  ### 1 ### If both male and female genes measured: 
  if(length(femGenes) > 0 & length(malGenes)>0){
    # Get the cenetroid expression of female genes
    femCenter = kmeansResults$centers[,femGenes]
    if(!is.null(ncol(femCenter))){
      femCenter = rowMeans(femCenter)
    }
    
    # Get the centroid expression of male genes
    malCenter = kmeansResults$centers[,malGenes]
    if(!is.null(ncol(malCenter))){
      malCenter = rowMeans(malCenter)
    }
    
    # Get the centroid with high female expression and low male expression
    highFemCentroid = names(femCenter)[which(femCenter == max(femCenter))]
    lowMalCentroid = names(malCenter)[which(malCenter == min(malCenter))]
    
    # Make sure it agrees
    if(highFemCentroid != lowMalCentroid){
      warning("Strange centroid assignment in kmeans. Try different genes")
      return(NULL)
    }
    
    # Assign the femaleCentroid
    femaleCentroid = highFemCentroid
    
  }
  
  ### 2 ### if there are no male genes measured in the dataset
  if(length(femGenes) ==0 ){
    # Get the cenetroid expression of male genes
    malCenter = kmeansResults$centers[,malGenes]
    if(!is.null(ncol(malCenter))){
      malCenter = rowMeans(malCenter)
    }
    
    # Assign the female centroid
    femaleCentroid = names(malCenter)[which(malCenter == min(malCenter))]
  }
  
  ### 3 ### if there are no female genes measured in the dataset
  if(length(malGenes) ==0 ){
    # Get the cenetroid expression of female genes
    femCenter = kmeansResults$centers[,femGenes]
    if(!is.null(ncol(femCenter))){
      femCenter = rowMeans(femCenter)
    }
    
    # Assign the female centroid
    femaleCentroid = names(femCenter)[which(femCenter == max(femCenter))]
  }
  
  
  # Create imputed sex
  imputedSex = ifelse(kmeansResults$cluster == as.numeric(femaleCentroid), yes = "female", no = "male")
  
  return(imputedSex)
}


#.............#
##### Age #####
#.............#

# sexAgeTtest
# Tests if there's a statistical difference between
# the ages between the two classes in a dataset object
#
# Inputs: 
#   myDataset: a Khatri Lab Dataset Object
#              Must have class listed with 1's and 0's

classAgeTtest <-function(myDataset){
  # Give warning if there's no age in pheno
  if(!"age" %in% colnames(myDataset$pheno)){
    warning("age not in Pheno")
    return(NULL)
  }
  
  # Give warning if class is all 0's or all 1's
  if(length(unique(myDataset$class))==1){
    warning("Class vector must have both 1's and 0's")
    return(NULL)
  }
  
  femAge = myDataset$pheno$age[which(myDataset$class == 1)]
  malAge = myDataset$pheno$age[which(myDataset$class == 0)]
  
  tTestResults = t.test(femAge, malAge)
  return(tTestResults$p.value)
}

#..........................#
##### Plotting Helpers #####
#..........................#
# metaHeatmatrix extracts the effect sizes from the datasets in a 
# validation or discovery meta-object, but puts them in the order 
# of pooled effect sizes from discovery
#
# Inputs:
#   metaObject: Khatri lab Meta object, intended to be Discovery
#               The Meta-object that the effect sizes are ordered by
#   filterObject: Khatri lab Filter object, contains genes to be plotted
#   valObject: Khatri Lab Meta object or NULL, defaults to NULL
#              If NULL, the effect sizes from metaObject (Discovery) will be plotted
#              If a Meta object, the effect sizes from valObject will be plotted, but in 
#                  the order of the pooled effect sizes within metaObject
#   colorRange: Gives the accepted most extreme positive and negative values for effect sizes
#               aka any ES greater than 1 gets convereted to 1, so that the weaker effect sizes
#               aren't drowned out. Defaults to c(-1,1)
#   sortByTissue: boolean, defaults to F. Whether to sort datasets based on tissue type
#                 There must be a 'tissue' column in $pheno in order for it to work
#   naToZero: boolean, defaults to F. Convert NAs to 0. Use when your heatmap function
#             can't handle NAs.
#
# Outputs:
#   heatmatrix: a matrix where rows are datasets and columns are genes 
#               Contains effect sizes of the genes from filterObject from
#               the datasets in metaObject or valObject
metaHeatmatrix = function(metaObject, filterObject, valObject = NULL, colorRange = c(-1, 1), sortByTissue =F, naToZero=F) {
  
  # Get gene names from Meta object
  geneNames <- c(filterObject$posGeneNames, filterObject$negGeneNames)
  geneNames <- geneNames[geneNames %in% rownames(metaObject$metaAnalysis$pooledResults)]
  
  # Get the gene orders from Meta object
  pooledES = metaObject$metaAnalysis$pooledResults[geneNames, "effectSize"]
  names(pooledES) = geneNames
  pooledES = sort(pooledES)
  
  # Determine which meta object we want to grab ES from
  if(!is.null(valObject)){
    # Make sure valObject passes checkDataObj
    if(!checkDataObject(valObject, "Meta", "Pre-Filter")){
      return()
    }
    
    # Make valObject the object to plot
    metaObjectToPlot = valObject
  }else{
    # If valObject isn't given, then plot "Discovery" metaObject
    metaObjectToPlot = metaObject
  }
  
  
  # Extract dataset-level effect sizes from metaObjectToPlot
  # Only extracts genes that were measured in metaObjecToPlot
  geneNames_measured = geneNames[geneNames %in% rownames(metaObjectToPlot$metaAnalysis$pooledResults)]
  heatmatrix <- cbind(data.frame(Pooled = metaObjectToPlot$metaAnalysis$pooledResults[geneNames_measured, "effectSize"]), 
                      metaObjectToPlot$metaAnalysis$datasetEffectSizes[geneNames_measured,]) 
  
  # If metaObjectToPlot is missing any genes, add them back into heatmatrix
  # give all their values as NAs
  if(length(geneNames) > length(geneNames_measured)){
    heatmatrix = heatmatrix[geneNames,]
    rownames(heatmatrix) = geneNames
  }
  
  
  ### Adjust dataset/column names
  
  # Create a function to grab formatted names                                                                                                                                                                                                                                                                                  ])
  nameFun <- function(obj) {
    obj$formattedName
  }
  
  
  # Rename datasets by their formatted names
  formattedNames <- c("Pooled", sapply(metaObjectToPlot$originalData, 
                                       nameFun))
  colnames(heatmatrix) <- formattedNames
  
  
  # If sortByTissue, then rearrange use formatted names
  if(sortByTissue){
    tissue = sapply(metaObjectToPlot$originalData, function(x) x$pheno$tissue[1])
    heatmatrix = heatmatrix[,c(1,(order(tissue)+1))]
  }
  
  
  ### Rearrange genes
  # reorder genes in heatmatrix by the pooledES from Discovery (metaObj)
  heatmatrix <- as.matrix(heatmatrix[names(pooledES),])
  
  
  # Take values that are greater than 1 or less than -1
  # and set to 1 or -1 (or color range thresholds)
  heatmatrix[which(heatmatrix < colorRange[[1]])] <- colorRange[[1]]
  heatmatrix[which(heatmatrix > colorRange[[2]])] <- colorRange[[2]]
  heatmatrix <- t(heatmatrix)
  
  # Set NAs to 0
  if(naToZero){
    heatmatrix[is.na(heatmatrix)] = 0
  }
  
  # Return heatmatrix, so I can make heatmap
  if(sortByTissue){
    return(list(heatmatrix = heatmatrix, tissue=tissue))
  }
  # If not, then just return heatmatrix
  return(heatmatrix)
  
  
}

# pToStars takes a p-value and converts it to the star symbol 
#
### What the stars mean
# ns --> P > 0.05
# * --> P ≤ 0.05
# ** --> P ≤ 0.01
# *** --> P ≤ 0.001
# **** --> P ≤ 0.0001 (For the last two choices only)
#
# Inputs: 
#   - p: a numeric, a p-value (e.g. 0.003)
# Outputs:
#   - A string denoting the * symbol of significance
pToStars = function(p){
  if(p > 0.05){
    return('N.S.')
  } else if(p> 0.01){
    return('*')
  } else if(p >0.001){
    return('**')
  } else if(p > 0.0001){
    return('***')
  }else{
    return('****')
  }
}


#.......................#
##### Deconvolution #####
#.......................#
# 9/27/16 Devonvolution using Dataset objects

# Purpose: 
# Create wrapper script that will take a gem and a basis matrix
# and perform immunoStates deconvolution from the public version of MetaIntegrator

# source Francesco's deconvolution code
library(MetaIntegrator)
library(data.table)


# gemDeconvo
# Takes a Dataset object and performs deconvolution on it using
# MetaIntegrator::immunoStatesDecov
#
# It's basically a wrapper for immunoStatesDecov that creates a single
# Dataset object
#
# Inputs: 
# gem: a Khatri lab Dataset object
#
# Outputs:
# ciberGEM: a Khatri lab Dataset object, but instead of gene symbols it as cell types
gemDeconvo <- function(gem){
  # Create Meta Object
  myMetaObj = list(originalData = list(gem = gem))
  
  decovResults = MetaIntegrator::immunoStatesDecov(myMetaObj)
  
  # Create an expr from cell prop
  cellProp = decovResults$immunoStates$gem
  sampleIDs = cellProp$rn
  cellPropStats = as.data.frame(cellProp[,c("P-value", "Correlation", "RMSE")])
  rownames(cellPropStats) = sampleIDs
  cellProp = as.matrix(cellProp[,!c("P-value", "Correlation", "RMSE", "rn")])
  rownames(cellProp) = sampleIDs
  
  # Create Dataset object with cell type proportions for $expr
  decovGEM = gem
  decovGEM$geneExpr = gem$expr
  decovGEM$cellProp = t(cellProp)
  decovGEM$cellPropStats = cellPropStats
  decovGEM$expr = t(cellProp)
  decovGEM$geneKeys = gem$keys
  decovGEM$cellKeys = rownames(decovGEM$expr)
  names(decovGEM$cellKeys) = rownames(decovGEM$expr)
  decovGEM$keys = decovGEM$cellKeys
  
  # Check if it passes our Dataset rules
  checkDataObject(decovGEM, "Dataset")
  return(decovGEM)
}


# summarizeDeconvoCorr
# Summarizes the distribution of deconvolution correlations
#
# Inputs: deconvoGEM: a Dataset object created from gemDeconvo()
#
# Outputs: A vector saying how many samples are in different bins
#          also, a historgram

summarizeDeconvoCorr <- function(deconvoGEM){
  myCor = deconvoGEM$cellPropStats[,"Correlation"]
  
  corDist <- c()
  corDist["lessThan0.5"] <- sum(myCor<0.5)
  corDist["lt0.7"] <- sum(myCor < 0.7 & myCor > 0.5)
  corDist["lt0.8"] <- sum(myCor < 0.8 & myCor > 0.7)
  corDist["lt0.9"] <- sum(myCor < 0.9 & myCor > 0.8)
  corDist["lt1"] <- sum(myCor < 1 & myCor > 0.9)
  
  hist(myCor, main = deconvoGEM$formattedName)
  
  return(corDist)
  
}

# gemDeconvoQC
# Takes a Khatri lab Dataset object, runs gemDeconvo to calculate
# deconvoluted cell proportions. Then, it removes samples whose
# correlation are below a threshold set by the user
# (default is 0.7)
#
# If all samples are below the correlation threshold, 
# it returns the entire deconvoluted dataset with a warning
#
# Inputs: 
#   myDataset: Khatri lab Dataset object
#   corThresh: float, the correlation threshold that all included samples must have
gemDeconvoQC <- function(myDataset, corThresh = 0.7){
  # Run gemDeconvo
  deconvoDataset = gemDeconvo(myDataset)
  
  # Summarize the quality of the results
  summarizeDeconvoCorr(deconvoDataset)
  
  # If any of the samples have a threshold below the threshold, remove them
  if (any(deconvoDataset$cellPropStats$Correlation < corThresh)){
    
    # If all samples are below the threshold, stop everything
    if (all(deconvoDataset$cellPropStats$Correlation < corThresh)){
      print(paste("All samples below correlation threshold of", corThresh))
      return(deconvoDataset)
    }
    
    # Print the number of samples below the threshold
    print(paste(sum(deconvoDataset$cellPropStats$Correlation < corThresh), "samples below threshold"))
    
    # Get GSMs of samples with correlation below the threshold
    sampsBelowThresh = rownames(deconvoDataset$pheno)[deconvoDataset$cellPropStats$Correlation < corThresh]
    # Remove those samples one by one from dataset
    for(mySamp in sampsBelowThresh){
      deconvoDataset = removeOneSample(deconvoDataset, mySamp)
    }
    
    # Print new class balance
    print("QC Pass Sample Distribution:")
    print(table(deconvoDataset$class))
  }
  # Return Dataset
  print(checkDataObject(deconvoDataset, "Dataset"))
  return(deconvoDataset)
}


#......................#
# Cell Freq Correction #
#......................#
# corrector takes an IMmunoState Dataset Object
# and returns a Dataset object with corrected gene expression
#
# Inputs: iSObj: a Dataset Object that is outputted by my ImmunoStates function
#          # $geneExpr: the original matrix of gene expression
#          # $ expr: the matrix of immune cell frequencies
#
# Outputs: 
#     correctedObj: a Dataset object where expr is immune cell frequency corrected gene expression
corrector <- function(iSObj, exprDatasetObj){
  expMat <- na.omit(iSObj$geneExpr)
  iS <- iSObj$expr
  
  res <- lm(t(expMat) ~ 1+t(iS))
  corrected <- t(sweep(residuals(res), 2L, coef(res)[1L,], '+'))
  
  # Create corrected expression Dataset object
  correctedObj = iSObj
  correctedObj$exprUncorrected = iSObj$geneExpr
  correctedObj$exprCorrected = corrected
  correctedObj$expr = corrected
  
  correctedObj$keys = exprDatasetObj$keys[rownames(corrected)]
  
  checkDataObject(correctedObj, "Dataset")
  
  return(correctedObj)
}

# switchGeneCellprop
# Takes a deconvoluted Dataset object and reverts it back 
# to a gene expression Dataset object
#
# Inputs: 
#   deconvoDataset: a Khatri lab dataset object with cell proportions instead of gene expr
#
# Outputs: 
#   deconvoDataset: the same dataset, except it's been converted back to a 
switchGeneCellprop <- function(deconvoDataset, switchTo = "gene"){
  # If switchTo is not gene or cellProp
  # let the user know
  if(!switchTo %in% c("gene", "cellProp")){
    warning('switchTo must be cellProp or gene')
    return(deconvoDataset)
  }
  
  # If switchTo is gene, set expr and keys to have to do with gene expression
  # If it's cellProp, make expr cellproportions and keys cell type names
  if (switchTo == "gene"){
    deconvoDataset$expr = deconvoDataset$geneExpr
    deconvoDataset$keys = deconvoDataset$geneKeys
  } else if(switchTo == "cellProp"){
    deconvoDataset$expr = deconvoDataset$cellProp
    deconvoDataset$keys = deconvoDataset$cellKeys
  }
  
  # Check to make sure it passes muster
  print(checkDataObject(deconvoDataset, "Dataset"))
  return(deconvoDataset)
}

# reorderGenesByCluster
# Takes a list of genes and reorders them based off of
# hierarchical clustering of the genes within an effect 
# size matrix (myMatrix)
#
# myMatrix is intended to be a matrix of immunoStates
# effect sizes
#
# Inputs: 
#   myGenes: vector of gene symbols, aka c('GAPDH', 'XIST', 'CD40LG')
#   myMatrix: matrix of values (usually effect sizes) where row names
#             are gene symbols and include the genes listed in myGenes
#
# Outputs:
#   newGenes: vector of gene symbols, the same genes as myGenes, except
#             reordered based off of the hierarchical clustering of myGenes
reorderGenesByCluster = function(myGenes, myMatrix){
  # Remove genes that are not in myMatrix
  myGenes = myGenes[myGenes %in% rownames(myMatrix)]
  
  # If no genes are present in my Matrix
  # then give warning and return nothing
  if(length(myGenes) ==0){
    warning('myGenes not present in myMatrix')
    return()
  }
  
  # If only one gene present, return that gene
  if(length(myGenes) == 1){
    return(myGenes)
  }
  
  #Calculate distance matrix
  d = dist(myMatrix[myGenes,])
  
  # Perform hierarchical clustering
  hc = hclust(d)
  
  # Get the new order of the genes by hierarchical clustering
  newGenes = hc$labels[hc$order]
  return(newGenes)
}



# cellTypeEnrichment
# Takes a list of genes and an immunoStates matrix and tells you 
# Whether this gene list is enriched for a particular cell type
#
# Takes the effect sizes for all genes, raises them all above zero. Then,
# it takes the geometric mean of the effect sizes of all the genes. It then
# centers and scales the geometric means to convert them into z-scores. 
# The z-scores can tell you whether it's statistically significant difference
#
# Inputs:
#   - myGenes: vecotr of gene names
#   - isMatrix: one of Francesco's immunoStates data tables
#
# Outputs:
#   - myData: a data frame with two columns:
#       - cell type names
#       - Z-score of the gene list in that cell type
cellTypeEnrichment = function(myGenes, isMatrix){
  # Extract immunoStates ES for each iSEXS gene
  isexsCellTypes = subset(isMatrix, gene %in% myGenes)
  
  # Create a matrix where rows are ISDS genes, column are cell types
  esCellTypes = isexsCellTypes[,c("gene", "g", "cellType"), with =F]
  esCellTypes = dcast(data = esCellTypes, formula = gene ~ cellType, value.var = "g")
  rownames(esCellTypes) = esCellTypes$gene
  esCellTypes$gene = NULL
  esCellTypes = as.matrix(esCellTypes)
  
  # Put all above zero for geometic mean to work
  range(esCellTypes, na.rm = T)
  esCellTypes = esCellTypes + abs(min(esCellTypes, na.rm = T)) + 1.01
  
  # Get geometric mean for each cell type (column)
  cellEnrichment = apply(esCellTypes, 2, function(x) geomMean(x, na.rm = T))
  cellEnrichment = scale(cellEnrichment, center = T, scale = T)
  cellEnrichment = cellEnrichment[,1]
  cellEnrichment = sort(cellEnrichment)
  
  
  # Put into a data frame I can use for plotting
  myData = cbind(names(cellEnrichment), cellEnrichment)
  myData = as.data.frame(myData)
  colnames(myData) = c('celltype', 'value')
  
  # Fix columns to work well with plotting
  myData$value = as.numeric(as.character(myData$value))
  myData$celltype = factor(myData$celltype, levels = names(sort(cellEnrichment)))
  
  return(myData)
  
}


#.............................................#
##### Correlate Filter with Pheno Columns #####
#.............................................#
# getPhenoCorr
# correlates a score with a particular column in pheno
#
# Inputs: 
#   myDataset: Dataset object
#   myColumn: string, name of a $pheno column you want to correlate
#             (e.g. "sledai", or "cd4_Tcell_freq")
#   myScore: numeric vector, meant to be the output of calculateScore()
#
# Outputs: 
#   the correlation coefficient (r) and p-value of the correlation
#   I think it's Pearson's correlation
getPhenoCorr = function(myDataset, myColumn, myScore){
  corResults = cor.test(myDataset$pheno[,myColumn], myScore)
  return(list(r = signif(corResults$estimate,2), p=signif(corResults$p.value,2)))
}

# getPhenoCorrMulti
# Correlates many columns at once and aggregates them into a matrix
#
# Inputs: 
#   myDataset: dataset object
#   myColumns: vector of column names in $pheno to be correlated with the score
#              e.g. c('age', 'sledai', 'cd4_t_cells', 'b_cells')
#   myFilter: Filter object, describing the score to be calculated
getPhenoCorrMulti = function(myDataset, myColumns, myFilter){
  # Calculate score for this dataset
  myScore = calculateScore(myFilter, myDataset)
  
  # Calculate the correlations for each column
  allMyCorr = sapply(myColumns, function(x) getPhenoCorr(myDataset = myDataset, myColumn = x, myScore = myScore))
  
  # Reformat to be a matrix that we can sort easily in View
  allMyCorr = cbind(unlist(allMyCorr[1,]), unlist(allMyCorr[2,]))
  colnames(allMyCorr) = c('r', 'p')
  
  # Add in a column for FDR
  allMyCorr = as.data.frame(allMyCorr)
  allMyCorr$fdr =  p.adjust(allMyCorr$p, method = 'fdr')
  
  return(allMyCorr)
}

#.............#
##### SAM #####
#.............#
#~ SAM_GEM takes a GEM object as input, and performs SAM analysis on it
#~ Look up the vignette of the function SAM for details
#~
#~ Inputs: 
#~ GEM: usual GEM object with
#~     $expr: gene expression matrix. Rownames are probeIDs. Colnames are sampleIDs
#~     $class: vector of 1's for cases, and 0's for controls. Corresponds to the samples in the orders they are in $expr
#~     $keys: character vector of gene symbols. The names of $keys are the probeIDs listed in the rownames of $expr
#~ fdrThreshold: the maximum False Discovery Rate that you care about
#~ 
#~ Outputs: Look at the vignette of the function SAM for details
#~    samr.object: the object returned by the function samr
#~    siggenes.table: a list containing a table of up genes, table of down genes, number of up genes, number of downgenes
#~            siggenes.table is the most interesting one!
#~ delta.table
#~ del: value of delta (distance from 45 degree line in SAM plot)
#~ call: calling sequence
SAM_GEM <- function(gem, fdrThreshold, nperms = 100, rseed = 0) {
  require(samr)
  
  results <- SAM(x= gem$expr, y = (gem$class + 1),resp.type= "Two class unpaired", genenames= as.character(gem$keys), 
                 nperms = nperms, logged2 = TRUE, fdr.output = fdrThreshold, random.seed = rseed)
  
  return(results)
}

# subsetToFilterGenes
# Takes a dataset object and removes all probes that do not correspond
# to the genes in a filter object
#
# Takes probes that map to multiple genes into account
#
# Inputs: 
#   myDataset: Dataset object
#   myFilter: Meta Filter object
#
# Outputs: 
#   myDataset: Same as input, but probes that do not correspond to
#              filter genes have been removed
subsetToFilterGenes = function(myDataset, myFilter){
  # Check that they're the right data types
  if(!checkDataObject(myDataset, "Dataset")){
    error("Incorrect formatting in Dataset object")
  }
  
  if(!checkDataObject(myFilter, "MetaFilter")){
    error("Incorrect formatting in Filter object")
  }
  
  # Extract the genes in the filter
  filterGenes = c(myFilter$posGeneNames, myFilter$negGeneNames)
  
  # Get the probes that correspond to filterGenes
  filterProbes = sapply(myDataset$keys, function(x) length(intersect(filterGenes, unlist(strsplit(x = x,split=',' )))) > 0)
  
  # Subset expr and keys to just those probes
  myDataset$expr = myDataset$expr[filterProbes,]
  myDataset$keys = myDataset$keys[filterProbes]
  
  print(checkDataObject(myDataset, "Dataset"))
  
  return(myDataset)
}

