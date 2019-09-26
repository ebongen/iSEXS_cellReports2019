# quickPlots: A toolkit of quick and easy graphs
# Author: Erika Bongen
# Date: May 2017

# Purpose: 
# Quickly make graphs for MetaIntegrator Dataset and Meta objects
# Contains code for: 
#    - Gene level violin plots
#    - Gene level heatmaps
#    - Timecourse of gene expression
#    - Scatterplot comparing two genes' expression patterns


#....................#
##### Load Stuff #####
#....................#

require(ggplot2)
require(stringr)
require(RColorBrewer)
require(gplots)
require(MetaIntegrator)
theme_set(theme_gray(base_size = 18))

femColor = "#F42C04"
malColor = "#190E4F"
#..................#
##### strsplit #####
#..................#
#......................#
# strsplit_emptyString #
#......................#
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

#................#
# strsplitVector #
#................#
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

#................................................#
##### Old Way - extractMultiGenesFromDataset #####
#................................................#
# extractGeneFromDataset takes a gene symbol (e.g. "GAPDH")
# and a Dataset object and returns a vector of gene expression 
#
# Where multiple probes mapped to the same gene, the arithmatic mean was taken
#
# Probes that map to multiple genes were ignored
#
# Inputs: 
# myDataset: Dataset object 
# myGene: string, the gene symbol of your gene of interest
# Returns: 
# geneExpr: a vector where values are expression values of that gene
# and names are the sample name it corresponds to
# extractGeneFromDataset = function(myDataset, myGene) {
#   # Establish the pieces of myDataset
#   expr = myDataset$expr
#   keys = myDataset$keys
#   
#   # Grab all probes that correspond to myGene
#   indices = which(keys == myGene)
#   geneExpr = expr[indices,]
#   
#   # If there's more than one probe for this gene, take the average
#   # expression of each probe
#   if(!is.null(dim(geneExpr))) {
#     geneExpr = colMeans(geneExpr)
#   }
#   
#   return(geneExpr)
# }
# 
# # extractMultiGenesFromDataset is a wrapper for extractGeneFromDataset
# # where instead of getting you the expression values of one gene it
# # gives you a matrix where cols are samples, rows are genes, and each
# # row tells you the expression of a different gene in each sample. 
# #
# # Where multiple probes mapped to the same gene, the arithmatic mean was taken
# #
# # Probes that map to multiple genes were ignored
# #
# # Inputs: 
# # myDataset: Dataset object 
# # myGene: string, the gene symbol of your gene of interest
# # Returns: 
# # geneExpr: a matrix where rows are genes and cols are samples
# # each row is the expr of a different gene in each sample
# extractMultiGenesFromDataset = function(myDataset, myGenes){
#   geneExpr = sapply(myGenes, function(x) extractGeneFromDataset(myDataset, x))
#   return(t(geneExpr))
# }

#......................................#
##### New Way - getSampleLevelData #####
#......................................#
#....................#
# extractDataFromGEM #
#....................#
# Extracts expression of genes from a dataset object (gem)
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
# getSampleLevelGeneData #
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

#.......................#
# gem_expand.df #
#.......................#
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
#..........................#
##### Gene-Level Plots #####
#..........................#
#.............#
# quickViolin #
#.............#
quickViolin = function(myDataset, myGene, columnName, mySize = 3, fontSize = 12, myPalette=NULL){
  # Create myData, a data frame with phenotype and expression data for your gene
  geneExpr = getSampleLevelGeneData(myDataset, myGene)
  myData = myDataset$pheno
  myData$gene = as.vector(as.matrix(geneExpr))
  myData$trait = myData[[columnName]]
  
  # If there's no palette listed, use the default colors
  if(is.null(myPalette)) {
    
    if(columnName == "sex"){
      # if the trait column is sex, make it pink and blue
      myPalette = c(femColor, malColor)
      } else if( length(unique(myData$trait)) == 2){
        myPalette = c(femColor, malColor)
        } else if (length(unique(myData$trait)) <= 8){
          # If the trait column has 8 or less categories, use the Dark2 pallet
          myPalette = brewer.pal(length(unique(myData$trait)), "Dark2")
        } else {
          myPalette = rep("black", length(unique(myData$trait)))
        }
    }
  
  # Create the violin plot
  ggplot() +
    theme(text = element_text(size = fontSize)) +
    geom_violin(data=myData, aes(x=trait, y=gene),fill='grey',trim=F) +
    geom_jitter(data=myData, aes(x=trait, y=gene, col=trait), size = mySize, shape=19, position=position_jitter(width=0.1)) +
    scale_colour_manual(values=myPalette) +
    xlab(columnName) +
    ggtitle(myDataset$formattedName) +
    ylab(myGene)
}



#..............#
# quickHeatmap #
#..............#
# quickHeatmap makes a heatmap of a vector of genes from a particular dataset
#
# It defaults to an orange/purple == positive/negative color scheme, but you can 
# change it to "rainbow" if you want more different colors
#
# Inputs: 
#   myGenes: a character vector of gene symbols of your genes of interest
#   gem: a Khatri lab Dataset object, must contain expr and class
#
# Outputs: 
#   a heatmap of the expression of myGenes in the dataset(gem)
quickHeatmap <- function(myGenes, gem, columnName, colorPalette = "orangePurple", symbreaks=F, scale=F){
  # Extract expression matrix of genes of intrest(myGenes) from dataset(gem)
  myExpr = as.matrix(getSampleLevelGeneData(gem, myGenes))
  myExpr = na.omit(myExpr)
  
  if(scale){
    myExpr = t(scale(t(myExpr)))
  }
  
  # establish color pallette
  if(colorPalette == "rainbow"){
    # try new pallette 
    my_palette <- colorRampPalette(c("purple4", "blue", "dodgerblue", "white", "yellow", "darkorange1", "red"))
  } else if(colorPalette == "oneColor"){
    my_palette <- colorRampPalette(c("white", "white","purple","purple3","purple4"))(n = 500)
    } else {
      my_palette <- colorRampPalette(c("purple4","purple3","purple", "white", "orange","orange3","orange4"))(n = 500)
      }
  
  # Set the colors of the classes in columnName
  if(columnName == "sex"){
    # if the trait column is sex, make it pink and blue
    classColors = c(femColor, malColor)
    } else if( length(unique(gem$pheno[,columnName])) == 2){
    classColors = c("salmon", "dodgerblue2")
    } else if (length(unique(gem$pheno[,columnName])) <= 8){
    # If the trait column has 8 or less categories, use the Dark2 pallet
    classColors = brewer.pal(length(unique(gem$pheno[,columnName])), "Dark2")
    } 
  
  # Create a vector of colors that each sample should be marked as
  sampleClass = c()
  for (i in 1:length(gem$pheno[,columnName])){
    sampleClass[i] = classColors[which(unique(as.character(gem$pheno[,columnName])) == as.character(gem$pheno[i,columnName]))]
  }
  
  
  #Symbreaks is important to center around zero
  heatmap.2(myExpr,  col=my_palette, symbreaks=symbreaks, 
            ColSideColors = sampleClass,
            trace = "none",
            main = gem$formattedName)
  
  # Add legend
  par(lend = 1)
  legend("topright",      # location of the legend on the heatmap plot
         legend = unique(gem$pheno[,columnName]), # category labels
         col = classColors,  # color key
         lty= 1,             # line style
         lwd = 10            # line width
  )
}


#...........................#
### ## quickTimecourse ## ###
#...........................#
# Inputs: 
# myDataset: Dataset object
# myGene: string, Gene Symbol, (e.g. "GAPDH")
# timeColumn: the name of the column in $pheno that has numeric time information, (e.g. "time")
# timeUnits: string, e.g. ("hours", "days")
# main: the title of the plot
quickTimecourse <- function(myDataset, myGene, timeColumn, groupColumn, timeUnits = "", main=myDataset$formattedName) {
  myData = myDataset$pheno
  myData$geneExpr = as.vector(as.matrix(getSampleLevelGeneData(myDataset, myGene)))
  myData[["time"]] = myData[[timeColumn]]
  myData[["group"]] = myData[[groupColumn]]
  
  # Get units prepped
  if(timeUnits != ""){
    timeUnits = paste("(", timeUnits, ")", sep ="")
  }
  
  a = ggplot(data=myData, aes(x=time, y=geneExpr, colour=as.factor(group))) + geom_point(size=3) + geom_smooth(size=2) + ggtitle(main) +
    xlab(paste("Time", timeUnits, sep = " ")) +ylab(myGene)
  
  if(groupColumn == "sex"){
    a = ggplot(data=myData, aes(x=time, y=geneExpr, colour=as.factor(group))) + geom_point(size=3) + geom_smooth(size=2) + ggtitle(main) +
      xlab(paste("Time", timeUnits, sep = " ")) +ylab(myGene) + scale_color_manual(values = c(femColor, malColor))
  }
  
  return(a)
}

#................#
## quickScatter ##
#................#
# quickScatter creates a scatter plot of gene expression where the x and y 
# axes can either be the expression of a single gene, or the geometric
# mean of a vector of genes. 
#
# I built this function specifically to look at separateion of samples
# according to X-genes (e.g. XIST) vs Y-genes (e.g. RPS4Y1). 
# If this is what you're plotting, I recommend having female genes be "geneY"
# This just helps me think of it in terms of the XIST violin plots
#
# Inputs: 
#   myDataset: a Khatri lab Dataset object
#   geneX: Gene that goes on the x-axis
#          either a single gene symbol, or a vector of gene symbols
#          e.g. "XIST" or c("XIST", "RORA", "ZFX")
#   geneY: Gene that goes on the Y axis
#          either a single gene symbol, or a vector of gene symbols
#          e.g. "XIST" or c("XIST", "RORA", "ZFX")
quickScatter = function(myDataset, geneX, geneY, columnName, geneXLab="", geneYLab =""){
  # If geneX is only one gene, get raw values
  if(length(geneX) == 1){
    x = as.vector(as.matrix(getSampleLevelGeneData(myDataset, geneX)))

    # If geneXLab is empty, replace it with geneX's gene symobl
    if(geneXLab == ""){
      geneXLab = geneX
    }
    
  } else{
    x = getSampleLevelGeneData(myDataset, geneX) # extract expression of genes in geneX
    x = apply(x, 2, .geomMean) # Take the geometric mean of the gene's expression
    
    # Get x label that combines all gene names
    if (geneXLab == ""){
      geneXLab = paste(geneX, collapse = ',')
    }
  }
  
  
  # If geneY is only one gene, get raw values
  if (length(geneY) ==1){
    y = as.vector(as.matrix(getSampleLevelGeneData(myDataset, geneY)))
    
    # If geneYLab is empty, replace it with geneY's gene symbol
    if(geneYLab == ""){
      geneYLab = geneY
    }
    
  } else{
    y = getSampleLevelGeneData(myDataset, geneY) # extract expression of genes in geneX
    y = apply(y, 2, .geomMean) # Take the geometric mean of the gene's expression
    
    # Get y label that combines all gene names
    if(geneYLab == ""){
      geneYLab = paste(geneY, collapse = ',')
    }
  }
  
  # Create data frame from which to create plots
  myData = myDataset$pheno
  myData$x = x
  myData$y = y
  myData$group = myData[[columnName]]
  
  
  # Scatterplot
  ggplot(myData, aes(x, y, colour = group)) + geom_point(size = 3) +
    ylab(geneYLab) + xlab(geneXLab) + ggtitle(myDataset$formattedName)
  
}

#quickPCAplot
# Creates a scatter plot of the PC's of your choice given a dataset
# and a prcomp object
#
# Inputs: 
#   myDataset: Khatri lab Dataset object that PCA has been run on
#   myPCA: prcomp object
#   colName: string, column name in pheno to color dots by
#   PCx: int, the PC number to plot on the x-axis
#   PCy: int, the PC number to plot on the y-axis
#
# Generates: 
#   a colored scatterplot of the values of PCx and PCy
quickPCAplot  = function(myDataset, myPCA, colName, PCx=1, PCy=2){
  myDataset$pheno$PCx = myPCA$x[,PCx]
  myDataset$pheno$PCy = myPCA$x[,PCy]
  myDataset$pheno$myGroup = myDataset$pheno[[colName]]
  
  if ("sex" %in% colnames(myDataset$pheno)){
    ggplot(data = myDataset$pheno) +geom_point(aes(x=PCx, y=PCy, col = myGroup, shape = sex), size = 3) + 
      ggtitle(myDataset$formattedName) + xlab(paste("PC", PCx, sep= "")) + ylab(paste("PC", PCy, sep=""))
  } else{
    ggplot(data = myDataset$pheno) +geom_point(aes(x=PCx, y=PCy, col = myGroup), size = 3) + 
      ggtitle(myDataset$formattedName) + xlab(paste("PC", PCx, sep= "")) + ylab(paste("PC", PCy, sep=""))
  }
  
}

# quickScatterPheno
# Creates a scatter plot of a gene and a numerical column in pheno
# e.g. you wanna correlate expression of a gene with age
#
# Inputs: 
#  - myDataset: Khatri lab Dataset object
#  - myGene: String, gene symbol, plotted on y-axis
#  - myCol: string, name of column in $Pheno you wanna on x-axis
#  - colorColumn: string, the name of column you want to color points by
#
# Generates: 
#  - scatter plot of pheno column vs gene expression
quickScatterPheno = function(myDataset, myGene, myCol, colorColumn){
  myData = myDataset$pheno
  myData$gene = unlist(getSampleLevelGeneData(myDataset, myGene))
  myData[["xColumn"]] = unlist(myData[[myCol]])
  myData[["colorColumn"]] = unlist(myData[[colorColumn]])
  
  ggplot(myData, aes(x=xColumn, y=gene, col = colorColumn))+
    geom_point(size=3) + theme_cowplot() + xlab(myCol) + ylab(myGene)
  
}

# quickForestPlot
# Adaptation of ggForestPlot from MetaIntegrator that makes it have
# a simpler cowplot style
# 
# It can make a grid of plots, but the study names won't be pretty
#
# Inputs: 
#   - metaObject: a meta object
#   - genes: a vector of gene symbols
#   

quickForestPlot = function (metaObject, genes, confLevel = 0.95, facetCols = NULL, 
                            facetScales = "free_x", boxScales = c(6, 16), fontSize = 16) {
  # If the number of columns isn't given, automatically set it to 4
  if (is.null(facetCols)) {
    facetCols = ifelse(length(genes) == 1, 1, 4)
  }
  
  # Set the value of confidence intervals
  ci.value <- -stats::qnorm((1 - confLevel)/2)
  
  # Extract formation for the gene
  studyDataMeans <- as.data.table(metaObject$metaAnalysis$datasetEffectSizes, 
                                  keep.rownames = T)[rn %in% genes]
  studyDataSEs <- as.data.table(metaObject$metaAnalysis$datasetEffectSizeStandardErrors, 
                                keep.rownames = T)[rn %in% genes]
  pooledMean <- metaObject$metaAnalysis$pooledResults[genes, 
                                                      "effectSize"]
  pooledSE <- metaObject$metaAnalysis$pooledResults[genes, 
                                                    "effectSizeStandardError"]
  
  # Add useful info into studyDataMeans so it can be plotted
  studyDataMeans$Summary = pooledMean
  studyDataMeans = studyDataMeans %>% data.table::melt(id.vars = "rn") %>% 
    data.table::setnames("value", "means")
  studyDataSEs$Summary = pooledSE
  studyDataSEs = studyDataSEs %>% data.table::melt(id.vars = "rn") %>% 
    data.table::setnames("value", "SEs")
  studyData = merge(studyDataMeans, studyDataSEs, by = c("rn", 
                                                         "variable")) %>% data.table::setnames(c("rn", "variable"), 
                                                                                               c("gene", "study"))
  # Get information for the summary diamond
  pooled = data.table(gene = genes, means = pooledMean, SEs = pooledSE)
  pooled$study = "Summary"
  pooled[, `:=`(cilb, means - ci.value * SEs)]
  pooled[, `:=`(ciub, means + ci.value * SEs)]
  pooled[, `:=`(RelConf, 1/(ciub - cilb))]
  studyData[, `:=`(ok, is.finite(means + SEs))]
  if (is.null(xlim)) {
    studyData[, `:=`(xlim.lb, min(means[ok] - ci.value * 
                                    SEs[ok], na.rm = TRUE))]
    studyData[, `:=`(xlim.ub, max(means[ok] + ci.value * 
                                    SEs[ok], na.rm = TRUE))]
  }
  studyData[, `:=`(cilb, means - ci.value * SEs)]
  studyData[, `:=`(ciub, means + ci.value * SEs)]
  studyData[, `:=`(RelConf, 1/((ciub - cilb)))]
  studyData$study = stats::relevel(studyData$study, ref = "Summary")
  studyData[study == "Summary", `:=`(means = NA, SEs = NA, 
                                     ok = FALSE, cilb = NA, ciub = NA, RelConf = NA)]
  
  # Get variables necessary for plotting
  nStudies = length(unique(studyData$study))
  labelstext = c("Summary", as.character(unique(studyData$study)))
  labelspos = seq_along(labelstext) - 1
  pooled_poly = cbind(pooled[, list(gene, cilb, means, ciub, 
                                    means)] %>% data.table::melt(id.vars = "gene") %>% data.table::setnames("value", 
                                                                                                            "x"), pooled[, list(gene, 1.02, 0.85, 1.02, 1.15)] %>% 
                        data.table::melt(id.vars = "gene") %>% data.table::setnames("value", 
                                                                                    "y"))[, .(gene, x, y)]
  # Get pretty names, so long as they match!
  if(all(names(metaObject$originalData)== studyData$study[1:length(studyData$study)-1])){
    prettyNames = sapply(metaObject$originalData, function(x) x$formattedName)
    studyData$study = factor(c(prettyNames, 'Summary'), levels = c('Summary', prettyNames))
  }
  
  # 
  ggplot(studyData, aes_string("means", "study")) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_errorbarh(color = "black", height = 0.1, 
                   aes_string(xmax = "ciub", xmin = "cilb")) + 
    geom_point(aes_string(size = "RelConf"), shape = 15, color = "black") + 
    scale_size(range = boxScales, guide = FALSE) + 
    theme_cowplot(font_size = fontSize)+
    theme(strip.background =element_rect(fill="white"), axis.title.y = element_blank()) + 
    xlab("Standardized Mean Difference (log2 scale)") + 
    facet_wrap(~gene, ncol = facetCols, scales = facetScales) + 
    geom_polygon(data = pooled_poly, aes_string(x = "x", 
                                                y = "y"), fill = "grey33") + 
    coord_cartesian(ylim = c(1.2, nStudies))
}
