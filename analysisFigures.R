# Create figures for iSEXS paper

# Purpose: 
#   - Sex and gender play noticable roles in disease susceptibility
#   - Can we identify transcriptional differences between healthy
#     female and male immune systems?


# Note: 
#   - My original iSEXS thresholds lead to ~20% garbage genes in validation
#   - Garbage == -0.1 < ES < 0.1
#   - So, we need to tweak the threshold to lead to fewer FP, without losing too many TP


# Outline: 
#  Fig 1: Schematic (to be made in PowerPoint or bullshit Photoshop)
#  Fig 2: Heatmap introducing iSEXS genes with chromosomal location annotated
#  Fig 3: Scores work in XX-female, XY-males, and XXY-males
#  Fig 4: immunoStates implicates CD4's and monocytes
#  Fig 5: sex differences due to immune cell dynamics during influenza infection
#  Fig 6: Auto-iSEXS and Ab response
#....................#
##### Load Stuff #####
#....................#
# Set working directory
myDir = "/labs/khatrilab/ebongen/sexDifferences2.0/isexs/"
setwd(myDir)

# Load libraries
library(rmeta)
library(MetaIntegrator)

library(ggplot2)
library(magick)
library(cowplot)
library(RColorBrewer)
library(colorspace)
library(ggsignif)
library(pheatmap)
library(ggbeeswarm)
library(car) # For ANOVA


# Source code
source("general_GEMfunx.R")
source("~/Tools/Graphing Scripts/quickPlots.R")

### Initialize constants
# The index number for a figure, so that they'll be listed in /plots in the order I make them
# Allows me to insert new figures without having to change plot names manually
# Does not correspond directly to "Figure 1" of the actual paper
figIndex = 0 
nReps = 800

# Load items 
load("sexMetaObj.RData")


#..........................#
##### Validation Stats #####
#..........................#
# Purpose:
#   - Get p-values for iSEXS genes in validation cohorts

# Female-assoc genes
valFem = valMetaObj$metaAnalysis$pooledResults[sexMetaObj$filterResults$isexs$posGeneNames,]
valFem = subset(valFem, effectSize > 0) # get only positive genes
nrow(valFem) # 85/94 have positive effect sizes
sum(valFem$effectSizePval < 0.05) # 47

# Male-assoc genes
valMal = valMetaObj$metaAnalysis$pooledResults[sexMetaObj$filterResults$isexs$negGeneNames,]
valMal = subset(valMal, effectSize < 0)
nrow(valMal) # 45
sum(valMal$effectSizePval < 0.05)


# Genes that show the same trend
nrow(valFem) + nrow(valMal) # 130

# Genes that are significant p-values: 80
sum(valFem$effectSizePval < 0.05) + sum(valMal$effectSizePval < 0.05)


#................................#
##### Heatmap with Chr Loc #####
#................................#
# Purpose: 
#   - Create a heatmap with:
#         - A) Discovery effect sizes sorted by summary effect size
#         - B) Validation effect sizes sorted by Discovery summary effect size
#         - Rows are datasets, columns are genes
#         - Columns are color coded by chromosome location


### Make country into the formatted name
# Discovery
for(i in 1:length(sexMetaObj$originalData)){
  sexMetaObj$originalData[[i]]$formattedName = sexMetaObj$originalData[[i]]$country
}

# Validation
for(i in 1:length(valMetaObj$originalData)){
  valMetaObj$originalData[[i]]$formattedName = valMetaObj$originalData[[i]]$country
}


### Chosen colors:
#   1) Autosomal is grey "#CDCFD0"
#   2) Y is dark blue, malColor "#190E4F"
#   3) PAR1 is light blue, "#7180AC"
#   4) X (non escape) is light red, "#FAA08E"
#   5) X escape is red, femColor "#F42C04"


# Create simplified location
ISEXSLoc_simple = ISEXSLoc
ISEXSLoc_simple[grepl(x=ISEXSLoc_simple, pattern = "[0-9]") & ISEXSLoc_simple != 'PAR1'] = "Autosomal" 
ISEXSLoc_simple[ISEXSLoc_simple == 'X_escape'] = 'Known X Escape Gene'

# Assign each location to a color
colorDict = list(Chromosome=c('PAR1' ="#7180AC", 'Y' = malColor, 'X' =  "#FAA08E", 'Known X Escape Gene'= femColor, 'Autosomal' ="#CDCFD0"))


### Get ES from Discovery and Validation cohorts
discES = metaHeatmatrix(metaObject = sexMetaObj, filterObject = sexMetaObj$filterResults$isexs, sortByTissue = T)
tissue_discovery = discES$tissue
discES = discES$heatmatrix
dim(discES) # 7 rows (6 Disc dataset + 1 Pooled); 144 columns 

valES = metaHeatmatrix(metaObject = sexMetaObj, filterObject = sexMetaObj$filterResults$isexs, valObject = valMetaObj, sortByTissue = T)
tissue_val = valES$tissue
valES = valES$heatmatrix
dim(valES) # 11 rows (10 val datasets + 1 Pooled); 144 columns

### Combine Discovery and Validation ES information
all(colnames(discES) == colnames(valES)) # The genes match, so ready to combine
comboES = rbind(discES, valES)
rownames(comboES)[1] = "Pooled Discovery"
rownames(comboES)[nrow(discES)+1] = 'Pooled Validation'

### Create tissue color code
colorDict$tissue = c('Pooled Validation' = 'white', 'Pooled Discovery' = 'white','PBMC' = 'black', 'Whole Blood' = 'white')
tissueLabels = c('Pooled Discovery', sort(tissue_discovery), 'Pooled Validation', sort(tissue_val))
tissueLabels = data.frame(tissue = factor(tissueLabels))
rownames(tissueLabels) = rownames(comboES)

### Set colors
paletteLength = 100
myBreaks <- c(seq(min(comboES, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(comboES, na.rm = T)/paletteLength, max(comboES, na.rm = T), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(c("purple", "white", "orange"))(paletteLength)
#myColor <- colorRampPalette(c(malColor, "white", femColor))(paletteLength)


### Make heatmap
fig2_heatmap = pheatmap(comboES,
                        annotation_col = data.frame(Chromosome = factor(ISEXSLoc_simple)), 
                        annotation_row = tissueLabels,
                        annotation_colors = colorDict,
                        annotation_names_col = F,
                        annotation_names_row = F,
                        annotation_legend = F,
                        cluster_cols = F, cluster_rows = F,
                        border_color = NA, treeheight_row = 0,
                        breaks = myBreaks, 
                        color = myColor,
                        cellwidth = 3, cellheight = 15,
                        fontsize = 14, show_rownames = T,
                        show_colnames = F, legend=T, 
                        gaps_row = c(7,7,7))

### Add annotation
fig2_heatmap = plot_grid(fig2_heatmap[[4]], nrow = 1, ncol=1)

# Variables for annotation location
colorBox_yMin = 0.945
colorBox_yMax = colorBox_yMin + 0.03

colorBox_xMin = 0.075
colorBox_xMax = colorBox_xMin + 0.015

fig2_heatmap_annot = fig2_heatmap + 
  geom_rect(aes(xmin=0.04, xmax=0.065, ymin=0, ymax=1), colour='white', fill='white')+
  draw_label('Y Chr.        PAR1       Autosome      X Chr.', x=0.255, y=0.96)+
  geom_rect(aes(xmin=colorBox_xMin, xmax=colorBox_xMax, ymin=colorBox_yMin, ymax=colorBox_yMax), colour=colorDict$Chromosome['Y'], fill=colorDict$Chromosome['Y']) +
  geom_rect(aes(xmin=colorBox_xMin + 0.087, xmax=colorBox_xMax + 0.087, ymin=colorBox_yMin, ymax=colorBox_yMax), colour=colorDict$Chromosome['PAR1'], fill=colorDict$Chromosome['PAR1']) +
  geom_rect(aes(xmin=colorBox_xMin + 0.165, xmax=colorBox_xMax + 0.165, ymin=colorBox_yMin, ymax=colorBox_yMax), colour=colorDict$Chromosome['Autosomal'], fill=colorDict$Chromosome['Autosomal'])+
  geom_rect(aes(xmin=colorBox_xMin + 0.271, xmax=colorBox_xMax + 0.271, ymin=colorBox_yMin, ymax=colorBox_yMax), colour=colorDict$Chromosome['X'], fill=colorDict$Chromosome['X'])+
  geom_rect(aes(xmin=colorBox_xMin+0.379, xmax=colorBox_xMax+0.379, ymin=colorBox_yMin, ymax=colorBox_yMax), colour=colorDict$Chromosome['Known X Escape Gene'], fill=colorDict$Chromosome['Known X Escape Gene']) +
  draw_label('X Escapee', x=0.515, y=0.96)+
  draw_label('PBMC', x=0.035, y=0.78) +
  draw_label(paste('Whole', '\n', 'Blood', sep=''), x=0.035, y=0.66) +
  draw_label('PBMC', x=0.035, y=0.41) +
  draw_label(paste('Whole', '\n', 'Blood', sep=''), x=0.035, y=0.19) +
  draw_label('a', fontface = 'bold', x=0.01, y= 0.945, size = 18)+
  draw_label('b', fontface = 'bold', x=0.01, y=0.525, size= 18)


save_plot(filename = 'plots/fig2_comboHeatmap.pdf',fig2_heatmap_annot, base_width = 11.25, base_height = 5)


#.....................................#
##### Milieu Interieur Validation #####
#.....................................#
# Purpose: 
#   - Milieu Interieur cohort measured 13/144 iSEXS genes via Nanostring
#   - See if the effect sizes of young people with iSEXS genes are comparable

### Chosen colors:
#   1) Autosomal is grey "#CDCFD0"
#   2) Y is dark blue, malColor "#190E4F"
#   3) PAR1 is light blue, "#7180AC"
#   4) X (non escape) is light red, "#FAA08E"
#   5) X escape is red, femColor "#F42C04"

# Load Milieu cohort
load("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/4_sortedCells/cellProp/patin/milieuMeta.RData")

# Get isexs genes that are shared
milieu_isexs = isexs[isexs %in% milieuMeta$originalData$milieuGEM_young$keys]
length(milieu_isexs) # 13 genes are shared
length(isexs) # 144 total
table(ISEXSLoc_simple[milieu_isexs]) # 10 Autosomal, 1 X-escape, 1 X, 1 PAR1

# Create matrix of effect sizes
milieuES = cbind(milieuMeta$metaAnalysis$datasetEffectSizes[milieu_isexs,"milieuGEM_young"],
                     sexMetaObj$metaAnalysis$pooledResults[milieu_isexs,'effectSize'],
                     valMetaObj$metaAnalysis$pooledResults[milieu_isexs,'effectSize'])

colnames(milieuES) = c('milieu', 'discovery', 'validation')
milieuES = as.data.frame(milieuES)
milieuES$chrLoc = ISEXSLoc_simple[rownames(milieuES)]

# Add gene symbols for the points we want to label
milieuES$geneSymbol = rownames(milieuES)
milieuES$geneSymbolMal = rownames(milieuES)
milieuES$geneSymbolFem = rownames(milieuES)

# Make a column for the male-assoc. genes we want to highlight
myGenesMal = c('CD99', 'CTSG', 'CYBB', 'PSMB7')
milieuES$geneSymbolMal[!milieuES$geneSymbolMal %in% myGenesMal] = NA 

# Make a column for the fem-assoc. genes we want to highlight
myGenesFem = c('CD40LG', 'CD27') 
milieuES$geneSymbolFem[!milieuES$geneSymbolFem %in% myGenesFem] = NA 
milieuES$geneSymbolIL23R = rownames(milieuES)
milieuES$geneSymbolIL23R[!milieuES$geneSymbolIL23R =='IL23R'] = NA

# Look at which genes ware where
ggplot(milieuES, aes(x=validation, y=milieu, col =chrLoc)) + 
  geom_hline(yintercept = 0, col='black', size=0.75, linetype='dashed') + 
  geom_vline(xintercept = 0, col='black', size=0.75, linetype='dashed')+
  geom_point(size=2)+
  geom_text(aes(label=geneSymbol, fontface='bold'),hjust=0, vjust=0) +
  theme_cowplot() + ggtitle('Milieu Cohort vs Discovery Effect Sizes') +
  xlab('Discovery Effect Sizes') + ylab('Milieu Interieur Effect Sizes') +
  scale_color_manual(values = c("grey48", femColor, "#7180AC",  "#FAA08E"),name = 'Chromosomal Location')+
  theme(legend.position = 'none')



# Create gene-specific horizontal shift
milieuES_hjust = rep(0, nrow(milieuES))
names(milieuES_hjust) = milieuES$geneSymbol
milieuES_hjust[c('CD40LG', 'IL23R', 'NT5E', 'PSMB7')] = 1 # Shift some to the left
milieuES_hjust['DPP4'] = 0.5


# Create gene-specific
milieuES_vjust = rep(0, nrow(milieuES))
names(milieuES_vjust) = milieuES$geneSymbol
milieuES_vjust[c('CD40LG', 'CD27', 'CEACAM6', 'PSMB7', 'NT5E')] =1
milieuES_vjust['IL23R'] = 0.3
# Set the font size
fontSize = 4

a = ggplot(milieuES, aes(x=discovery, y=milieu, col =chrLoc)) + 
  geom_hline(yintercept = 0, col='black', size=0.75, linetype='dashed') + 
  geom_vline(xintercept = 0, col='black', size=0.75, linetype='dashed')+
  geom_label(aes(label=geneSymbol, fontface='plain'),hjust=milieuES_hjust, 
             vjust= milieuES_vjust, size=fontSize) +
  geom_point(size=2)+ theme_cowplot() + 
  xlab('Discovery Effect Sizes') + ylab('Milieu Interieur Effect Sizes') +
  scale_color_manual(values = c("grey48", femColor, "#7180AC",  "#FAA08E"),name = 'Chromosomal Location')

# Readjust the vertical and horizontal shift
milieuES_hjust[c('CD27')] = 1 # 1 shifts it left
milieuES_hjust['CD96'] = 1
milieuES_hjust[c('IL23R')] =0.6 # 0 shifts it right
milieuES_hjust['CEACAM6'] = 0.7
milieuES_hjust['DPP4'] = 0.7
milieuES_vjust[c('DPP4')] = 1 # 1 shifts it down
milieuES_vjust['NT5E'] = 0.5
milieuES_vjust['IL23R']=0.1

b = ggplot(milieuES, aes(x=validation, y=milieu, col=chrLoc)) + 
  geom_hline(yintercept = 0, col='black', size=0.75, linetype='dashed') + 
  geom_vline(xintercept = 0, col='black', size=0.75, linetype='dashed')+
  geom_label(aes(label=geneSymbol, fontface='plain'),hjust=milieuES_hjust, 
             vjust= milieuES_vjust, size=fontSize) +
  geom_point(size=2)+ theme_cowplot()+ 
  xlab('Validation Effect Sizes') + ylab('Milieu Interieur Effect Sizes')+
  scale_color_manual(values = c("grey48", femColor, "#7180AC",  "#FAA08E"), name = 'Chromosomal Location')




# milieuPlot = plot_grid(a +theme(legend.position = 'none'),b, labels = c('a','b'), 
#                        nrow = 1, ncol = 2, rel_widths = c(1,1.6))
# save_plot(filename = 'plots/milieuPlot.pdf',milieuPlot,ncol = 1, nrow=1, base_width = 12, base_height = 5)

# Save plot of Discovery ES vs Milieu ES
save_plot('plots/fig2c_milieuScatter.pdf', a + theme_cowplot(font_size = 18), base_aspect_ratio = 1.8, base_height = 4.6)


### ugg, p-values
milieuExpr = getSampleLevelGeneData(milieuMeta$originalData$milieuGEM_young, milieu_isexs)
milieuSex = milieuMeta$originalData$milieuGEM_young$pheno$sex

milieuPvals = apply(milieuExpr, 1, function(x) t.test(x[milieuSex == 'Female'], x[milieuSex == 'Male'])$p.value)
sort(milieuPvals < 0.05)

#.............................#
##### Milieu Forest Plots #####
#.............................#
# Purpose: 
#   - Showcase some of the very consistent genes
#   - Pick from the Milieu cohort

# Make forest plots
forestPlots = list()
forestPlots$CD40LG = quickForestPlot(valMetaObj, 'CD40LG', fontSize = 24) 
forestPlots$CD99 = quickForestPlot(valMetaObj, 'CD99', fontSize = 24) 
forestPlots$CTSG = quickForestPlot(valMetaObj, 'CTSG', fontSize = 24) 

# Make a grid
forestGrid = plot_grid(plotlist = forestPlots, nrow=1, ncol=3)

save_plot('plots/fig2def_forestPlots.pdf',forestGrid, nrow=1, ncol =3, base_height = 8, base_width = 8)

#.............................................#
##### Discovery - XY- and Auto- ROC Plots #####
#.............................................#
# Purpose:
#   - ROC plots of Discovery cohorts

# Make a version of sexMetaObj where we can 
# make ROC plot specific annotation
sexMetaObj_sorted = sexMetaObj



# Add tissue formatted name
for(i in 1:length(sexMetaObj_sorted$originalData)){
  tissueType = unique(sexMetaObj_sorted$originalData[[i]]$pheno$tissue)
  if(tissueType == 'PBMC'){
    # Add tissue to formatted name
    sexMetaObj_sorted$originalData[[i]]$formattedName = paste('PBMC', '-', sexMetaObj_sorted$originalData[[i]]$formattedName)
    # change name within list so that it'll sort by p's and w's
    names(sexMetaObj_sorted$originalData)[i] = paste('p', names(sexMetaObj_sorted$originalData)[i], sep='')
  }else{
    # Create pretty tissue formatted name "WB - GSE0000"
    sexMetaObj_sorted$originalData[[i]]$formattedName = paste('WB', '-', sexMetaObj_sorted$originalData[[i]]$formattedName)
    # Add w in front of names e.g. originalData$wgse0000
    names(sexMetaObj_sorted$originalData)[i] = paste('w', names(sexMetaObj_sorted$originalData)[i], sep='')
    
  }
}

# Check new names
names(sexMetaObj_sorted$originalData) # have w's and p's to sort them

# Check formatted names
for(myDataset in sexMetaObj_sorted$originalData){
  print(myDataset$formattedName)
}


# Sort the order of the Dataset objects
sexMetaObj_sorted$originalData = sexMetaObj_sorted$originalData[sort(names(sexMetaObj_sorted$originalData))]
checkDataObject(sexMetaObj_sorted, 'Meta', 'Pre-Filter') #True


### XY Score
sFig1a_core = MetaIntegrator::summaryROCPlot(sexMetaObj_sorted, xySig, bootstrapReps = nReps)
sFig1 = list()

sFig1$a_XYsig = sFig1a_core + theme_cowplot() + 
  theme(legend.background = element_rect(fill=NA), 
        legend.position = c(0.15, 0.25),
        legend.text=element_text(size=10),
        plot.title = element_text(size=18)) +
  labs(color = '') + ggtitle('XY-iSEXS Score in Discovery Cohorts')+
  theme(plot.title = element_text(size = 15, face = "bold"))


### Autosomal Score
sFig1b_core = MetaIntegrator::summaryROCPlot(sexMetaObj_sorted, autoSig, bootstrapReps = nReps)

sFig1$b_autoSig = sFig1b_core + theme_cowplot() + 
  theme(legend.background = element_rect(fill=NA), 
        legend.position = c(0.10, 0.25),
        legend.text=element_text(size=10),
        plot.title = element_text(size=18)) + 
  labs(color='') + ggtitle('Autosomal-iSEXS Score in Discovery Cohorts') +
  theme(plot.title = element_text(size = 15, face = "bold"))


# Save plot
sFig1_grid = plot_grid(plotlist = sFig1, labels = LETTERS[1:2])
save_plot('plots/sFig1_discoveryROC.pdf', sFig1_grid, ncol=2, base_height = 5)



#...................................#
##### Validation -  XY- and Auto- ROC Plots #####
#...................................#
# Purpose: 
#   - ROC plots of Validation cohorts
#   - Show accuracy of XY-iSEXS and auto-iSEXS

# Create a version of valMetaObj where we can make the annotation
# specific to the ROC plot needs
valMetaObj_sorted = valMetaObj


# Add tissue formatted name
for(i in 1:length(valMetaObj_sorted$originalData)){
  tissueType = unique(valMetaObj_sorted$originalData[[i]]$pheno$tissue)
  if(tissueType == 'PBMC'){
    # Add tissue to formatted name
    valMetaObj_sorted$originalData[[i]]$formattedName = paste('PBMC', '-', valMetaObj_sorted$originalData[[i]]$formattedName)
    # change name within list so that it'll sort by p's and w's
    names(valMetaObj_sorted$originalData)[i] = paste('p', names(valMetaObj_sorted$originalData)[i], sep='')
  }else{
    # Create pretty tissue formatted name "WB - GSE0000"
    valMetaObj_sorted$originalData[[i]]$formattedName = paste('WB', '-', valMetaObj_sorted$originalData[[i]]$formattedName)
    # Add w in front of names e.g. originalData$wgse0000
    names(valMetaObj_sorted$originalData)[i] = paste('w', names(valMetaObj_sorted$originalData)[i], sep='')
    
  }
}

# Check object names
names(valMetaObj_sorted$originalData)

# Look at formatted Names
for(myDataset in valMetaObj_sorted$originalData){
  print(myDataset$formattedName)
}

# Sort the order of the Dataset objects
valMetaObj_sorted$originalData = valMetaObj_sorted$originalData[sort(names(valMetaObj_sorted$originalData))]
checkDataObject(valMetaObj_sorted, 'Meta', 'Pre-Filter') #True

### XY Score
fig3a_core = MetaIntegrator::summaryROCPlot(valMetaObj_sorted, xySig, bootstrapReps = nReps)

fig3a = fig3a_core + theme_cowplot() + 
  theme(legend.background = element_rect(fill=NA), 
        legend.position = c(0.30, 0.35),
        legend.text=element_text(size=10),
        plot.title = element_text(size=18)) +
  labs(color = '') + ggtitle('XY-iSEXS Score in Validation Datasets')


### Autosomal Score
fig3b_core = MetaIntegrator::summaryROCPlot(valMetaObj_sorted, autoSig, bootstrapReps = nReps)

fig3b = fig3b_core + theme_cowplot() + 
  theme(legend.background = element_rect(fill=NA), 
        legend.position = c(0.24, 0.35),
        legend.text=element_text(size=10),
        plot.title = element_text(size=18)) +
  labs(color = '') + ggtitle('Autosomal-iSEXS Score in Validation Datasets')


# Summary AUC .76 (.7-.82)


#..........................#
##### iSEXS Stat Table #####
#..........................#
# Purpose:
#   - Table with dataset effect sizes
#     - Discovery summary effect sizes
#     - Discovery p-val
#     - Validation summary effect sizes
#      - Validation p-val

# Get summary stats from discovery
discSumm = sexMetaObj_sorted$metaAnalysis$pooledResults[isexs, c('effectSize', 'effectSizePval')]
colnames(discSumm) = paste('discovery_',colnames(discSumm), sep='')

discES = sexMetaObj_sorted$metaAnalysis$datasetEffectSizes[isexs,]
for(i in 1:length(sexMetaObj_sorted$originalData)){
  colnames(discES)[i] = sexMetaObj_sorted$originalData[[i]]$formattedName
}

# Make sure they all match!
all(rownames(discSumm) == rownames(discES))
discStats = cbind(discES, discSumm)


# Get summary stats for validation
valSumm = valMetaObj_sorted$metaAnalysis$pooledResults[isexs, c('effectSize', 'effectSizePval')]
colnames(valSumm) = paste('validation_', colnames(valSumm), sep='')
sum(rownames(valSumm) == isexs) # only 143 match
rownames(valSumm) = isexs


# Get validation effect sizes
isexs[!isexs %in% rownames(valMetaObj_sorted$metaAnalysis$datasetEffectSizes)] # "ZC3H11B" missing

valES = valMetaObj_sorted$metaAnalysis$datasetEffectSizes[isexs[isexs %in% rownames(valMetaObj_sorted$metaAnalysis$datasetEffectSizes)],]
for(i in 1:length(valMetaObj_sorted$originalData)){
  colnames(valES)[i] = valMetaObj_sorted$originalData[[i]]$formattedName
}

# Combine valSumm and valES
nrow(valSumm) # full 144
nrow(valES) # only 143
valES = rbind(valES, rep(NA, ncol(valES)))
rownames(valES)[nrow(valES)] = "ZC3H11B"
valES = valES[isexs,]
all(rownames(valES) == rownames(valSumm)) # all match!
valStats = cbind(valES, valSumm)

# Combine Discovery and Validation
isexsStats = cbind(discStats, valStats)
isexsStats = cbind(ISEXSLoc_simple[isexs], isexsStats)
colnames(isexsStats)[1] = 'Chr Location'

write.csv(isexsStats,file = 'reference/isexsStats.csv', quote=F)


#......................................#
##### Klinefelter's Beeswarm Plots #####
#......................................#
# Purpose: 
#   - Create figure showing how autosomal and XY sig work
#     in Klinefelter's
#
# Outline:
#    4a) XY-iSEXS in Klinefelter's
#    4b) Auto-iSEXS in Klinefelter's
#    4c) Heatmap of Klinefetler's associated genes
#
# Color for XXY: 'grey48'

# Load Klinefelter's Datasets
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')
load("0_datasets/3_humanPerturbations/gse42331_klinefelter/gse42331_XXY.RData")
load("0_datasets/3_humanPerturbations/gse47584_klinefelter/gse47584_XXY.RData")
setwd(myDir)


### Calculate scores
myDataXXY = gse42331$pheno
myDataXXY$xyScore = calculateScore(xySig, gse42331)
myDataXXY$autoScore = calculateScore(autoSig, gse42331)


### XXY - XY-sig

# XXY vs XX-Females  p < 2.2E-16 --> ****
t.test(myDataXXY$xyScore[myDataXXY$group == 'XXY Male'],
       myDataXXY$xyScore[myDataXXY$group == 'XX Female'])
# XXY vs XY-Males   p = 0.002 --> **
t.test(myDataXXY$xyScore[myDataXXY$group == 'XXY Male'],
       myDataXXY$xyScore[myDataXXY$group == 'XY Male'])

### XXY beeswarm+boxplot
fig3c_XY_XXY = ggplot(myDataXXY, aes(x=group, y=xyScore)) +
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5, priority = 'random') + 
  scale_color_manual(values = c(femColor, 'grey48', malColor))+
  theme_cowplot(font_size = 18)+
  ggtitle("XY-iSEXS   GSE42331") + ylab("XY-iSEXS Score") + xlab('') +
  geom_signif(annotations = '****', y_position = 2, textsize = 7, 
              xmin =1.1, xmax=1.9)+
  geom_signif(xmin=2.1, xmax=2.8, 
              annotations = '**', y_position = 0, textsize = 7) 




### XXY - Auto-Sig
# XXY vs XX-Females  NS
t.test(myDataXXY$autoScore[myDataXXY$group == 'XXY Male'],
       myDataXXY$autoScore[myDataXXY$group == 'XX Female'])
# XXY vs XY-Males p = 0.002
t.test(myDataXXY$autoScore[myDataXXY$group == 'XXY Male'],
       myDataXXY$autoScore[myDataXXY$group == 'XY Male'])


fig3d_auto_XXY = ggplot(myDataXXY, aes(x=group, y=autoScore)) + 
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5, priority = 'random') + 
  scale_color_manual(values = c(femColor, 'grey48', malColor))+
  theme_cowplot(font_size = 18)+
  ggtitle("Autosomal-iSEXS   GSE42331") + ylab("Autosomal-iSEXS Score") + xlab('') +
  geom_signif(annotations = 'N.S.', y_position = 1.2, textsize = 5, 
              xmin =1.1, xmax=1.9)+
  geom_signif(xmin=2.1, xmax=2.8, 
              annotations = '**', y_position = 1.2, textsize = 7) 


### Quick check of XXY-males with X-chr genes higher expressed in females
# Question: Is male-like XY-sig driven by Y-chromosome genes?
# Answer: Yes, the fem-assoc. X-chromosome genes are expressed at fem-like levels in XXY-males
xSig = xySig
xSig$negGeneNames = c('POOPS')

myDataXXY$xScore = calculateScore(xSig, gse42331)


ggplot(myDataXXY, aes(x=group, y=xScore)) + geom_boxplot(fill='grey') + 
  theme_cowplot(font_size = 18)+
  ggtitle("X-Chr Higher in Females - GSE42331") + ylab("X Fem Genes") + xlab('') +
  geom_signif(annotations = 'N.S.', y_position = 1.5, textsize = 5, 
              xmin =1.1, xmax=1.9)+
  geom_signif(xmin=2.1, xmax=2.8, 
              annotations = '****', y_position = 1.5, textsize = 7) 

# XXY vs XX-Females  NS
t.test(myDataXXY$xScore[myDataXXY$group == 'XXY Male'],
       myDataXXY$xScore[myDataXXY$group == 'XX Female'])
# XXY vs XY-Males p = 2E-6
t.test(myDataXXY$xScore[myDataXXY$group == 'XXY Male'],
       myDataXXY$xScore[myDataXXY$group == 'XY Male'])

pToStars(0.0000027)

#...........................#
##### Fig 3 - ASSEMBLE! #####
#...........................#
# Purpose: 
#  - Assemble subparts of Figure 3
#  - Will show that autosomal-iSEXS score also works
#  - Showcases particularly consistent X genes, and autosomal genes

# Assemble rows
fig3 = plot_grid(fig3a, fig3b, 
                 fig3c_XY_XXY + theme(legend.position = 'none'), 
                 fig3d_auto_XXY+theme(legend.position = 'none'), 
                 nrow = 2, ncol = 2, labels = LETTERS[c(1:4)], label_size = 18)


# Save plot!
save_plot('plots/fig3_autosomeGenes.pdf', plot = fig3, nrow = 2, ncol = 2, 
          base_height = 5.2)


#.....................................#
##### Klinefelter's Meta-Analysis #####
#.....................................#
# Purpose: 
#  - run a meta-analysis across two Klinefelter's datasets
#  - Identify iSEXS genes that are differentially expressed between
#    XXY-males and XY-males

### Prep gse42331
phenoXXY = subset(gse42331$pheno, sex == "male")
gse42331_xy = subsetGEMFromPheno(gse42331, phenoXXY)

# Create class vector
gse42331_xy$class = createClassVector(groupColumn = gse42331_xy$pheno$group, casesAre = "XXY Male", pheno = gse42331_xy$pheno)

# This Dataset does not have probes that map to multiple genes, so I don't have to expand it
any(sapply(gse42331_xy$keys,function(x) grepl(',', x)))

# Make it only have iSEXS genes
keys = gse42331_xy$keys[which(gse42331_xy$keys %in% isexs)]
expr = gse42331_xy$expr[names(keys), ]
gse42331_xy$expr = expr
gse42331_xy$keys = keys
checkDataObject(gse42331_xy, "Dataset") # True!

### Prep gse47584
# It's already XY (0's) vs XXY (1's)
table(gse47584$class, gse47584$pheno$group)

# Limit to ISEXS genes
keys = gse47584$keys[which(gse47584$keys %in% isexs)]
expr = gse47584$expr[names(keys), ]
gse47584$expr = expr
gse47584$keys = keys
checkDataObject(gse47584, "Dataset") # True!


# Run meta-analysis
xxyDatasets = list(GSE42331 = gse42331_xy, GSE47584 = gse47584)
xxyMeta = runMetaAnalysis(list(originalData = xxyDatasets), runLeaveOneOutAnalysis = F)

# Examine gene Significance
pooledResults = xxyMeta$metaAnalysis$pooledResults
#View(pooledResults)

# Significant by both effect size FDR and Fisher FDR
# I manually checked Fisher's FDR
sigGenes = subset(pooledResults, effectSizeFDR < 0.05)
nrow(sigGenes)

# Chromosome location
sort(ISEXSLoc_simple[rownames(sigGenes)])

#### Effect size plot!
# get dataset effect sizes
esPlotData= xxyMeta$metaAnalysis$datasetEffectSizes[rownames(sigGenes),]
esPlotData[is.na(esPlotData)] = 0 # convert NA to 0


# get sex labels
sexLabel = sexMetaObj$metaAnalysis$pooledResults[rownames(sigGenes), "effectSize"] > 0
sexLabel = ifelse(sexLabel, yes = 'Female-Associated', no = 'Male-Associated')
sexLabel = data.frame(`iSEXS.Group` = as.factor(sexLabel))
rownames(sexLabel) = rownames(esPlotData)

# Sex label color key
sexLabelColorKey = list()
sexLabelColorKey$iSEXS.Group = c(`Female-Associated` =femColor, `Male-Associated` =malColor)



### Set colors
paletteLength = 100
myBreaks <- c(seq(min(esPlotData), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(esPlotData)/paletteLength, max(esPlotData), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(c("purple", "white", "orange"))(paletteLength)

#### Proper heatmap

pdf('plots/klinefelters_heatmap.pdf')
pheatmap(esPlotData,
         cluster_cols = F, border_color = NA, treeheight_row = 0,
         annotation_row = sexLabel,
         breaks = myBreaks, 
         color = myColor, annotation_colors = sexLabelColorKey, annotation_names_row = F,
         fontsize = 16, fontsize_row = 20)
dev.off()


p = pheatmap(esPlotData,
         cluster_cols = F, border_color = NA, treeheight_row = 0,
         annotation_row = sexLabel,
         breaks = myBreaks, 
         color = myColor, annotation_colors = sexLabelColorKey, annotation_names_row = F,
         fontsize = 12, fontsize_row = 14, 
         cellheight = 20, cellwidth = 30)

xxy_heatmap = p[[4]]
dev.off()

#................................#
##### Klinefelter's Supp Fig #####
#................................#
# Purpose: 
#   - Create supp fig with extra XXY vs XY dataset GSE47584
#   - Show what genes are consistentl diff between XXY and XY


### GEt data from bonus dataset
gse47584$pheno$autoSig = calculateScore(autoSig, gse47584)
gse47584$pheno$xySig = calculateScore(xySig, gse47584)

# Autosomal-iSEXs is signficant (bartely)
t.test(gse47584$pheno$autoSig[gse47584$pheno$group == 'XXY Male'],
       gse47584$pheno$autoSig[gse47584$pheno$group == 'XY Male']) # p=0.048
t.test(gse47584$pheno$xySig[gse47584$pheno$group == 'XXY Male'],
       gse47584$pheno$xySig[gse47584$pheno$group == 'XY Male']) # p=0.012

suppXXY = list()


suppXXY$xySig = ggplot(gse47584$pheno, aes(x=group, y=xySig)) +
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5) + 
  scale_color_manual(values = c('grey47', malColor))+
  theme_cowplot(font_size = 18)+
  ggtitle("XY-iSEXS   GSE47584") + ylab("XY-iSEXS Score") + xlab('') +
  geom_signif(annotations = 'p=0.012', y_position = 1, textsize = 5, 
              xmin =1.1, xmax=1.9) + theme(legend.position = 'none')

suppXXY$autoSig = ggplot(gse47584$pheno, aes(x=group, y=autoSig)) +
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5) + 
  scale_color_manual(values = c('grey47', malColor))+
  theme_cowplot(font_size = 18)+
  ggtitle("Autosomal-iSEXS   GSE47584") + ylab("Autosomal-iSEXS Score") + xlab('') +
  geom_signif(annotations = 'p=0.048', y_position = 1.1, textsize = 5, 
              xmin =1.1, xmax=1.9) + theme(legend.position = 'none')


#.................................#
##### Demo XXY Genes Boxplots #####
#.................................#
# Purpose: 
#   - Make boxplots that showcase some of the sig
#     genes differentially expressed between XXY and XY

# Pick genes to showcase
quickViolin(gse42331, 'CD40LG', 'group')+ theme_cowplot() # Looking good!
quickViolin(gse42331, 'CA5BP1', 'group')+ theme_cowplot() # Looking good!
quickViolin(gse42331, 'MPO', 'group') # also nice!
quickViolin(gse42331, 'BPI', 'group')

### Get data for boxplots
myDataXXY = gse42331$pheno
expr = getSampleLevelGeneData(gse42331, c('CD40LG', 'CA5BP1', 'MPO','BPI'))
all(colnames(expr) == rownames(myDataXXY)) # True
myDataXXY = cbind(myDataXXY, t(expr))

xxyGenePlots =list()

### CD40LG
t.test(myDataXXY$CD40LG[myDataXXY$group == 'XX Female'],
       myDataXXY$CD40LG[myDataXXY$group == 'XXY Male']) # NS between XX and XXY
t.test(myDataXXY$CD40LG[myDataXXY$group == 'XY Male'],
       myDataXXY$CD40LG[myDataXXY$group == 'XXY Male'])# p=0.0063

xxyGenePlots$CD40LG = ggplot(myDataXXY, aes(x=group, y=CD40LG)) + 
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5, priority = 'random', size=1.5) + 
  scale_color_manual(values = c(femColor, 'grey48', malColor))+
  theme_cowplot(font_size = 16)+theme(legend.position = 'none') +
  ggtitle("GSE42331") + ylab("CD40LG") + xlab('') +
  geom_signif(annotations = 'N.S.', y_position = 8.3, textsize = 4.5, 
              xmin =1.1, xmax=1.9)+
  geom_signif(xmin=2.1, xmax=2.8, 
              annotations = 'p=0.0063', y_position = 8.3, textsize = 4.5) +
  scale_x_discrete(labels=c('XX Female' = paste('XX', 'Female', sep='\n'),
                            'XXY Male' = paste('XXY', 'Male', sep='\n'),
                            'XY Male' = paste('XY', 'Male', sep='\n')))


### CA5BP1
t.test(myDataXXY$CA5BP1[myDataXXY$group == 'XX Female'],
       myDataXXY$CA5BP1[myDataXXY$group == 'XXY Male']) # NS between XX and XXY
t.test(myDataXXY$CA5BP1[myDataXXY$group == 'XY Male'],
       myDataXXY$CA5BP1[myDataXXY$group == 'XXY Male'])# p=0.0057

xxyGenePlots$CA5BP1 = ggplot(myDataXXY, aes(x=group, y=CA5BP1)) + 
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5, priority = 'random', size=1.5) + 
  scale_color_manual(values = c(femColor, 'grey48', malColor))+
  theme_cowplot(font_size = 16)+theme(legend.position = 'none') +
  ggtitle("GSE42331") + ylab("CA5BP1") + xlab('') +
  geom_signif(annotations = 'N.S.', y_position = 8.19, textsize = 5, 
              xmin =1.1, xmax=1.9)+
  geom_signif(xmin=2.1, xmax=2.8, 
              annotations = 'p=0.0057', y_position = 8.19, textsize = 5) +
  scale_x_discrete(labels=c('XX Female' = paste('XX', 'Female', sep='\n'),
                            'XXY Male' = paste('XXY', 'Male', sep='\n'),
                            'XY Male' = paste('XY', 'Male', sep='\n')))


### MPO
t.test(myDataXXY$MPO[myDataXXY$group == 'XX Female'],
       myDataXXY$MPO[myDataXXY$group == 'XXY Male']) # NS between XX and XXY
t.test(myDataXXY$MPO[myDataXXY$group == 'XY Male'],
       myDataXXY$MPO[myDataXXY$group == 'XXY Male'])# p=0.011

xxyGenePlots$MPO = ggplot(myDataXXY, aes(x=group, y=MPO)) + 
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5, priority = 'random', size=1.5) + 
  scale_color_manual(values = c(femColor, 'grey48', malColor))+
  theme_cowplot(font_size = 16)+theme(legend.position = 'none') +
  ggtitle("GSE42331") + ylab("MPO") + xlab('') +
  geom_signif(annotations = 'N.S.', y_position = 7.7, textsize = 5, 
              xmin =1.1, xmax=1.9)+
  geom_signif(xmin=2.1, xmax=2.8, 
              annotations = 'p=0.011', y_position = 7.7, textsize = 5) +
  scale_x_discrete(labels=c('XX Female' = paste('XX', 'Female', sep='\n'),
                            'XXY Male' = paste('XXY', 'Male', sep='\n'),
                            'XY Male' = paste('XY', 'Male', sep='\n')))

### BPI
t.test(myDataXXY$BPI[myDataXXY$group == 'XX Female'],
       myDataXXY$BPI[myDataXXY$group == 'XXY Male']) # NS between XX and XXY
t.test(myDataXXY$BPI[myDataXXY$group == 'XY Male'],
       myDataXXY$BPI[myDataXXY$group == 'XXY Male'])# p=0.0023

xxyGenePlots$BPI = ggplot(myDataXXY, aes(x=group, y=BPI)) + 
  geom_boxplot() + geom_beeswarm(aes(col=group), cex=1.5, priority = 'random', size=1.5) + 
  scale_color_manual(values = c(femColor, 'grey48', malColor))+
  theme_cowplot(font_size = 16)+theme(legend.position = 'none') +
  ggtitle("GSE42331") + ylab("BPI") + xlab('') +
  geom_signif(annotations = 'N.S.', y_position = 9.2, textsize = 5, 
              xmin =1.1, xmax=1.9)+
  geom_signif(xmin=2.1, xmax=2.8, 
              annotations = 'p=0.0023', y_position = 9.2, textsize = 5) +
  scale_x_discrete(labels=c('XX Female' = paste('XX', 'Female', sep='\n'),
                            'XXY Male' = paste('XXY', 'Male', sep='\n'),
                            'XY Male' = paste('XY', 'Male', sep='\n')))





#...................................#
##### Assemble Supplemental XXY #####
#...................................#

### Assemble plot
suppFig_XXY_top = plot_grid(suppXXY$xySig, suppXXY$autoSig, xxy_heatmap, 
                        nrow= 1, ncol=3, labels = LETTERS[1:3])
suppFig_XXY_bottom = plot_grid(plotlist = xxyGenePlots, nrow=1, ncol=4, labels = LETTERS[4:7])
suppFig_XXY = plot_grid(suppFig_XXY_top, suppFig_XXY_bottom, nrow=2, ncol=1)


save_plot(filename = 'plots/supp2_XXY.pdf', suppFig_XXY, ncol=3, nrow=2, base_height = 4.3)


#..................................#
##### immunoStates ES Heatmaps #####
#..................................#
# Purpose:
#   - Create heatmaps of the immunoStates effect sizes for 
#     each group of iSEXS genes
# 
# iSEXS groups:
#   - Female X
#   - Female Auto
#   - Male Y
#   - Male X
#   - Male Autosome

# Load FRancesco's mid-range immunoStates matrix
load("/labs/khatrilab/ebongen/sexDifferences2.0/isexs/reference/immunoStates_es.RData")
isMatrix = ISfreshMqnESmidCell

# Create lists of interesting genes
gene_femX = xySig$posGeneNames
gene_femA = autoSig$posGeneNames
gene_malA = autoSig$negGeneNames
gene_malY = xySig$negGeneNames[ISEXSLoc[xySig$negGeneNames] == 'Y']
gene_malX = xySig$negGeneNames[ISEXSLoc[xySig$negGeneNames] != 'Y']

# Remove genes that are missing in isMatrix
gene_femA = gene_femA[ gene_femA %in% isMatrix$gene]


# Extract immunoStates ES for each iSEXS gene
isexsCellTypes = subset(isMatrix, gene %in% isexs)

# Create a matrix where rows are ISDS genes, column are cell types
esCellTypes = isexsCellTypes[,c("gene", "g", "cellType"), with =F]
esCellTypes = dcast(data = esCellTypes, formula = gene ~ cellType, value.var = "g")
myRowNames= esCellTypes$gene
esCellTypes$gene = NULL
esCellTypes = as.matrix(esCellTypes)
rownames(esCellTypes) = myRowNames


# Set max values to 1
esCellTypes[esCellTypes > 1 & !is.na(esCellTypes)] = 1
esCellTypes[esCellTypes < -1 & !is.na(esCellTypes)] = -1

### Pallete and color information
paletteLength = 11
myBreaks <- c(seq(min(esCellTypes, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(esCellTypes, na.rm = T)/paletteLength, max(esCellTypes, na.rm = T), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(c("purple", "white", "orange"))(paletteLength)

### Get Female genes 
femGenes = c(reorderGenesByCluster(gene_femX, esCellTypes), 
             reorderGenesByCluster(gene_femA, esCellTypes))
femGapsAt = rep(length(gene_femX), 4)

### Get male genes
malGenes = c(reorderGenesByCluster(gene_malY, esCellTypes), 
             reorderGenesByCluster(gene_malX, esCellTypes), 
             reorderGenesByCluster(gene_malA, esCellTypes))
malGapsAt = c(rep(length(gene_malY), 4), rep((length(gene_malY) + length(gene_malX)), 4))

# Left-hand female column
p1 = pheatmap(esCellTypes[femGenes,names(cellTypesPrettyNames)],
              cluster_cols = F, cluster_rows = F,
              border_color = NA, treeheight_row = 0,
              breaks = myBreaks, 
              color = myColor,
              cellwidth = 15, cellheight = 3,
              fontsize = 14, show_rownames = F,
              show_colnames = T, legend=T, 
              gaps_row = femGapsAt,
              labels_col = c("Hematopoietic Progenitor", "Neutrophil",
                             "Basophil", "Eosinophil", "MAST Cell",
                             "CD14+ Monocyte","CD16+ Monocyte","Macrophage M0","Macrophage M1","Macrophage M2",
                             "mDC","pDC",
                             "Plasma Cell","Naive B Cell","Memory B Cell",
                             "CD56bright NK Cell" ,"CD56dim NK Cell",
                             expression(paste(gamma, delta, ' T cell')),
                             "CD8+ T Cell","CD4+ T Cell"))[[4]]

# Right-hand male column
p2 = pheatmap(esCellTypes[malGenes,names(cellTypesPrettyNames)],
              cluster_cols = F, cluster_rows = F,
              border_color = NA, treeheight_row = 0,
              breaks = myBreaks, 
              color = myColor,
              cellwidth = 15, cellheight = 3,
              fontsize = 14, show_rownames = F,
              show_colnames = T, legend=T, 
              gaps_row = malGapsAt,
              labels_col = c("Hematopoietic Progenitor", "Neutrophil",
                                          "Basophil", "Eosinophil", "MAST Cell",
                                          "CD14+ Monocyte","CD16+ Monocyte","Macrophage M0","Macrophage M1","Macrophage M2",
                                          "mDC","pDC",
                                          "Plasma Cell","Naive B Cell","Memory B Cell",
                                          "CD56bright NK Cell" ,"CD56dim NK Cell",
                                          expression(paste(gamma, delta, ' T cell')),
                                          "CD8+ T Cell","CD4+ T Cell"))[[4]]


fig4_cells_hmap = plot_grid(p1, p2, nrow=1, ncol=2, labels = c('b','c'))

fig4_cells_hmap = fig4_cells_hmap + draw_label('Female-Associated iSEXS Genes', x=0.25, y=0.98) +
  draw_label('Male-Associated iSEXS Genes', x=0.75, y=0.98) +
  draw_label('X-Chr', x=0.024, y=0.87) + 
  draw_label('Autosomal', x=0.036, y=0.6, angle=90) +
  draw_label('Y-Chr', x=0.522, y=0.815) +
  draw_label('X-Chr', x=0.522, y=0.75) +
  draw_label('Autosomal', x=0.5357, y=0.6, angle=90)


### Save the heatmaps as a pdf
save_plot('plots/immunoStates_heatmap.pdf', fig4_cells_hmap, nrow=1, ncol=2, base_height = 7, base_width = 6)

#................................................#
##### Statistical Enrichment of immune cells #####
#................................................#
# Purpose: 
#   - Attach statististics to observation of CD4's in females 
#     and myeloid cells in males

# Calculate enrichment stats
cellEnrich_fem = cellTypeEnrichment(sexMetaObj$filterResults$isexs$posGeneNames, isMatrix)
cellEnrich_mal = cellTypeEnrichment(sexMetaObj$filterResults$isexs$negGeneNames, isMatrix)

### Get CD4 stats
#   - Z-score is 2.34
z_cd4 = cellEnrich_fem$value[grepl(pattern = 'CD4_pos', x = cellEnrich_fem$celltype)] # z-score is 2.34
pnorm(-abs(z_cd4)) # p = 0.0098


fig_cellEnrich_fem = ggplot(cellEnrich_fem, aes(celltype, value)) + 
  geom_hline(yintercept = 1.645, linetype='dashed', colour='grey47')+
  geom_hline(yintercept = 2.328, linetype='dashed', colour='black')+
  geom_col(fill = femColor) + theme_cowplot() +
  xlab('') + ylab('immunoStates Enrichment') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('Female-Associated iSEXS') 





### Get myeloid stats - best os M0 macrophage
z_mph = cellEnrich_mal$value[grepl(pattern = 'm0', x=cellEnrich_mal$celltype)]
pnorm(-abs(z_mph)) # p= 0.042

# What about 2nd runner up: Basophil
z_baso = cellEnrich_mal$value[grepl(pattern = 'basophil', x=cellEnrich_mal$celltype)]
pnorm(-abs(z_baso)) # p= 0.098


fig_cellEnrich_mal = ggplot(cellEnrich_mal, aes(celltype, value)) + 
  geom_hline(yintercept = 1.645, linetype='dashed', colour='grey47')+
  geom_hline(yintercept = 2.328, linetype='dashed', colour='black')+
  geom_col(fill = malColor) + theme_cowplot() +
  xlab('') + ylab('immunoStates Enrichment') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('Male-Associated iSEXS') 


### Combine into a plot
fig_cellEnrich = plot_grid(fig_cellEnrich_fem, fig_cellEnrich_mal, nrow=1, ncol=2)
save_plot('plots/cellEnrichment.pdf',fig_cellEnrich, ncol=2, nrow=1, base_height = 6, base_width = 5.5)



#.........................................#
##### immunoStates - demo Forest Plot #####
#.........................................#
# Purpose: 
#   - Create a forest plot for immunoStates to be a demo
#   - Used to help explain immunoStates values

### Get data for forest plot
cd19Data = subset(isMatrix, gene == 'CD19')
myCells = c('CD14_positive_monocyte', 'neutrophil', 
            'CD4_positive_alpha_beta_T_cell', 'CD8_positive_alpha_beta_T_cell', 
            'naive_B_cell', "plasma_cell", "memory_B_cell")
cd19Data = subset(cd19Data, cellType %in% myCells)
cd19Data = cd19Data[order(cd19Data$g, decreasing = T)]

cellLabels = as.character(cd19Data$cellType)
cellLabels[cellLabels == 'CD4_positive_alpha_beta_T_cell'] = 'CD4+ T Cell'
cellLabels[cellLabels == 'CD8_positive_alpha_beta_T_cell'] = 'CD8+ T Cell'
cellLabels[cellLabels == 'CD14_positive_monocyte'] = 'CD14+ Monocyte'
cellLabels = gsub(pattern = '_', replacement = ' ', x = cellLabels)
cellLabels = sapply(cellLabels, simpleCap)
cellLabels = paste(cellLabels, 's', sep='')

### Forest plot
dev.off() # For some reason, recordPlot doesn't work, does this fix it?
par(mar=c(9,0,7,7), cex = 1)
rmeta::metaplot(mn = cd19Data$g, se = cd19Data$se.g, labels ='' ,
                xlab = 'Immune Cell Specificity Score', ylab = '', lwd =3, main = '')
cd19Plot = recordPlot()
dev.off()

#..................................#
##### immunoStates - Schematic #####
#..................................#
# Purpose: 
#   - Add immunoStates schematic in as an image


# can't look at it alone, just plots the hex code for each pixel
schematic = ggdraw() + draw_image('reference/superImage.001.jpeg')


#fish = plot_grid(p2, nrow = 1, ncol=1)
#fish = fish + draw_label(paste(expression(gamma), expression(delta), ' T cell', sep = ''))
#fish = fish+ draw_label(expression(paste(gamma, delta, ' T cell')), x = 0.75, y = 0.75)
#save_plot('plots/fish.pdf', fish)


# Combine with CD14 forest plot
fig4_top = plot_grid(schematic, cd19Plot, rel_widths = c(3,1))

### Add labels for CD14 forest plot
# CD14, neut, macrophage, CD4, CD8, B
cellLabel = paste('CD14+ Monocyte', '\n', 'Neutrophil', '\n',
                  'Macrophage', '\n', 'CD4+ T cell', '\n', 'CD8+ T cell', 
                  '\n', 'B cell', sep = '')

# Create Pretty labels!
cellLabels2 = c()
for(i in 1:length(cellLabels)){
  cellLabels2 = c(cellLabels2, cellLabels[i], '\n')
}
cellLabels2 = cellLabels2[1:length(cellLabels2)-1] # remove last \n

fig4_top = fig4_top + 
  draw_label('CD19', x = 0.87, y=0.9, fontface = 'bold', size=14)+
  draw_label(paste(cellLabels2, collapse = ''), x=0.75, y = 0.54, size=11)

save_plot('plots/immunoStates_schematic.pdf', fig4_top, base_height = 3, base_width = 12)

#................................#
##### 4 Timecourse of Scores #####
#................................#
# Purpose: 
#   - Combine autosomal and XY scores of multiple aging cohorts
#   - Plot them over age
# 
# Future: 
#   - Look at that Japanese dataset

# Load timecourse datasets
setwd("/labs/khatrilab/ebongen/sexDifferences2.0/")
load("0_datasets/2_age/postMenopause/gse58137/gse58137_timecourse.RData")
load("0_datasets/1_validation/gse21311/gse21311_email.RData")
load("0_datasets/1_validation/gse38484/gse38484_validation.RData")
setwd(myDir)

# Grab datasets with before and after 40 well populated
timeDatasets = list(gse21311 = gse21311_email, gse38484 = gse38484_healthy, gse58137 = gse58137gpl6)
timeColumn = "age"
groupColumn = "sex"

### XY plot
myScoreData = list()

for(ds in names(timeDatasets)){
  datasetObject = timeDatasets[[ds]]
  myData = datasetObject$pheno
  myData$score = calculateScore(sexMetaObj$filterResults$xySig, datasetObject, suppressMessages = TRUE)
  myData[["time"]] = myData[[timeColumn]]
  myData[["group"]] = myData[[groupColumn]]
  
  myData = subset(myData, age < 75)
  myData = subset(myData, age >= 18)
  
  myScoreData[[ds]] = myData
}


comboData = rbind(myScoreData[[1]][,c("age", "sex", "score")],
                  myScoreData[[2]][,c("age", "sex", "score")],
                  myScoreData[[3]][,c("age", "sex", "score")])


pointSize = 1
lineSize = 1

timecourseXY = ggplot() + 
  geom_point(data = myScoreData$gse21311, aes(x = age, y = score, col = sex), size=pointSize) +
  geom_point(data = myScoreData$gse38484, aes(x = age, y = score, col = sex), size = pointSize) +
  geom_point(data = myScoreData$gse58137, aes(x = age, y = score, col = sex), size = pointSize) +
  geom_smooth(data = comboData, aes(x = age, y= score, col = sex), size = lineSize) +
  scale_color_manual(values=c("#F42C04","#190E4F")) + ggtitle("GSE21311, GSE38484, GSE58137") +
  xlab("Age (years)") + ylab("XY-iSEXS Score") +theme_cowplot(font_size = 12) +
  theme(legend.position = "none")


### Autosomal score aging plot
myScoreData = list()

for(ds in names(timeDatasets)){
  datasetObject = timeDatasets[[ds]]
  myData = datasetObject$pheno
  myData$score = calculateScore(sexMetaObj$filterResults$autoSig, datasetObject, suppressMessages = TRUE)
  myData[["time"]] = myData[[timeColumn]]
  myData[["group"]] = myData[[groupColumn]]
  
  myData = subset(myData, age < 75)
  myData = subset(myData, age >= 18)
  
  myScoreData[[ds]] = myData
}


comboData = rbind(myScoreData[[1]][,c("age", "sex", "score")],
                  myScoreData[[2]][,c("age", "sex", "score")],
                  myScoreData[[3]][,c("age", "sex", "score")])

# Capitalize the sex labels for prettiness
comboData$sex = sapply(as.character(comboData$sex), simpleCap)


# Plot Autosomal-iSEXS score with age
timecourseAuto = ggplot(comboData, aes(x=age, y=score, col=sex)) +
    geom_smooth(size=2) + geom_point(size=1.5)+
    scale_color_manual(values=c(femColor,malColor), name='') + ylab("Autosomal-iSEXS Score") +
    xlab("Age (years)") + ggtitle("GSE58137, GSE21311, GSE38484") +theme_cowplot(font_size = 16)+
    theme(legend.position = c(0.8, 0.1))
  
  
# 435 individuals total
for(i in 1:length(myScoreData)){
  print(nrow(myScoreData[[i]]))
}

#..............................#
##### Monocytes Timecourse #####
#..............................#
# Purpose: 
#   - Make timecourse of Monocytes (% of CD45+ cells)
#   - Patin et al. data (Millieu Interier cohort)
#   - Also see if CD4's  could play a role in changing cell proportions
#
# Results; 
#   - Monocytes only have diffrence in proportions before 40yo, becaome same after that
#   - CD4's go low in males after 60 years old, could explain separation later
#
# Conclusion:
#   - Monocyte difference is bigger, so showcase that
#   - CD4's can go into supplementary 

# Load Patin data
library(mmi)

# Data frames of interest
dim(facs) # flow data
dim(ecrf)# demographic data

# Subjects are in the same order 
all(as.character(facs$SUBJID) == as.character(ecrf$SUBJID)) # True


# Get cell types of interest
cellTypes = colnames(facs)[2:ncol(facs)]
monoVariables = cellTypes[grepl(pattern = 'mono', x = cellTypes, ignore.case = T)]
monoVariables = c(cellTypes[grepl(pattern = "CD45", x = cellTypes)], monoVariables)
monoVariables = c(monoVariables, "N_CD4pos.panel1", "N_CD4pos_naive.panel1")

# Create data frame of data to plot
patin = facs[,c('SUBJID', monoVariables)]
all(patin$SUBJID == ecrf$SUBJID) # true
patin$sex = ecrf$Sex
patin$age = ecrf$Age

# Examine patin data frame and sanity check
table(patin$sex) # 399 males, 417 females
hist(patin$age) # good coverage of ages 20-70
all(patin$age == ecrf$Age) # true
all(colnames(patin$SUBJID) == colnames(ecrf$SUBJID)) # true

# 36 samples don't have flow measured
table(is.na(patin$N_mono.panel5), is.na(patin$N_CD45pos.panel5))
patin = subset(patin, !is.na(N_mono.panel5))
nrow(patin) # 780 people total
table(patin$sex)

### Monocytes ###

# Create figure5g for monocyte proportions with age
monoTimecourse = ggplot(patin, aes(x=age, y=(N_mono.panel5/N_CD45pos.panel5*100), col=sex)) + 
  geom_point(size = 1.5) + geom_smooth(size =1.5)+
  scale_color_manual(values = c(malColor, femColor), name ='') +theme_cowplot(font_size = 16) + 
  xlab("Age (years)") + ylab("Monocytes (% of CD45+ cells)") + ggtitle("Milieu Interieur Cohort - Flow Cytometry ") +
  ylim(c(1.2,12)) + theme(legend.position = c(0.86,0.1))


### CD4's ###

# Not different before 40, after 40 females go a little higher
# Not strong like monocytes are
ggplot(patin, aes(x=age, y=(N_CD4pos.panel1/N_CD45pos.panel5*100), col=sex)) + 
  geom_point(size = 2) + geom_smooth(size =2)+
  scale_color_manual(values = c(malColor, femColor)) +theme_cowplot()+
  ylim(c(0,50))

# Same pattern as total CD4's
ggplot(patin, aes(x=age, y=(N_CD4pos_naive.panel1/N_CD45pos.panel5*100), col=sex)) + 
  geom_point(size = 2) + geom_smooth(size =2)+
  scale_color_manual(values = c(malColor, femColor)) +theme_cowplot() + ylim(c(0, 30))


### Raw counts of CD4's in males decline with age, but less so in females 
ggplot(patin, aes(x=age, y=(N_CD4pos.panel1), col=sex)) + 
  geom_point(size = 2) + geom_smooth(size =2)+
  scale_color_manual(values = c(malColor, femColor)) +theme_cowplot()
ggplot(patin, aes(x=age, y=(N_CD4pos_naive.panel1), col=sex)) + 
  geom_point(size = 2) + geom_smooth(size =2)+
  scale_color_manual(values = c(malColor, femColor)) +theme_cowplot()


cor.test(patin$N_CD4pos.panel1, patin$N_CD4pos_naive.panel1) # r = 0.78


## Monocytes in young people
patin_young = subset(patin, age <=40)
range(patin_young$age) # 20-40
patin_young$monoProp = patin_young$N_mono.panel5/patin_young$N_CD45pos.panel5
ggplot(patin_young, aes(x=sex, y=monoProp)) + geom_violin(trim=F, fill='grey')+
  geom_jitter(width=0.1, size=2)
t.test(patin_young$monoProp[patin_young$sex == 'Male'],
       patin_young$monoProp[patin_young$sex == 'Female']) #

#.....................#
##### 2-way Anova #####
#.....................#
# Purpose:
#   - Do 2-way ANOVA to look at how age group and sex affect iSEXS/monocytes

# Prepare patin's data for ANOVA
monoAnova = patin
monoAnova = subset(monoAnova, age <=40 | age >=50)
hist(monoAnova$age) # No one in their 40s
monoAnova$monoProp = monoAnova$N_mono.panel5/monoAnova$N_CD45pos.panel5
monoAnova$ageGroup = factor(ifelse(monoAnova$age > 45, yes = 'Older', no = 'Younger'))
levels(monoAnova$sex)


# Perform ANOVA on monocytes
# Age, sex, and interaction of age and sex significant
colnames(monoAnova)
mono_anova <- aov(monoProp ~ sex * ageGroup, data = monoAnova)
Anova(mono_anova, type = "III")


### Autosomal-iSEXS score
myScoreData = list()

for(ds in names(timeDatasets)){
  datasetObject = timeDatasets[[ds]]
  myData = datasetObject$pheno
  myData$score = calculateScore(sexMetaObj$filterResults$autoSig, datasetObject, suppressMessages = TRUE)
  myData[["time"]] = myData[[timeColumn]]
  myData[["group"]] = myData[[groupColumn]]
  
  myData = subset(myData, age < 75)
  myData = subset(myData, age >= 18)
  
  myScoreData[[ds]] = myData
}

# Get data needed for Autosomal-iSEXS Anova
autoAnova = rbind(myScoreData[[1]][,c("age", "sex", "score")],
                  myScoreData[[2]][,c("age", "sex", "score")],
                  myScoreData[[3]][,c("age", "sex", "score")])

# Make sure the data's what we think it is - yup, it is!
ggplot(autoAnova, aes(x=age, y=score, col=sex)) + 
  geom_point() + geom_smooth() + scale_color_manual(values = c(femColor, malColor))



# Subset to remove people in their 40s
autoAnova = subset(autoAnova, age <=40 | age >=50)
hist(autoAnova$age)# Nobody in their 40s

# Create age group
autoAnova$ageGroup = factor(ifelse(autoAnova$age <=40, yes = 'Younger', no = 'Older'))
table(autoAnova$ageGroup, autoAnova$age <=40)

# Anova results
# Sex and Age individually are significant
# The interaction term is not significant
auto_anovaResults <- aov(score ~ sex * ageGroup, data = autoAnova)
Anova(auto_anovaResults, type = "III")

# Exlore what non-significance means
autoAnova$group = paste(autoAnova$ageGroup, autoAnova$sex)

ggplot(autoAnova, aes(x=group, y=score)) +
  geom_boxplot() +theme_cowplot()

# Probably because they're both still significant. 
# Okay, so males and females are always statistically different
# But, younger females are sig higher than older females
# There's no difference between older and younger males
t.test(autoAnova$score[autoAnova$group == 'Older female'],
       autoAnova$score[autoAnova$group == 'Older male'])# p = 0.00126

t.test(autoAnova$score[autoAnova$group == 'Younger female'],
       autoAnova$score[autoAnova$group == 'Younger male']) # p =8.32E-8

t.test(autoAnova$score[autoAnova$group == 'Older female'],
       autoAnova$score[autoAnova$group == 'Younger female']) # p=4.08E-5

t.test(autoAnova$score[autoAnova$group == 'Older male'],
       autoAnova$score[autoAnova$group == 'Younger male'])# p = 0.15


# What about monocytes? 
monoAnova$group = paste(monoAnova$ageGroup, monoAnova$sex)
t.test(monoAnova$monoProp[monoAnova$group == 'Older Female'],
       monoAnova$monoProp[monoAnova$group == 'Older Male']) #p=0.028
t.test(monoAnova$monoProp[monoAnova$group == 'Younger Female'],
       monoAnova$monoProp[monoAnova$group == 'Younger Male']) #p=2E-11
#.............................#
##### Fig 4 - ASSEMBLE!!! #####
#.............................#

# Get bottom row
fig4_bottom = plot_grid(monoTimecourse, timecourseAuto,
                 nrow=1, ncol = 2)

save_plot('plots/fig4de_monoTimecourse.pdf', plot=fig4_bottom, ncol=2, nrow=1, base_width = 6, base_height = 3.75)

fig4_grid = plot_grid(fig4_top, fig4_cells_hmap, fig4_bottom, nrow = 3, ncol = 1, rel_heights = c(1,3,1.5))
save_plot('plots/fig4_immunoStates.pdf', plot = fig4_grid, ncol = 2, nrow = 3, base_width = 5.5, base_height = 4.5)


#.............................................#
##### Correlate Monocytes with Auto-iSEXS #####
#.............................................#
# Purpose: 
#   - Correlate Monocyte proportion with Autosomal-iSEXS score
#
# Background:
#    - gse65133: Newman et al. 
#         - Flow on 20 healthy ~young people
#         - According to supplemental methods, monocytes were identified forward/side scatter
#         - They report %, so I'm assuming % of CD45+ live singlets, but I can't find
#           any sentences explicity saying that
#    - gse47353: Tsang et al. 
#         - It's a Discovery dataset, so hard to interpret
#         - % CD14+ of viable CD45+ cells (Total Monocytes)
#         - After digging around in their supplemental information, each
#           cell type is given as the % of their parent population. 
#         - So, Total Monocytes = % CD14+ of CD45+ cells
# Load cohorts
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')
load("0_datasets/4_sortedCells/cellProp/gse65133/gse65133.RData")
load('0_datasets/4_sortedCells/cellProp/SDY80/gse47353.RData')
setwd(myDir)

# Get cell type columns
myColumns6 = colnames(gse65133$pheno)[14:22]
myColumns4 = colnames(gse47353_baseline$pheno[30:58])


### Get data for gse65133 - CIBERSORT paper data
myData6 = gse65133$pheno
myData6$scoreXY = calculateScore(sexMetaObj$filterResults$xySig, gse65133)
myData6$scoreAuto = calculateScore(sexMetaObj$filterResults$autoSig, gse65133)

# Plot with autosomal score
autoStats = cor.test(myData6$scoreAuto, myData6$Monos)
autoLabel = paste("r=", signif(autoStats$estimate,2), "\n", "p=", signif(autoStats$p.value,2), sep='')

# Make sex be capitalized
myData6$sex = sapply(myData6$sex, simpleCap)

# List of supFig2 figures
sFig2_list = list()


sFig2_list$c_gse65133 = ggplot(myData6, aes(x=Monos, y=scoreAuto, col=sex)) + geom_smooth(method='lm',formula=y~x, colour='black') +
  geom_point(size=3)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes (% of Mononuclear Cells)") +
  ylab("Autosomal-iSEXS Score") +
  draw_label(label = autoLabel, x= 30, 0.5) +
  ggtitle("GSE65133 - PBMCs") +
  theme_cowplot()+
  theme(legend.position=c(0.08, 0.18))+
  labs(color='Sex')


# Plot Monocytes vs Autosomal score, without outlier!
myData6_noOutlier = subset(myData6, Monos != max(myData6$Monos, na.rm = T))

noOutlierTest = cor.test(myData6_noOutlier$Monos, myData6_noOutlier$scoreAuto)
noOutlierLabel = paste('r=', signif(noOutlierTest$estimate,2), '\n', 'p=', signif(noOutlierTest$p.value,2), sep='')

ggplot(myData6_noOutlier, aes(x=Monos, y=scoreAuto, col=sex)) + geom_smooth(method='lm',formula=y~x, colour='black') +
  geom_point(size=3)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes (% of CD45+ Cells)") +
  ylab("Autosomal-iSEXS Score") +
  draw_label(label = noOutlierLabel, x= 20, 1) +
  ggtitle("GSE65133 - PBMCs - Outlier Removed") +
  theme_cowplot()+
  theme(legend.position=c(0.08, 0.18))+
  labs(color='Sex')

### Plots form validation cohort
myData4 = gse47353_baseline$pheno
myData4$sex = sapply(myData4$sex, function(x) simpleCap(x))
myData4$scoreXY = calculateScore(sexMetaObj$filterResults$xySig,gse47353_baseline)
myData4$scoreAuto = calculateScore(sexMetaObj$filterResults$autoSig, gse47353_baseline)

# Get stats
autoStats4 = cor.test(myData4$scoreAuto, myData4$Total_Monocytes)
autoLabel4 = paste("r=", signif(autoStats4$estimate,2), "\n", "p=", 
                   signif(autoStats4$p.value,2), sep='')


sFig2_list$d_gse47353 = ggplot(myData4, aes(x=Total_Monocytes, y=scoreAuto, col=sex)) + geom_smooth(method='lm',formula=y~x, colour='black') +
  geom_point(size=3)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes (% of CD45+ Cells)") +
  ylab("Autosomal-iSEXS Score") +
  draw_label(label = autoLabel4, x= 35, 0.5) +
  ggtitle('GSE47353 - PBMCs') +
  theme_cowplot() +
  theme(legend.position = c(0.08, 0.18)) +
  labs(color ='Sex')
#...........................#
##### Gabi's Cytof Data #####
#...........................#
# Purpose: 
#    - Show that according to Cytof, there's also no sex diff in monocytes in old people

# Load Gabi's data:
gabi = read.delim('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/4_sortedCells/cellProp/gabi/cytof_cellProp.txt')

# Rearrange group, so that it's young first, then older
group = as.character(gabi$group)
group = factor(group, levels=unique(group)[c(2,3,1,4)])
gabi$group = group

# Create group2
group2 = gsub(pattern = ' ', replacement = '\n', x = as.character(gabi$group))
gabi$group2 = factor(group2,levels = unique(group2)[c(2,3,1,4)] )

# Get stats
youngPval = t.test(gabi$monocyte[gabi$group == 'Younger Female'],
                   gabi$monocyte[gabi$group == 'Younger Male'])$p.value
oldPval = t.test(gabi$monocyte[gabi$group == 'Older Female'],
                 gabi$monocyte[gabi$group == 'Older Male'])$p.value

sFig2_list$a_gabi = ggplot(gabi, aes(x=group2, y=monocyte*100)) + 
  geom_violin(trim=F, fill='grey') + geom_jitter(width=0.1, size=2, aes(col=sex))+
  scale_color_manual(values = c(femColor, malColor)) + ylab('Monocytes (% of Mononuclear Cells)') +
  xlab('') + theme_cowplot() + theme(legend.position = 'none') + 
  geom_signif(annotations = paste('p =', signif(youngPval,2)), y_position = 40, xmin= 1.2, xmax=1.8, textsize=4)+
  geom_signif(annotations = paste('p =', signif(oldPval,2)), y_position = 30, xmin= 3.2, xmax=3.8, textsize=4) +
  ggtitle('CyTOF - Fragiadakis & Bjornson et al.') 


# But what about young females vs older females? - No difference
t.test(gabi$monocyte[gabi$group == 'Younger Female'],
       gabi$monocyte[gabi$group == 'Older Female' & gabi$age >=50])

# Females over 50 do not cluster at the top of the violin
ggplot(gabi, aes(x=group2, y=monocyte)) + 
  geom_violin(trim=F, fill='grey') + 
  geom_jitter(width=0.1, size=2, aes(col=age >=50))

#.............................................#
##### SDY212 - Validate againg Auto-iSEXS #####
#.............................................#
# Purpose:
#   - Mark Davis HIPC vaccination, from 2008
#   - Show Auto-iSEXS only different between young male and female
#   - Younger is 20-30 years old
#   - Older is 60-90 years old

load('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/2_age/postMenopause/sdy212/sdy212.RData')

# Prep pheno for plotting
sdy212$pheno$autoScore = calculateScore(autoSig, sdy212)
sdy212$pheno$group = sapply(sdy212$pheno$group, simpleCap)
sdy212$pheno$group = factor(sdy212$pheno$group, levels = unique(sdy212$pheno$group))

# Create pretty group
group = sdy212$pheno$group
group = gsub(pattern = ' ', replacement = '\n ', x = group)
sdy212$pheno$group2 = factor(group, levels = unique(group))

# Get p-values
youngPval = t.test(sdy212$pheno$autoScore[sdy212$pheno$group == 'Younger Female'],
       sdy212$pheno$autoScore[sdy212$pheno$group == 'Younger Male'])$p.value
oldPval = t.test(sdy212$pheno$autoScore[sdy212$pheno$group == 'Older Female'],
                   sdy212$pheno$autoScore[sdy212$pheno$group == 'Older Male'])$p.value

sFig2_list$b_sdy212 = ggplot(sdy212$pheno, aes(x=group2, y=autoScore)) + 
  geom_violin(fill='grey', trim=F) + geom_jitter(width=0.1, size=2, aes(col=sex)) +
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Autosomal-iSEXS Score') + theme(legend.position = 'none')+
  theme_cowplot() + 
  geom_signif(annotations = paste('p =', signif(youngPval,2)), y_position = 3.8, xmin= 1.2, xmax=1.8, textsize=4)+
  geom_signif(annotations = paste('p =', signif(oldPval,2)), y_position = 2.2, xmin= 3.2, xmax=3.8, textsize=4) +
  theme(legend.position = 'none') + ggtitle('SDY212 - Whole Blood')



#.............................#
##### Asssemble sFig2 !!! #####
#.............................#
# Purpose: 
#   - Make  monocytes and aging supp fig!

sfig2_monoAge = plot_grid(sFig2_list$a_gabi, 
                          sFig2_list$b_sdy212,
                          sFig2_list$c_gse65133, sFig2_list$d_gse47353,
                          nrow = 2, ncol=2, labels = LETTERS[1:4])
save_plot('plots/supp3_monocytes.pdf', sfig2_monoAge, ncol=2, nrow=2)


#......................................................#
##### sFig3: Genes that change with age in Females #####
#......................................................#
# Purpose: 
#   - Are the genes that increase with age in females assoc. with monocytes?
#   - Would support the "decrease in iSEXS is due to monocytes proportion" idea

# Get the significant genes from agingMetaObj
oldGenes = c(agingMetaObj$filterResults$FDR0.05_es0_nStudies2_looaFALSE_hetero0$posGeneNames,
             agingMetaObj$filterResults$FDR0.05_es0_nStudies2_looaFALSE_hetero0$negGeneNames)
length(oldGenes) # 6 genes

# All the genes in the datasets of intrest?
all(oldGenes %in% gse65133$keys) # True
all(oldGenes %in% gse47353_baseline$keys) # True

# Get data to correlate monocytes with
myData4 = gse47353_baseline$pheno
myData6 = gse65133$pheno
for(myGene in oldGenes){
  myData4[[myGene]] = unlist(getSampleLevelGeneData(gse47353_baseline, myGene))
  myData6[[myGene]] = unlist(getSampleLevelGeneData(gse65133, myGene))
}


# Get a list of plots
oldLadyPlots = list()
for(myGene in oldGenes){
  # Forest plot
  oldLadyPlots[[paste(myGene, '_fp', sep='')]] = quickForestPlot(agingMetaObj, myGene, fontSize = 12)
  
  # Correlation with gse4
  myData4$myGene = myData4[,myGene]
  oldLadyPlots[[paste(myGene,'_gse4', sep='')]] = 
    ggplot(myData4, aes(x=Total_Monocytes, y=myGene)) + 
    geom_smooth(method='lm', col='black') + geom_point(aes(col=sex), size=2) +
    theme_cowplot(font_size = 12) + scale_color_manual(values = c(femColor, malColor)) +
    xlab('Monocytes (% of CD45+ Mononuclear Cells)') + ggtitle('GSE47353') +
    theme(legend.position = 'none') +ylab(myGene)
  
  myData6$myGene = myData6[,myGene]
  oldLadyPlots[[paste(myGene, '_gse6', sep='')]] = 
    ggplot(myData6, aes(x=Monos, y=myGene))+
    geom_smooth(method='lm', col='black') + geom_point(aes(col=sex), size=2) +
    theme_cowplot(font_size = 12) + scale_color_manual(values = c(femColor, malColor)) +
    xlab('Monocytes (% of Mononuclear Cells)') + ggtitle('GSE65133') +
    theme(legend.position = 'none') +ylab(myGene)
}

oldLadyGrid = plot_grid(plotlist = oldLadyPlots, ncol=3)

save_plot('plots/oldLady_fp_monoCor.pdf', oldLadyGrid, ncol=3, nrow=length(oldGenes))



### ZNF827
oldLadyPlots$ZNF827_fp = quickForestPlot(agingMetaObj, 'ZNF827')

cor.test(myData4$Total_Monocytes, myData4$ZNF827)

ggplot(myData4, aes(x=Total_Monocytes, y=ZNF827)) + 
  geom_smooth(method='lm', col='black') + geom_point(aes(col=sex), size=2) +
  theme_cowplot() + scale_color_manual(values = c(femColor, malColor)) +
  xlab('Monocytes (% of CD45+ Mononuclear Cells)') + ggtitle('GSE47353') +
  theme(legend.position = 'none')

ggplot(myData6, aes(x=Monos, y=ZNF827))+
  geom_smooth(method='lm', col='black') + geom_point(aes(col=sex), size=2) +
  theme_cowplot() + scale_color_manual(values = c(femColor, malColor)) +
  xlab('Monocytes (% of Mononuclear Cells)') + ggtitle('GSE65133') +
  theme(legend.position = 'none')

#.........................................#
##### 6 GSE68310 Auto-iSEXS vs Titers #####
#.........................................#
# Purpose: 
#   - plot baseline auto-iSEXS vs change in Ab titers

# Load Dataset
load('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/6_infection/infection/gse68310_WT_flu/gse68310_flu.RData')


# Get data
myDataFlu = gse68310Flu$pheno
myDataFlu$xyScore = calculateScore(xySig, gse68310Flu)
myDataFlu$autoScore = calculateScore(autoSig, gse68310Flu)
myDataFlu$Sex = sapply(myDataFlu$sex, simpleCap)

# Get baseline baseline data
myDataFluBase = subset(myDataFlu, time_point %in% c('Baseline'))
table(myDataFluBase$sex, myDataFluBase$time_point) # Males and females, only baseline

# Make xlab with delta
delta_xlab = expression(paste(Delta, 'Titers (anti-H1N1 Ab)'))


fig6a_titers= ggplot(myDataFluBase, aes(x=springTiter-baseTiter, y=autoScore, col=Sex)) +
  geom_smooth(method = 'lm')+
  geom_point(size=2.5) + theme_cowplot() +
  scale_color_manual(values = c(femColor, malColor)) +
  #ylab('Autosomal-iSEXS at Baseline') + xlab('Spring Titers - Baseline Titers (anti-H1N1 Ab)') +
  ylab('Autosomal-iSEXS at Baseline') + xlab(delta_xlab) +
  draw_label('p<0.001, r=0.74', x=5, y=-1, colour=malColor) +
  draw_label('N.S.', x=1, y=1.3, colour = femColor) +
  theme(legend.position = c(0.75, 0.15))

# Calculate the correlation stats
myDataFluBase$diff = myDataFluBase$springTiter - myDataFluBase$baseTiter
cor.test(myDataFluBase$autoScore[myDataFluBase$sex == 'male'],
         myDataFluBase$diff[myDataFluBase$sex == 'male']) # r = 0.74, p = 0.00064
cor.test(myDataFluBase$autoScore[myDataFluBase$sex == 'female'],
         myDataFluBase$diff[myDataFluBase$sex == 'female']) # zilch


#### Make a ROC plot for it
# Use > instead of >= to make the cutoff lineup wiht male median
# 1 = high resopnder, 0 = low responder
titerClass = ifelse(myDataFluBase$diff > median(myDataFluBase$diff), yes = 1, no = 0)
myDataFluBase$titerClass = titerClass

# subset by sexes
myDataFluBase_mal = subset(myDataFluBase, sex == 'male')
myDataFluBase_fem = subset(myDataFluBase, sex == 'female')

# Run ROC calculation
rocResults = list()
rocResults$male = calculateROC(myDataFluBase_mal$titerClass, myDataFluBase_mal$autoScore)
rocResults$female = calculateROC(myDataFluBase_fem$titerClass, myDataFluBase_fem$autoScore)

# reorganize the ploting data
rocResults$female$roc = cbind(rocResults$female$roc, sex = rep('female', nrow(rocResults$female$roc)))
rocResults$male$roc = cbind(rocResults$male$roc, sex = rep('male', nrow(rocResults$male$roc)))
rocPlotData = as.data.frame(rbind(rocResults$male$roc, rocResults$female$roc))

# Get labels
malLabel = paste('Males AUC=', signif(rocResults$male$auc,2), ' (95% CI ',signif(rocResults$male$auc.CI,2)[1],
                 '-', signif(rocResults$male$auc.CI,2)[2],')', sep='' )
femLabel = paste('Females AUC=', signif(rocResults$female$auc,2), ' (95% CI ',signif(rocResults$female$auc.CI,2)[1],
                 '-', signif(rocResults$female$auc.CI,2)[2],')', sep='' )

# ROC plot!
abROC = ggplot(rocPlotData, aes(x=x, y=y, col=sex))+
  geom_abline(slope=1, intercept=0, linetype= 'longdash', colour = 'darkgrey', size=1.2)+
  geom_line(size=1.2)+
  xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")+
  theme_cowplot() +
  scale_color_manual(values = c(malColor, femColor), labels = c(malLabel, femLabel), name = '')+
  theme(legend.position = c(0.1,0.15))

#.........................................#
##### 6 Auto-iSEXS Baseline vs Spring #####
#.........................................#
# Purpose: 
#   - Create correlation plots showing how consistent
#     Autosomal-iSEXS is 
#   - Show that hey, people have a pretty consistent baseline Auto-iSEXS

# Get Autosomal score form gse68310
fluBaseSpring = gse68310Flu$pheno
fluBaseSpring$autoScore = calculateScore(autoSig, gse68310Flu)

# Subset to only baseline and spring
unique(fluBaseSpring$time_point)
fluBaseSpring = subset(fluBaseSpring, time_point %in% c('Baseline', 'Spring'))

# dcast to make it autosomal scores at baseline and spring
fluBaseSpring = data.table::dcast(fluBaseSpring, subject_id + sex ~ time_point, value.var = 'autoScore')

# Subset to males and females 
fluBaseSpring_fem = subset(fluBaseSpring, sex == 'female')
fluBaseSpring_mal = subset(fluBaseSpring, sex == 'male')

# Get stats for males and females
malTest = cor.test(fluBaseSpring_mal$Baseline, fluBaseSpring_mal$Spring)
femTest = cor.test(fluBaseSpring_fem$Baseline, fluBaseSpring_fem$Spring)

# Extract stats to make captions
malLabel = paste('r=', signif(malTest$estimate,2), '  p=', signif(malTest$p.value,2), sep='')
femLabel = paste('r=', signif(femTest$estimate,2), '  p=', signif(femTest$p.value,2), sep='')

# Plot of females in spring vs Baseline
femSpring = ggplot(fluBaseSpring_fem, aes(x=Baseline, y=Spring)) +
  geom_smooth(method = 'lm', colour=femColor)+
  geom_point(size=2.5, colour=femColor) +
  theme_cowplot() + 
  xlab('Baseline Autosomal-iSEXS Score') +
  ylab(paste('Post-Flu Season', '\n', 'Autosomal-iSEXS Score')) + 
  ggtitle('Females') +
  draw_label(femLabel, x=1,y=-1, colour=femColor)


malSpring = ggplot(fluBaseSpring_mal, aes(x=Baseline, y=Spring)) +
  geom_smooth(method = 'lm', colour=malColor)+
  geom_point(size=2.5, colour=malColor) +
  theme_cowplot() + 
  xlab('Baseline Autosomal-iSEXS Score') +
  ylab(paste('Post-Flu Season', '\n', 'Autosomal-iSEXS Score')) + 
  ggtitle('Males')+
  draw_label(malLabel, x=0.45,y=-1.5, colour=malColor)



#................................#
##### Assemble 6! Flu Titers #####
#................................#
fig6 = plot_grid(fig6a_titers, abROC, malSpring, femSpring, nrow = 2, ncol = 2, labels = LETTERS[1:4])

save_plot('plots/fig6_fluTiters.pdf', fig6, ncol=2, nrow=2, base_height = 4.11)

#....................................#
##### Auto-iSEXS and Vaccination #####
#....................................#
# Purpose: 
#   - Show relationship between Autosomal-iSEXS and vaccine response

# Load Datasets
stevenData = readRDS(file = '/labs/khatrilab/ebongen/friends/steven/fluCombo.rdat')

### Males ###

# Get Male data
gse48018 = stevenData$GSE48018

# Subset to baseline
phenoBase = subset(gse48018$pheno, timepoint == 0)
dim(phenoBase) # 111 rows
gse48018_baseline = subsetGEMFromPheno(gse48018, phenoBase)


# Get data for plotting
vaccMal = gse48018_baseline$pheno
vaccMal$autoScore = calculateScore(autoSig, gse48018_baseline)
vaccMal$deltaTiter = vaccMal$D28 - vaccMal$D0
vaccMal = subset(vaccMal, !is.na(deltaTiter))

# Get stats
malStats = cor.test(vaccMal$deltaTiter, vaccMal$autoScore)
malLabel = paste('r=', signif(malStats$estimate, 2), ', p=', signif(malStats$p.value,2), sep='')

# Plot for males!
vaccinePlots = list()

vaccinePlots$males = ggplot(vaccMal, aes(x=deltaTiter, y=autoScore)) +
  geom_smooth(method = 'lm', colour=malColor)+
  geom_point(size=2, colour=malColor) + ggtitle('GSE48018 - Males') +
  theme_cowplot() + ylab('Autosomal-iSEXS Score') + xlab(expression(paste(Delta, 'Titers'))) +
  draw_label(malLabel, x=2000, y=-3.2, colour = malColor)


### Females
gse48023 = stevenData$GSE48023

phenoBase = subset(gse48023$pheno, timepoint == '0')
gse48023_baseline = subsetGEMFromPheno(gse48023, phenoBase)

# Cloud of values, all female
quickScatter(gse48023_baseline, 'RPS4Y1', 'XIST', 'seroClass')

# Get data for plotting
vaccFem = gse48023_baseline$pheno
vaccFem$autoScore = calculateScore(autoSig, gse48023_baseline)
vaccFem$deltaTiter = vaccFem$D28 - vaccFem$D0
vaccFem = subset(vaccFem, !is.na(deltaTiter))

# Get stats
femStats = cor.test(vaccFem$deltaTiter, vaccFem$autoScore)
femLabel = paste('r=', signif(femStats$estimate, 2), ', p=', signif(femStats$p.value,2), sep='')


vaccinePlots$females = ggplot(vaccFem, aes(x=deltaTiter, y=autoScore)) +
  geom_smooth(method = 'lm', colour=femColor)+
  geom_point(size=2, colour=femColor) + ggtitle('GSE48023 - Females') +
  theme_cowplot() + ylab('Autosomal-iSEXS Score') + xlab(expression(paste(Delta, 'Titers'))) +
  draw_label(femLabel, x=200, y=-3, colour = femColor)


# Assemble!
sfig3_vacc = plot_grid(plotlist = vaccinePlots, nrow = 1, ncol=2, labels = LETTERS[1:2])
save_plot('plots/supp3_vaccine.pdf', sfig3_vacc, ncol = 2, nrow = 1)

#.................................................#
##### 6 gse73072 Timecourses: auto, mono, CD4 #####
#.................................................#
# Purpose: 
#   - See how Auto-iSEXS, monocytes, and CD4+ T cells change
#     in males and females during influenza infection
#   - Decide whether or not to use the 0.7 cutoff
#
# 0.7 Cutoff
#   - Pros: 
#        - Makes pattern that exists in full data readily visible
#        - Eliminates potential messy garbage samples
#        - Uses same methods as viral challenge paper
#   - Cons: 
#        - Potentially a cherry-picking p-hacking way of getting
#          results that look nice
#        - No proof that low correlation samples are garbage
#   - Conclusions:
#        - 0.7 is okay becausue:
#             - That's what I did for the viral challenge paper
#             - I will validate each plot in WT challenge (supplemental)
#             - Trends that do not validate will be discussed approporiately in Discussion
#      

load('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/6_infection/infection/gse73072/gse73072_flu_ComBat.RData')

# Run immunoStates Deconvolution
gse73072_flu_symptShed$key_comment = 'Created manually'
gse73072_flu_symptShed$exp_comment = 'Created manually'
comboTime_sympt_cells = MetaIntegrator::immunoStatesDecov(list(originalData = list(comboFlu = gse73072_flu_symptShed)))

# Get all data combined for plotting
myData = gse73072_flu_symptShed$pheno
all(rownames(myData) == comboTime_sympt_cells$immunoStates$comboFlu$rn) # T
myData$monocytes = comboTime_sympt_cells$immunoStates$comboFlu$monocyte
myData$CD4 = comboTime_sympt_cells$immunoStates$comboFlu$CD4_positive_alpha_beta_T_cell
myData$correlation = comboTime_sympt_cells$immunoStates$comboFlu$Correlation
myData$autoScore = calculateScore(autoSig, gse73072_flu_symptShed)
myData$Sex = sapply(myData$sex, simpleCap)

# Remove weird extreme point
myData = subset(myData, time_point < 200)

# Make high correlation subset
myData_hiCorr = subset(myData, correlation >= 0.75)

# Set line sizes and point sizes
lineSize = 1.5
pointSize = 1.5

### Auto-iSEXS
fig6b_autoiSEXS = ggplot(myData, aes(x=time_point, y=autoScore, col=Sex)) + 
  geom_smooth(size = lineSize) + geom_point(size=pointSize)+
  theme_cowplot() + ggtitle("Autosomal-iSEXS Score") + 
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('Time (hours)') + ylab('Autosomal-iSEXS Score') +
  theme(legend.position = c(0.05, 0.15))

### Monocytes
ggplot(myData, aes(x=time_point, y=monocytes, col=sex)) + 
  geom_smooth() + geom_point(size=2.5)+
  theme_cowplot() + ggtitle("GSE73072 Influenza Challenge") + 
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('Time (hours)') + ylab('Estimated Monocyte Proportion')

fig6c_mono= ggplot(myData_hiCorr, aes(x=time_point, y=monocytes, col=Sex)) + 
  geom_smooth(size = lineSize) + geom_point(size=pointSize)+
  theme_cowplot() + ggtitle("Monocytes") + 
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('Time (hours)') + ylab('Estimated Monocyte Proportion') +
  theme(legend.position = c(0.75, 0.85))


### CD4s
ggplot(myData, aes(x=time_point, y=CD4, col=sex)) + 
  geom_smooth() + geom_point(size=2.5)+
  theme_cowplot() + ggtitle("GSE73072 Influenza Challenge") + 
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('Time (hours)') + ylab('Estimated CD4+ T cell Proportion')

fig6d_CD4 = ggplot(myData_hiCorr, aes(x=time_point, y=CD4, col=Sex)) + 
  geom_smooth(size=lineSize) + geom_point(size=pointSize)+
  theme_cowplot() + ggtitle("CD4+ T cells") + 
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('Time (hours)') + ylab('Estimated CD4+ T cell Proportion') +
  theme(legend.position = c(0.75, 0.15))



#..............................#
##### Fig 5 - ASSEMBLE!!!  #####
#..............................#


fig5 = plot_grid(fig6b_autoiSEXS, fig6d_CD4,fig6c_mono,
                 nrow=1, ncol=3, labels=LETTERS[c(1:3)])

save_plot(filename = 'plots/fig5_infection.pdf', plot = fig5, ncol = 3, nrow=1)



#....................................#
##### GSE68310 - WT Flu Boxplots #####
#....................................#
# Purpose: 
#   - Look for statistical significance in changes in:
#         - Auto-iSEXS
#         - Estimated CD4's
#         - Monocytes
#
# Background:
#   - Monocytes: Early decrease, potentially more dramatic in females
#   - CD4's: Early decrease, late increase more dramatic in females
#   - Auto: Same as CD4's


# Load from viral challenge project
load('/labs/khatrilab/ebongen/viralChallenge_clean/0_data/wildType/gse68310/gse68310_gem_flu_rhino.RData')

### Flu - perfect separation
quickViolin(gse68310Flu, 'RPS4Y1', 'sex') # perfect separation
quickViolin(gse68310Flu, 'XIST', 'sex') # perfect separation

### fix formatted names
gse68310Flu$formattedName = 'GSE68310 - Influenza'

### Put the timepoint factor in the correct order
# flu
time_point2= factor(gse68310Flu$pheno$time_point, levels = unique(gse68310Flu$pheno$time_point))
all(time_point2 ==gse68310Flu$pheno$time_point) # true
gse68310Flu$pheno$time_point = time_point2

### Make the titer data numeric
# flu - spring titer
head(gse68310Flu$pheno$`spring_anti-a/h1n1_ab_titers`)
springTiter = as.numeric(as.character(gse68310Flu$pheno$`spring_anti-a/h1n1_ab_titers`))
gse68310Flu$pheno$springTiter = springTiter

# flu base titer
baseTiter = as.numeric(as.character(gse68310Flu$pheno$`baseline_anti-a/h1n1_ab_titers`))
table(baseTiter, gse68310Flu$pheno$`baseline_anti-a/h1n1_ab_titers`)
gse68310Flu$pheno$baseTiter = baseTiter


### Limit down to Day 6
phenoShort = subset(gse68310Flu$pheno, !time_point %in% c('Day21', 'Spring'))
unique(phenoShort$time_point)
gse68310Flu_short = subsetGEMFromPheno(gse68310Flu, phenoShort)

### Get data for plotting
myDataFlu = gse68310Flu_short$pheno
myDataFlu$autoScore = calculateScore(autoSig, gse68310Flu_short)
group = paste(myDataFlu$time_point, '\n',myDataFlu$sex, sep='')
myDataFlu$group = factor(group, levels=sort(unique(group)))

### Get p-values

pValuesAuto = c()
for(myTime in unique(myDataFlu$time_point)){
  myTimeData = subset(myDataFlu, time_point == myTime)
  pVal = t.test(myTimeData$autoScore[myTimeData$sex == 'female'],
                myTimeData$autoScore[myTimeData$sex == 'male'])$p.value
  pValuesAuto[myTime] = pVal
}

starsAuto = sapply(pValuesAuto, pToStars)


### Compare sexes across infection
ggplot(myDataFlu, aes(x=group, y=autoScore)) + 
  geom_boxplot()+
  geom_rect(aes(xmin= 2.5, xmax=4.5, ymin=-3, ymax=2.4 ), colour='grey95', fill='grey95')+
  geom_rect(aes(xmin= 6.5, xmax=8.5, ymin=-3, ymax=2.4 ), colour='grey95', fill='grey95')+
  geom_boxplot(aes(col=sex)) + theme_cowplot() +
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Autosomal-iSEXS Score') +
  ggtitle('GSE68310 - Influenza A') +
  geom_signif(annotations = starsAuto['Baseline'], y_position = 2, xmin= 1.01, xmax=1.9, textsize=8) +
  geom_signif(annotations = starsAuto['Day0'], y_position = 2, xmin= 3.01, xmax=3.9, textsize=8) +
  geom_signif(annotations = starsAuto['Day2'], y_position = 2, xmin= 5.01, xmax=5.9, textsize=8) +
  geom_signif(annotations = starsAuto['Day4'], y_position = 2, xmin= 7.01, xmax=7.9, textsize=8) +
  geom_signif(annotations = starsAuto['Day6'], y_position = 2, xmin= 9.01, xmax=9.9, textsize=8)
  
  



### Get data for plotting - separate by sexes first
myDataFlu = gse68310Flu_short$pheno
myDataFlu$autoScore = calculateScore(autoSig, gse68310Flu_short)
group = paste(myDataFlu$time_point, '\n',myDataFlu$sex, sep='')
myDataFlu$group = factor(group, levels=unique(group))

ggplot(myDataFlu, aes(x=group, y=autoScore)) + 
  geom_boxplot(aes(col=sex)) + theme_cowplot() +
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Autosomal-iSEXS Score') +
  ggtitle('GSE68310 - Influenza A') 

#.......................................................#
##### Timecourse Differences in Monocytes and CD4's #####
#.......................................................#
# Purpose: 
#   - GSE73072 implies some sex differences in monocyte and CD4 
#     proportions during influenza
#   - Need to confirm in GSE68310


# Load cell proportion data
load('/labs/khatrilab/ebongen/viralChallenge_clean/0_data/wildType/gse68310/gse68310_cellProp_flu_rhino.RData')

# GEt data to plot
cellFlu = cellPropFlu$pheno
cellFlu$CD4 = unlist(getSampleLevelGeneData(cellPropFlu, "CD4_positive_alpha_beta_T_cell"))
cellFlu$monocyte = unlist(getSampleLevelGeneData(cellPropFlu, 'monocyte'))

# Remove Day 21 and Spring
cellFlu = subset(cellFlu, !time_point %in% c('Day21', 'Spring'))
unique(cellFlu$time_point)

# Create group
group = paste(cellFlu$time_point,'\n', cellFlu$sex, sep = '')
cellFlu$group = factor(group, levels=sort(unique(group)))

# Get p-values
cd4Pvals = c()
monoPvals = c()
for(myTime in unique(cellFlu$time_point)){
  timeData = subset(cellFlu, time_point == myTime)
  pval = t.test(timeData$CD4[timeData$sex == 'female'],
                timeData$CD4[timeData$sex == 'male'])$p.value
  cd4Pvals[myTime] = pval
  
  pval = t.test(timeData$monocyte[timeData$sex == 'female'],
                timeData$monocyte[timeData$sex =='male'])$p.value
  monoPvals[myTime] = pval
}

### CD4+ T cells
ggplot(cellFlu, aes(x=group, y=CD4)) + 
  geom_boxplot()+
  geom_rect(aes(xmin= 2.5, xmax=4.5, ymin=0.12, ymax=0.36 ), colour='grey95', fill='grey95')+
  geom_rect(aes(xmin= 6.5, xmax=8.5, ymin=0.12, ymax=0.36), colour='grey95', fill='grey95')+
  #geom_rect(aes(xmin= 10.5, xmax=12.5, ymin=0.12, ymax=0.36), colour='grey95', fill='grey95')+
  geom_boxplot(aes(col=sex)) + theme_cowplot() +
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Estimated CD4+ T cells') +
  ggtitle('GSE68310 - Influenza A') +
  geom_signif(annotations = 'N.S.', y_position = 0.32, xmin= 1.01, xmax=1.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.32, xmin= 3.01, xmax=3.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.32, xmin= 5.01, xmax=5.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.32, xmin= 7.01, xmax=7.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.32, xmin= 9.01, xmax=9.9, textsize=6)

  
### Monocytes
ggplot(cellFlu, aes(x=group, y=monocyte)) + 
  geom_boxplot()+
  geom_rect(aes(xmin= 2.5, xmax=4.5, ymin=0.15, ymax=0.49 ), colour='grey95', fill='grey95')+
  geom_rect(aes(xmin= 6.5, xmax=8.5, ymin=0.15, ymax=0.49), colour='grey95', fill='grey95')+
  #geom_rect(aes(xmin= 10.5, xmax=12.5, ymin=0.12, ymax=0.36), colour='grey95', fill='grey95')+
  geom_boxplot(aes(col=sex)) + theme_cowplot() +
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Estimated Monocytes') +
  ggtitle('GSE68310 - Influenza A') +
  geom_signif(annotations = 'N.S.', y_position = 0.45, xmin= 1.01, xmax=1.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.45, xmin= 3.01, xmax=3.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.45, xmin= 5.01, xmax=5.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.45, xmin= 7.01, xmax=7.9, textsize=6) +
  geom_signif(annotations = 'N.S.', y_position = 0.45, xmin= 9.01, xmax=9.9, textsize=6)
  
#..................................#
##### Supp Fig - Scores in SLE #####
#..................................#
# Purpose:
#   = One of the reviewers asked what the iSEXS score looks like in SLE
#   = Let's show them how it follows lymphocyte reduction
#   = We'll use gse


# Load SLE dataset with immune cell proportions paired with gene
# gse49454
load("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/5_autoImmunity/SLE/blood/gse49454/gse49454.RData")
load("/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/5_autoImmunity/SLE/blood/gse39088/gse39088.RData")



# Quick check - yup expression is as expected
violinPlot(autoSig, gse49454_baseline, "group2")
violinPlot(xySig, gse49454_baseline, "group2")

# Calcuate proportions
gse49454_baseline$pheno$prop_CD4 = gse49454_baseline$pheno$`cd4 tcells`/gse49454_baseline$pheno$leukocytes
gse49454_baseline$pheno$prop_CD8 = gse49454_baseline$pheno$`cd8 t cells`/gse49454_baseline$pheno$leukocytes
gse49454_baseline$pheno$prop_CD8 = gse49454_baseline$pheno$`cd8 t cells`/gse49454_baseline$pheno$leukocytes
gse49454_baseline$pheno$prop_nk = gse49454_baseline$pheno$`nk cells`/gse49454_baseline$pheno$leukocytes
gse49454_baseline$pheno$prop_neut = gse49454_baseline$pheno$neutrophils/gse49454_baseline$pheno$leukocytes
gse49454_baseline$pheno$prop_lymph = gse49454_baseline$pheno$lymphocytes/gse49454_baseline$pheno$leukocytes

# Add Auto- and XY-iSEXS to pheno
gse39088$pheno$xyISEXS = calculateScore(xySig, gse39088)
gse39088$pheno$autoISEXS = calculateScore(autoSig, gse39088)
gse49454_baseline$pheno$xyISEXS = calculateScore(xySig, gse49454_baseline)
gse49454_baseline$pheno$autoISEXS = calculateScore(autoSig, gse49454_baseline)

# Create a version that's just SLE
phenoSLE = subset(gse49454_baseline$pheno, group == "SLE")
gse49454_baseline_sle = subsetGEMFromPheno(gse49454_baseline, phenoSLE)

# Initiate lists for storing the plots of the SLE figures
# one list for the top row, and one list for the bottom
sleFigListTop = list()
sleFigListBottom = list()

# Constants for the slePlots
beeswarmPointSize = 0.7
signifFontSize = 3

### Box plot - gse39088 - XY-iSEXS
sleFigListTop$gse3_xy = ggplot(gse39088$pheno, aes(x=group2, y=xyISEXS)) + geom_boxplot()+ 
  geom_beeswarm(aes(col=sex), cex=1.5, priority = 'random', size=beeswarmPointSize) +
  scale_color_manual(values = c(femColor, malColor)) +theme_cowplot() +
  theme(legend.position = 'none') +xlab('') + ylab('XY-iSEXS Score') +
  ggtitle("GSE39088\nXY-iSEXS") + 
  geom_signif(xmin=1.2, xmax=1.8, annotations = 'p=6.70E-4', y_position = 2.1, textsize = signifFontSize)+
  ylim(c(min(gse39088$pheno$xyISEXS), 2.5))

t.test(gse39088$pheno$xyISEXS[gse39088$pheno$group == "SLE"],
       gse39088$pheno$xyISEXS[gse39088$pheno$group == 'healthy'])


### Box plot - gse39088 - Auto-iSEXS
sleFigListTop$gse3_auto = ggplot(gse39088$pheno, aes(x=group2, y=autoISEXS)) + geom_boxplot()+ 
  geom_beeswarm(aes(col=sex), cex=1.5, priority = 'random', size=0.8) +
  scale_color_manual(values = c(femColor, malColor)) +theme_cowplot() +
  theme(legend.position = 'none') +xlab('') + ylab('Autosomal-iSEXS Score') +
  ggtitle("GSE39088\nAutosomal-iSEXS") + 
  geom_signif(xmin=1.2, xmax=1.8, annotations = 'p=9.19E-10', y_position = 2.5, textsize = signifFontSize)+
  ylim(c(min(gse39088$pheno$autoISEXS), 3))


t.test(gse39088$pheno$autoISEXS[gse39088$pheno$group == "SLE"],
       gse39088$pheno$autoISEXS[gse39088$pheno$group == 'healthy'])

### Box plots - gse49454_baseline - XY-iSEXS
sleFigListTop$gse4_xy = ggplot(gse49454_baseline$pheno, aes(x=group2, y=xyISEXS)) + geom_boxplot()+ 
  geom_beeswarm(aes(col=sex), cex=1.5, priority = 'random', size=beeswarmPointSize, alpha=0.7) +
  scale_color_manual(values = c(femColor, malColor)) +theme_cowplot() +
  theme(legend.position = 'none') +xlab('') + ylab('XY-iSEXS Score') +
  ggtitle("GSE49454  XY-iSEXS") +  
  geom_signif(xmin=1.1, xmax=1.8, annotations = 'p=0.0014', y_position = 1.1, textsize = signifFontSize) +
  geom_signif(xmin=3.1, xmax=3.8, annotations = 'N.S.', y_position = -1.5, textsize = signifFontSize) +
  scale_x_discrete(labels= c("Female\nHealthy", "Female\nSLE", "Male\nHealthy", "Male\nSLE")) +
  ylim(c(-2.5, 1.3))


t.test(gse49454_baseline$pheno$xyISEXS[gse49454_baseline$pheno$group2 == "Female Healthy"],
       gse49454_baseline$pheno$xyISEXS[gse49454_baseline$pheno$group2 == "Female SLE"]) #P= 0.0014
t.test(gse49454_baseline$pheno$xyISEXS[gse49454_baseline$pheno$group2 == "Male Healthy"],
       gse49454_baseline$pheno$xyISEXS[gse49454_baseline$pheno$group2 == "Male SLE"]) # NS


### Box Plot - gse49454 - Auto-iSEXS
sleFigListTop$gse4_auto = ggplot(gse49454_baseline$pheno, aes(x=group2, y=autoISEXS)) + geom_boxplot()+ 
  geom_beeswarm(aes(col=sex), cex=1.5, priority = 'random', size=beeswarmPointSize, alpha=0.7) +
  scale_color_manual(values = c(femColor, malColor)) +theme_cowplot() +
  theme(legend.position = 'none') +xlab('') + ylab('Autosomal-iSEXS Score') +
  ggtitle("GSE49454  Autosomal-iSEXS") +
  geom_signif(xmin=1.1, xmax=1.8, annotations = 'p=9.41E-5', y_position = 2.2, textsize = signifFontSize) +
  geom_signif(xmin=3.1, xmax=3.8, annotations = 'p=0.10', y_position = 1.3, textsize = signifFontSize) +
  scale_x_discrete(labels= c("Female\nHealthy", "Female\nSLE", "Male\nHealthy", "Male\nSLE")) +
  ylim(c(-2.1, 2.4))

t.test(gse49454_baseline$pheno$autoISEXS[gse49454_baseline$pheno$group2 == "Female Healthy"],
       gse49454_baseline$pheno$autoISEXS[gse49454_baseline$pheno$group2 == "Female SLE"]) #P= 9.41E-5
t.test(gse49454_baseline$pheno$autoISEXS[gse49454_baseline$pheno$group2 == "Male Healthy"],
       gse49454_baseline$pheno$autoISEXS[gse49454_baseline$pheno$group2 == "Male SLE"]) # p=0.10




### Scatter plot - gse49454 - SLEDAI vs XY
sleFigListBottom$gse4_xySLEDAI = ggplot(gse49454_baseline_sle$pheno, aes(x=sledai, y=xyISEXS, col=group2)) + 
  geom_smooth(method="lm") + geom_point(size=0.75) + 
  scale_color_manual(values = c(femColor, malColor)) + theme_cowplot() +
  xlab("SLEDAI") + ylab("XY-iSEXS Score") +labs(col = "") +
  theme(legend.position = "none") + ggtitle("GSE49454") +
  draw_label("p=N.S.", x=12, y=-0.5, colour=femColor, fontface="bold", hjust = 0, size=10)+
  draw_label("p=N.S.", x=12, y=-1.7, colour=malColor, fontface="bold", hjust = 0, size=10)



cor.test(gse49454_baseline_sle$pheno$sledai[gse49454_baseline_sle$pheno$sex == "female"],
         gse49454_baseline_sle$pheno$xyISEXS[gse49454_baseline_sle$pheno$sex == "female"]) # NS
cor.test(gse49454_baseline_sle$pheno$sledai[gse49454_baseline_sle$pheno$sex == "male"],
         gse49454_baseline_sle$pheno$xyISEXS[gse49454_baseline_sle$pheno$sex == "male"]) # N.S.

### Scatter plot - gse49454 - SLEDAI vs Auto

sleFigListBottom$gse4_autoSLEDAI = ggplot(gse49454_baseline_sle$pheno, aes(x=sledai, y=autoISEXS, col=group2)) + 
  geom_smooth(method="lm") + geom_point(size=0.75) + 
  scale_color_manual(values = c(femColor, malColor)) + theme_cowplot() +
  xlab("SLEDAI") + ylab("Autosomal-iSEXS Score") +labs(col = "") +
  theme(legend.position = 'none') + ggtitle("GSE49454") +
  draw_label("p=N.S.", x=23, y=1.2, colour=femColor, fontface="bold", size=10)+
  draw_label("p=N.S.", x=23, y=0.8, colour=malColor, fontface="bold", size=10)
  
cor.test(gse49454_baseline_sle$pheno$sledai[gse49454_baseline_sle$pheno$sex == "female"],
         gse49454_baseline_sle$pheno$autoISEXS[gse49454_baseline_sle$pheno$sex == "female"]) # NS
cor.test(gse49454_baseline_sle$pheno$sledai[gse49454_baseline_sle$pheno$sex == "male"],
         gse49454_baseline_sle$pheno$autoISEXS[gse49454_baseline_sle$pheno$sex == "male"]) # NS




### Scatter plot - gse49454 - XY vs % of CD4's
sleFigListBottom$gse4_cd4_xy = ggplot(gse49454_baseline_sle$pheno, aes(x=prop_CD4*100, y=xyISEXS, col=group2)) + 
  geom_smooth(method="lm") + geom_point(size = 0.75) + 
  scale_color_manual(values = c(femColor, malColor)) + theme_cowplot() +
  xlab("CD4+ T cells\n(% of Leukocytes)") + ylab("XY-iSEXS Score") +labs(col = "") +
  theme(legend.position = "none") + ggtitle("GSE49454") +
  draw_label("r=0.60,\np=5.52E-5", x=15, y=-0.4, colour=femColor, fontface="bold", hjust = 0, size=10)+
  draw_label("r=0.52,\np=0.22", x=10, y=-1.7, colour=malColor, fontface="bold", hjust = 0, size=10)


# XY-iSEXS score also correlates with CD4 T cells
cor.test(gse49454_baseline_sle$pheno$prop_CD4[gse49454_baseline_sle$pheno$sex == "female"],
         gse49454_baseline_sle$pheno$xyISEXS[gse49454_baseline_sle$pheno$sex == "female"]) # r= 0.60, p=5.521e-05
cor.test(gse49454_baseline_sle$pheno$prop_CD4[gse49454_baseline_sle$pheno$sex == "male"],
         gse49454_baseline_sle$pheno$xyISEXS[gse49454_baseline_sle$pheno$sex == "male"]) # r=0.52, p=0.22


### Scatter plot- gse49454 - Auto vs% of CD4's

sleFigListBottom$gse4_cd4_auto = ggplot(gse49454_baseline_sle$pheno, aes(x=prop_CD4*100, y=autoISEXS, col=group2)) + 
  geom_smooth(method="lm") + geom_point(size=0.75) + 
  scale_color_manual(values = c(femColor, malColor)) + theme_cowplot() +
  xlab("CD4+ T cells\n(% of Leukocytes)") + ylab("Autosomal-iSEXS Score") +
  labs(col = "") + ggtitle("GSE49454") +
  draw_label("r=0.35,\np=0.03", x=17, y=-0.5, colour=femColor, fontface="bold", hjust = 0, size=10)+
  draw_label("N.S.", x=12, y=-1.5, colour=malColor, fontface="bold", hjust = 0, size=10)

cor.test(gse49454_baseline_sle$pheno$prop_CD4[gse49454_baseline_sle$pheno$sex == "female"],
         gse49454_baseline_sle$pheno$autoISEXS[gse49454_baseline_sle$pheno$sex == "female"]) # r=0.35, p=0.03
cor.test(gse49454_baseline_sle$pheno$prop_CD4[gse49454_baseline_sle$pheno$sex == "male"],
         gse49454_baseline_sle$pheno$autoISEXS[gse49454_baseline_sle$pheno$sex == "male"]) # Not significant







### Combine in to suppFig
sleFigTop = plot_grid(plotlist = sleFigListTop, nrow=1, ncol=4, 
                      rel_widths = c(0.6, 0.6, 1, 1), labels = LETTERS[1:4])
sleFigBottom = plot_grid(plotlist = sleFigListBottom, nrow=1, ncol=4, 
                         rel_widths = c(1,1,1,1.7),labels = LETTERS[5:8])

sleFigFull = plot_grid(sleFigTop, sleFigBottom, nrow=2, ncol=1)

save_plot("plots/sleFigure.pdf", sleFigFull, ncol = 2, nrow = 2, base_height = 3, base_width = 6)


#.......................#
##### XY- Supp Figs #####
#.......................#
# Purpose:
#  = One of the reviewers wanted to see what diff plots look like with XY-iSEXS
#  = Now, I think splitting them into X-iSEXS and Y-iSEXS are most useful
#  = So, I'll try both!


# Create lists of figures:
xyiSEXS_figList = list()

# Constants for the figures
fontSize = 8
pointSize = 1
lineSize= 1

### Correlate with monocytes

# Get cell type columns
myColumns6 = colnames(gse65133$pheno)[14:22]
myColumns4 = colnames(gse47353_baseline$pheno[30:58])

### Get data for gse65133 - CIBERSORT paper data
myData6 = gse65133$pheno
myData6$scoreXY = calculateScore(sexMetaObj$filterResults$xySig, gse65133)

# Make sex be capitalized
myData6$sex = sapply(myData6$sex, simpleCap)

# Calculate stats - Monocytes
cor.test(myData6$Monos[myData6$sex=='Male'], myData6$scoreXY[myData6$sex == "Male"]) # N.S.
cor.test(myData6$Monos[myData6$sex=='Female'], myData6$scoreXY[myData6$sex == "Female"]) # N.S.

xyiSEXS_figList$a_monoCor6 = ggplot(myData6, aes(x=Monos, y=scoreXY, col=sex)) + 
  geom_smooth(method='lm') +
  geom_point(size=pointSize)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes\n(% of Mononuclear Cells)") +
  ylab("XY-iSEXS Score") +
  draw_label(label = "N.S.", x= 30, 0.7, col=femColor, size = fontSize) +
  draw_label(label = "N.S.", x= 35, -0.5, col=malColor, size=fontSize) +
  ggtitle("GSE65133 - PBMCs") +
  theme_cowplot(font_size = fontSize)+
  theme(legend.position="none")+
  labs(color='')

# Create CD4+ T cell column, right now split into subsets
CD4_Tcells = c("CD4.naïve.T.cells", "CD4.memory.T.cells.resting", "CD4.memory.T.cells.activated")
myData6$CD4_Tcells = apply(myData6[,CD4_Tcells], 1, sum)

# Calculate stats - CD4 T cells
cor.test(myData6$CD4_Tcells[myData6$sex=='Male'], myData6$scoreXY[myData6$sex == "Male"]) # N.S.
cor.test(myData6$CD4_Tcells[myData6$sex=='Female'], myData6$scoreXY[myData6$sex == "Female"]) # N.S.


xyiSEXS_figList$bb_Tcor6 = ggplot(myData6, aes(x=CD4_Tcells, y=scoreXY, col=sex)) + 
  geom_smooth(method='lm') +
  geom_point(size=pointSize)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("CD4+ T cells\n(% of Mononuclear Cells)") +
  ylab("XY-iSEXS Score") +
  draw_label(label = "N.S.", x= 20, y=1, col=femColor, size=fontSize) +
  draw_label(label = "N.S.", x= 50, -0.5, col=malColor, size = fontSize) +
  ggtitle("GSE65133 - PBMCs") +
  theme_cowplot(font_size = fontSize)+
  theme(legend.position="none")+
  labs(color='')





### b) Plot monocyte correlation in gse4
myData4 = gse47353_baseline$pheno
myData4$sex = sapply(myData4$sex, function(x) simpleCap(x))
myData4$scoreXY = calculateScore(sexMetaObj$filterResults$xySig,gse47353_baseline)


# Get stats
cor.test(myData4$Total_Monocytes[myData4$sex=='Male'], myData4$scoreXY[myData4$sex == "Male"]) # p=0.02
cor.test(myData4$Total_Monocytes[myData4$sex=='Female'], myData4$scoreXY[myData4$sex == "Female"]) # p=0.06


xyiSEXS_figList$b_monoCor4 = ggplot(myData4, aes(x=Total_Monocytes, y=scoreXY, col=sex)) + 
  geom_smooth(method='lm') +
  geom_point(size=pointSize)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes\n(% of CD45+ Cells)") +
  ylab("XY-iSEXS Score") +
  draw_label(label = "r=-0.50, p=0.02", x= 28, -0.9, col=malColor, size = fontSize) +
  draw_label(label = "r=-0.34, p=0.06", x= 28, 0.25, col= femColor, size = fontSize) +
  ggtitle('GSE47353 - PBMCs') +
  theme_cowplot(font_size = fontSize) +
  theme(legend.position = "none") +
  labs(color ='')


xyiSEXS_figList$bb_Tcor4 = ggplot(myData4, aes(x=Total_T_cells, y=scoreXY, col=sex)) + 
  geom_smooth(method='lm') +
  geom_point(size=pointSize)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("T cells\n(% of CD45+ Cells)") +
  ylab("XY-iSEXS Score") +
  draw_label(label = "N.S.", x= 50, -1, col=malColor, size=fontSize) +
  draw_label(label = "N.S.", x= 50, 0.6, col= femColor, size=fontSize) +
  ggtitle('GSE47353 - PBMCs') +
  theme_cowplot(font_size = fontSize) +
  theme(legend.position = "none") +
  labs(color ='')


cor.test(myData4$Total_T_cells[myData4$sex == "Female"],
         myData4$scoreXY[myData4$sex == "Female"])# not sig
cor.test(myData4$Total_T_cells[myData4$sex == "Male"],
         myData4$scoreXY[myData4$sex == "Male"]) # Not sig





### XY-iSEXS with age
xyiSEXS_figList$c_xyTimecourse = timecourseXY


### XY-iSEXS with flu infection
# Get all data combined for plotting
myData = gse73072_flu_symptShed$pheno
all(rownames(myData) == comboTime_sympt_cells$immunoStates$comboFlu$rn) # T
myData$correlation = comboTime_sympt_cells$immunoStates$comboFlu$Correlation
myData$xyScore = calculateScore(xySig, gse73072_flu_symptShed)
myData$Sex = sapply(myData$sex, simpleCap)

# Remove weird extreme point
myData = subset(myData, time_point < 200)

# Make high correlation subset
myData_hiCorr = subset(myData, correlation >= 0.75)




### XY-iSEXS
xyiSEXS_figList$d_aging = ggplot(myData, aes(x=time_point, y=xyScore, col=Sex)) + 
  geom_smooth(size = lineSize) + geom_point(size=pointSize)+
  theme_cowplot(font_size = fontSize) + ggtitle("GSE73072 - Flu Infection") + 
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('Time (hours)') + ylab('XY-iSEXS Score') +
  theme(legend.position = "none")


### XY-iSEXS and anti-flu titers
# Get data
myDataFlu = gse68310Flu$pheno
myDataFlu$xyScore = calculateScore(xySig, gse68310Flu)
myDataFlu$Sex = sapply(myDataFlu$sex, simpleCap)

# Get baseline baseline data
myDataFluBase = subset(myDataFlu, time_point %in% c('Baseline'))
table(myDataFluBase$sex, myDataFluBase$time_point) # Males and females, only baseline

# Make xlab with delta
delta_xlab = expression(paste(Delta, 'Titers (anti-H1N1 Ab)'))

# Get the stats
myDataFluBase$deltaTiters = myDataFluBase$springTiter - myDataFluBase$baseTiter
cor.test(myDataFluBase$deltaTiters[myDataFluBase$Sex == "Male"],
         myDataFluBase$xyScore[myDataFluBase$Sex == "Male"]) # Not sig
cor.test(myDataFluBase$deltaTiters[myDataFluBase$Sex == "Female"],
         myDataFluBase$xyScore[myDataFluBase$Sex == "Female"]) # Not sig

xyiSEXS_figList$e_deltaTiters= ggplot(myDataFluBase, aes(x=springTiter-baseTiter, y=xyScore, col=Sex)) +
  geom_smooth(method = 'lm', size=lineSize)+
  geom_point(size=pointSize) + theme_cowplot(font_size = fontSize) + ggtitle("GSE68310")+
  scale_color_manual(values = c(femColor, malColor)) +
  ylab('XY-iSEXS at Baseline') + xlab(delta_xlab) +
  draw_label('N.S.', x=6, y=-0.7, colour=malColor, size=fontSize) +
  draw_label('N.S.', x=5, y=0.5, colour = femColor, size=fontSize) +
  theme(legend.position = "none") +labs(col='')

### is the XY-iSEXS score consistent?

# Get Autosomal score form gse68310
fluBaseSpring = gse68310Flu$pheno
fluBaseSpring$xyScore = calculateScore(xySig, gse68310Flu)

# Subset to only baseline and spring
unique(fluBaseSpring$time_point)
fluBaseSpring = subset(fluBaseSpring, time_point %in% c('Baseline', 'Spring'))

# dcast to make it autosomal scores at baseline and spring
fluBaseSpring = data.table::dcast(fluBaseSpring, subject_id + sex ~ time_point, value.var = 'xyScore')

# Subset to males and females 
fluBaseSpring_fem = subset(fluBaseSpring, sex == 'female')
fluBaseSpring_mal = subset(fluBaseSpring, sex == 'male')

# Get stats for males and females
malTest = cor.test(fluBaseSpring_mal$Baseline, fluBaseSpring_mal$Spring)
femTest = cor.test(fluBaseSpring_fem$Baseline, fluBaseSpring_fem$Spring)

# Extract stats to make captions
malLabel = paste('r=', signif(malTest$estimate,2), '  p=', signif(malTest$p.value,2), sep='')
femLabel = paste('r=', signif(femTest$estimate,2), '  p=', signif(femTest$p.value,2), sep='')

# Plot of females in spring vs Baseline
xyiSEXS_figList$f_femYearCor = ggplot(fluBaseSpring_fem, aes(x=Baseline, y=Spring)) +
  geom_smooth(method = 'lm', colour=femColor, size=lineSize)+
  geom_point(size=pointSize, colour=femColor) +
  theme_cowplot(font_size=fontSize) + 
  xlab('Baseline XY-iSEXS Score') +
  ylab(paste('Post-Flu Season', '\n', 'XY-iSEXS Score')) + 
  ggtitle('GSE68310 Females') +
  draw_label(femLabel, x=0.7,y=1.2, colour=femColor, size=fontSize)


xyiSEXS_figList$g_malYearCor = ggplot(fluBaseSpring_mal, aes(x=Baseline, y=Spring)) +
  geom_smooth(method = 'lm', colour=malColor, size=lineSize)+
  geom_point(size=pointSize, colour=malColor) +
  theme_cowplot(font_size = fontSize) + 
  xlab('Baseline XY-iSEXS Score') +
  ylab(paste('Post-Flu Season', '\n', 'XY-iSEXS Score')) + 
  ggtitle('GSE68310 Males')+
  draw_label(malLabel, x=-0.9,y=-0.9, colour=malColor, size=fontSize)


### Assemble into a grid

# Top row: correlation plots with monocytes and CD4's
xyPlot_top = plot_grid(xyiSEXS_figList$b_monoCor4, xyiSEXS_figList$a_monoCor6, 
                       #xyiSEXS_figList$bb_Tcor4, xyiSEXS_figList$bb_Tcor6,
                       nrow=1, ncol=2, labels = LETTERS[1:2])
xyPlot_middle = plot_grid(xyiSEXS_figList$c_xyTimecourse, xyiSEXS_figList$d_aging, nrow=1, ncol=2, labels=c(LETTERS[3:4]))

xyPlot_bottom = plot_grid(xyiSEXS_figList$e_deltaTiters, xyiSEXS_figList$f_femYearCor, xyiSEXS_figList$g_malYearCor,
                          nrow=1, ncol=3, labels = LETTERS[5:8])

xyPlot = plot_grid(xyPlot_top, xyPlot_middle, xyPlot_bottom, 
                   nrow=3, ncol=1)

save_plot("plots/xyiSEXS_plot.pdf", xyPlot, ncol=2, nrow=3, base_height = 2, base_width = 3.4)

#........................................#
##### Fig s3: Monocyte Rearrangement #####
#........................................#
# Purpose:
#   = Purvesh had some good suggestsions for tweaking
#     the supplemental figure order, so I'm trying to make
#     a prettier version of the figure he made
# 
# Subfigures:
#   a) Gabi's cohort violin plot - sFig2_list$a_gabi
#   b) SDY 212 violin plot - sFig2_list$b_sdy212
#   c) gse6 monocyte correlation - Autosomal iSEXS
#   d) gse4 monocyte correlation - Autosomal iSEXS
#   e) gse6 monocyte correlation - XY-iSEXS
#   f) gse4 monocyte correlation -XY-isexs
#   g) XY-iSEXS and aging
#
# Order:
#   a b c d
#   e f  g


sFig3_mono = list()
fontSize = 10
pointSize = 1
lineSize = 1

### a - Gabi 
# Get stats
youngPval = t.test(gabi$monocyte[gabi$group == 'Younger Female'],
                   gabi$monocyte[gabi$group == 'Younger Male'])$p.value
oldPval = t.test(gabi$monocyte[gabi$group == 'Older Female'],
                 gabi$monocyte[gabi$group == 'Older Male'])$p.value

sFig3_mono$a_gabi = ggplot(gabi, aes(x=group2, y=monocyte*100)) + 
  geom_violin(trim=F, fill='grey') + geom_jitter(width=0.1, size= pointSize, aes(col=sex))+
  scale_color_manual(values = c(femColor, malColor)) + ylab('Monocytes (% of Mononuclear Cells)') +
  xlab('') + theme_cowplot(font_size = fontSize) + theme(legend.position = 'none') + 
  geom_signif(annotations = paste('p =', signif(youngPval,2)), y_position = 40, xmin= 1.2, xmax=1.8, textsize=2.75)+
  geom_signif(annotations = paste('p =', signif(oldPval,2)), y_position = 30, xmin= 3.2, xmax=3.8, textsize=2.75) +
  ggtitle('CyTOF - Fragiadakis & Bjornson') 



### b - SDY212
# Get p-values
youngPval = t.test(sdy212$pheno$autoScore[sdy212$pheno$group == 'Younger Female'],
                   sdy212$pheno$autoScore[sdy212$pheno$group == 'Younger Male'])$p.value
oldPval = t.test(sdy212$pheno$autoScore[sdy212$pheno$group == 'Older Female'],
                 sdy212$pheno$autoScore[sdy212$pheno$group == 'Older Male'])$p.value

sFig3_mono$b_sdy212 = ggplot(sdy212$pheno, aes(x=group2, y=autoScore)) + 
  geom_violin(fill='grey', trim=F) + geom_jitter(width=0.1, size=pointSize, aes(col=sex)) +
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Autosomal-iSEXS Score') + theme(legend.position = 'none')+
  theme_cowplot(font_size = fontSize) + 
  geom_signif(annotations = paste('p =', signif(youngPval,2)), y_position = 3.7, xmin= 1.2, xmax=1.8, textsize=2.75)+
  geom_signif(annotations = paste('p =', signif(oldPval,2)), y_position = 2.2, xmin= 3.2, xmax=3.8, textsize=2.75) +
  theme(legend.position = 'none') + ggtitle('SDY212 - Whole Blood')


### c - gse6 monocyte correlation
# Get cell type columns
myColumns6 = colnames(gse65133$pheno)[14:22]
myColumns4 = colnames(gse47353_baseline$pheno[30:58])

### Get data for gse65133 - CIBERSORT paper data
myData6 = gse65133$pheno
myData6$scoreXY = calculateScore(sexMetaObj$filterResults$xySig, gse65133)
myData6$scoreAuto = calculateScore(sexMetaObj$filterResults$autoSig, gse65133)

# Make sex be capitalized
myData6$sex = sapply(myData6$sex, simpleCap)

# Calculate stats
cor.test(myData6$Monos, myData6$scoreAuto) # r=-0.79, p=2.9E-5

sFig3_mono$c_monoCor6auto = ggplot(myData6, aes(x=Monos, y=scoreAuto)) + 
  geom_smooth(method='lm', col='black') +
  geom_point(size=pointSize, aes(col=sex))+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes\n(% of Mononuclear Cells)") +
  ylab("XY-iSEXS Score") +
  draw_label(label = "r=-0.79,\np=2.9E-5", x= 30, 0.7, size = fontSize) +
  ggtitle("GSE65133 - PBMCs") +
  theme_cowplot(font_size = fontSize)+
  theme(legend.position="none")+
  labs(color='')


# Calculate stats - Monocytes
cor.test(myData6$Monos[myData6$sex=='Male'], myData6$scoreXY[myData6$sex == "Male"]) # N.S.
cor.test(myData6$Monos[myData6$sex=='Female'], myData6$scoreXY[myData6$sex == "Female"]) # N.S.

sFig3_mono$e_monoCor6xy = ggplot(myData6, aes(x=Monos, y=scoreXY, col=sex)) + 
  geom_smooth(method='lm') +
  geom_point(size=pointSize)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes\n(% of Mononuclear Cells)") +
  ylab("XY-iSEXS Score") +
  draw_label(label = "N.S.", x= 30, 0.7, col=femColor, size = fontSize) +
  draw_label(label = "N.S.", x= 35, -0.5, col=malColor, size=fontSize) +
  ggtitle("GSE65133 - PBMCs") +
  theme_cowplot(font_size = fontSize)+
  theme(legend.position="none")+
  labs(color='')



### d - Plot monocyte vs Auto-iSEXS in gse4
myData4 = gse47353_baseline$pheno
myData4$sex = sapply(myData4$sex, function(x) simpleCap(x))
myData4$scoreXY = calculateScore(sexMetaObj$filterResults$xySig,gse47353_baseline)
myData4$scoreAuto = calculateScore(sexMetaObj$filterResults$autoSig,gse47353_baseline)



# Get stats
cor.test(myData4$Total_Monocytes, myData4$scoreAuto) # r=-0.59, p=4.4E-6

sFig3_mono$d_monoCor4auto = ggplot(myData4, aes(x=Total_Monocytes, y=scoreAuto)) + 
  geom_smooth(method='lm', col='black') +
  geom_point(size=pointSize, aes(col=sex))+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes\n(% of CD45+ Cells)") +
  ylab("Autosomal-iSEXS Score") +
  draw_label(label = "r=-0.59,\np=4.4E-6", x= 35, 0.9, size = fontSize) +
  ggtitle('GSE47353 - PBMCs') +
  theme_cowplot(font_size = fontSize) +
  theme(legend.position = "none") +
  labs(color ='')


# Get stats
cor.test(myData4$Total_Monocytes[myData4$sex=='Male'], myData4$scoreXY[myData4$sex == "Male"]) # p=0.02
cor.test(myData4$Total_Monocytes[myData4$sex=='Female'], myData4$scoreXY[myData4$sex == "Female"]) # p=0.06


sFig3_mono$f_monoCor4xy = ggplot(myData4, aes(x=Total_Monocytes, y=scoreXY, col=sex)) + 
  geom_smooth(method='lm') +
  geom_point(size=pointSize)+scale_color_manual(values = c(femColor,malColor)) + 
  xlab("Monocytes\n(% of CD45+ Cells)") +
  ylab("XY-iSEXS Score") +
  draw_label(label = "r=-0.50, p=0.02", x= 28, -0.9, col=malColor, size = fontSize) +
  draw_label(label = "r=-0.34, p=0.06", x= 28, 0.25, col= femColor, size = fontSize) +
  ggtitle('GSE47353 - PBMCs') +
  theme_cowplot(font_size = fontSize) +
  theme(legend.position = "none") +
  labs(color ='')

# Aging - xyiSEXS

sFig3_mono$g_age = timecourseXY




topRow = plot_grid(sFig3_mono$a_gabi, sFig3_mono$b_sdy212, sFig3_mono$c_monoCor6auto, 
                   sFig3_mono$d_monoCor4auto, nrow = 1, ncol = 4, labels= LETTERS[1:4])

bottomRow = plot_grid(sFig3_mono$e_monoCor6xy, sFig3_mono$f_monoCor4xy, sFig3_mono$g_age, 
                      rel_widths = c(1,1,2), nrow = 1, ncol=3, labels = LETTERS[5:7])


fullPlot = plot_grid(topRow, bottomRow, nrow = 2, ncol=1)

save_plot("plots/sFig3_monocytev2.0.pdf", fullPlot, base_width = 10, base_height = 6)


#...........................................#
##### Rearranged sFig4 - Flu Timecourse #####
#...........................................#

### Flu infection
# Get all data combined for plotting
myData = gse73072_flu_symptShed$pheno
all(rownames(myData) == comboTime_sympt_cells$immunoStates$comboFlu$rn) # T
myData$correlation = comboTime_sympt_cells$immunoStates$comboFlu$Correlation
myData$xyScore = calculateScore(xySig, gse73072_flu_symptShed)
myData$Sex = sapply(myData$sex, simpleCap)

# Remove weird extreme point
myData = subset(myData, time_point < 200)

# Make high correlation subset
myData_hiCorr = subset(myData, correlation >= 0.75)


### Auto-iSEXS
sFig4_xyFlu = ggplot(myData, aes(x=time_point, y=xyScore, col=Sex)) + 
  geom_smooth(size = 1.5) + geom_point(size=1.5)+
  theme_cowplot(font_size = 14) + ggtitle("GSE73072 - Flu Infection") + 
  scale_color_manual(values = c(femColor, malColor)) +
  xlab('Time (hours)') + ylab('XY-iSEXS Score') +
  theme(legend.position = "none")

save_plot("plots/sFig4_xyFlu.pdf", sFig4_xyFlu, base_height = 4, base_width = 6)
