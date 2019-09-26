# May 11, 2018 - Download and preprocess flow cytometry data from 
# Patin et al. 2018 
# Natural variation in the parameters of innate immune cells is preferentially driven by genetic factors
# Nature Immunology 

#....................#
##### Load Stuff #####
#....................#
# Load libraries
library(ggplot2)
library(cowplot)

# Load funx code
source("~/Tools/Graphing Scripts/quickPlots.R")

setwd('/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/4_sortedCells/cellProp/patin/')

#..............................#
##### Download the Package #####
#..............................#
# Purpose: 
#   - Patin made their flow data available as a package on github
#   - package is called: mmi
#   - Download here using devtools

library(devtools)

devtools::install_github("JacobBergstedt/mmi")

library(mmi)

#...........................#
##### Explore Flow Data #####
#...........................#


dim(facs) # flow data
dim(ecrf)# demographic data

# Explore facs data structure
is.data.frame(facs)
View(facs[1:5, 1:5])

# Get cell type names
colnames(facs)
ncol(facs) -1 # 166 cell types
cellTypes = colnames(facs)[2:ncol(facs)]

# Look at cell types including 'mono'
# Probably most useful to look at MFI or N of monocytes
monoVariables = cellTypes[grepl(pattern = 'mono', x = cellTypes, ignore.case = T)]

# Get number of leukocytes
monoVariables = c(cellTypes[grepl(pattern = "CD45", x = cellTypes)], monoVariables)

#..................................#
##### Explore demographic data #####
#..................................#

#44 variables
dim(ecrf)

# Really interesting variables!
colnames(ecrf)

# 399 males, 417 females
table(ecrf$Sex)

demoVariables = c("Age", "Sex")


#................................#
##### Combine and graph data #####
#................................#

# Create data frame of data to plot
myData = facs[,c('SUBJID', monoVariables)]
all(myData$SUBJID == ecrf$SUBJID) # true
myData$sex = ecrf$Sex
myData$age = ecrf$Age


# Examine data
View(head(myData))
table(myData$sex) # 399 males, 417 females
hist(myData$age) # good coverage of ages 20-70

#....................#
##### Raw Values #####
#....................#
### N_monocytes
# Consistently higher in males, but hard to say significant
ggplot(myData, aes(x=age, y=N_mono.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor))

### N CD14hi Monocytes
# slightly more separation 
ggplot(myData, aes(x=age, y=N_CD14hi_mono.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor))

### N CD16hi Monocytes
# entertwined twisty
ggplot(myData, aes(x=age, y=N_CD16hi_mono.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor))

### MFI CD16 in CD16hi monocytes
# trend mroe in men, then entertwined after 50
ggplot(myData, aes(x=age, y=MFI_CD16_in_CD16hi_mono.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor))

### MFI CD16 in CD14hi monocytes
# entertwined
ggplot(myData, aes(x=age, y=MFI_CD16_in_CD14hi_mono.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor))


#.........................#
##### Relative Values #####
#.........................#

###  Total monocytes
#Relative to CD45 huge difference!
ggplot(myData, aes(x=age, y=N_mono.panel5/N_CD45pos.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor)) +theme_cowplot()

### CD14+ monocytes
# Similar, slighlty more squished difference, I say report total monocytes
ggplot(myData, aes(x=age, y=N_CD14hi_mono.panel5/N_CD45pos.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor))

### CD16+ monocytes
# weak difference and all entertwiny
ggplot(myData, aes(x=age, y=N_CD16hi_mono.panel5/N_CD45pos.panel5, col=sex)) + geom_smooth()+
  geom_point() +scale_color_manual(values = c(malColor, femColor))

#...........................#
##### Save demographics #####
#...........................#

write.table(ecrf, file = "demographics.txt", sep='\t', quote=F)

