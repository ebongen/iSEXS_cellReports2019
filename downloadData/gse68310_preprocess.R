# July 3, 2018 - Prep GSE68310 for sex differences work

# Purpose: 
#   - I already have GSE68310 from viral challenge project
#   - But, there's some preprocessing that isn't saved within the original object
#   - It'll be easier to use if I do that preprocessing here

#....................#
##### Load Stuff #####
#....................#
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

# Load dataset 
load('/labs/khatrilab/ebongen/viralChallenge_clean/0_data/wildType/gse68310/gse68310_gem_flu_rhino.RData')

# Source useful code
source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

setwd('0_datasets/6_infection/infection/gse68310_WT_flu/')

#..........................#
##### Check sex labels #####
#..........................#

# Check sex labels
'RPS4Y1' %in% gse68310Flu$keys # TRUE
'KDM5D' %in% gse68310Flu$keys # FALSE
'XIST' %in% gse68310Flu$keys # TRUE

### Flu - perfect separation
quickViolin(gse68310Flu, 'RPS4Y1', 'sex') # perfect separation
quickViolin(gse68310Flu, 'XIST', 'sex') # perfect separation

# Fix formatted name
gse68310Flu$formattedName = 'GSE68310 - Influenza'

#............................#
##### Fix up annotations #####
#............................#

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


#..............#
##### Save #####
#..............#

save(gse68310Flu, file = 'gse68310_flu.RData')
