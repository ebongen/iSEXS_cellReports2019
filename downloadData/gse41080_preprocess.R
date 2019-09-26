# July 6th 2018 - Download Furman's sex differences in vaccination dataset

# Purpose: 
#   - Does A-iSEXS correlate with T levels in males?

#....................#
##### Load Stuff #####
#....................#
setwd('/labs/khatrilab/ebongen/sexDifferences2.0/')

library(MetaIntegrator)

source('00_tools/general_GEMfunx.R')
source('~/Tools/Graphing Scripts/quickPlots.R')

setwd('0_datasets/6_infection/vaccine/gse41080/')

#................................#
##### Load Purvesh's Version #####
#................................#
# Purpose:
#   - No preprocessed data available
#   - Purvesh has already normalized it, so I'll use his version

# Load all of Purvesh's influenza datasets
load("/labs/khatrilab/pk_projects/Influenza/data/input/Influenza.RData")

### Keys - Look good!
gse41080GEM$keys[1:10]

### Class
myClass = rep(0, nrow(gse41080GEM$pheno))
names(myClass) = rownames(gse41080GEM$pheno)
gse41080GEM$class = myClass

### Add formatted name
gse41080GEM$formattedName = 'GSE41080'

### Make expr be a matrix
is.data.frame(gse41080GEM$expr)
gse41080GEM$expr = as.matrix(gse41080GEM$expr)

checkDataObject(gse41080GEM, 'Dataset') # Passes!
#..............#
##### Expr #####
#..............#
#   - No negatvives
#   - No NAs
#   - No obvious batch effect
range(gse41080GEM$expr) # 6-16
ncol(gse41080GEM$expr) # 91 people

# No obvious batch effect
png('gse41080.png', width = 1000)
boxplot(gse41080GEM$expr, main = 'GSE41080')
dev.off()
#...............#
##### Pheno #####
#...............#
# Sex lables
table(gse41080GEM$pheno$sex)
sex = ifelse(gse41080GEM$pheno$sex == 'Female',yes = 'female', no = 'male')
table(gse41080GEM$pheno$sex, sex) # all match
gse41080GEM$pheno$sex = sex

# Age labels
# bimodal
# one group <30, another group 60-100
hist(gse41080GEM$pheno$age.years) 
gse41080GEM$pheno$age = gse41080GEM$pheno$age.years

# Incorperate testosterone data
testo = read.csv('testo_data.csv')
View(testo)

# Create a testosterone dictionary
testoDict = as.numeric(testo$T_nmol_l)
names(testoDict) = as.character(testo$xID)

# Make sure names match
sum(names(testoDict) %in% as.character(gse41080GEM$pheno$id)) # 78 are present
names(testoDict)[!names(testoDict) %in% as.character(gse41080GEM$pheno$id)]

testosterone = testoDict[as.character(gse41080GEM$pheno$id)]
cbind(names(testosterone), as.character(gse41080GEM$pheno$id))
gse41080GEM$pheno$testosterone = testosterone

#..........................#
##### Check sex Labels #####
#..........................#
# Purpose: 
#   - Make sure sex labels are correct
#
# Results: 
#    - 2 samples with incorrect labels
#    - Probably 1 male and 1 female swapped
#    - I removed those two
'RPS4Y1' %in% gse41080GEM$keys #T
'KDM5D' %in% gse41080GEM$keys # F
'XIST' %in% gse41080GEM$keys # T

# Clearly 1 male and 1 female got swapped
quickViolin(gse41080GEM, 'RPS4Y1', 'sex') # 1 male and 1 female swapped
quickViolin(gse41080GEM, 'XIST', 'sex') # one obvious male is female, one extra low female
quickScatter(gse41080GEM, 'RPS4Y1', 'XIST', 'sex')

# Remove weirdos
imputedSex = imputeSex(gse41080GEM)
weirdos = rownames(gse41080GEM$pheno)[gse41080GEM$pheno$sex != imputedSex]
weirdos # 2 samples

for(mySamp in weirdos){
  gse41080GEM = removeOneSample(gse41080GEM, mySamp)
}

# Clean separation!
quickViolin(gse41080GEM, 'RPS4Y1', 'sex')

#...................................#
##### Check Testosterone Levels #####
#...................................#
# Purpose:
#   - Make sure Testosterone is higher in males than females
#   - Sanity check for Furman's Testosterone level measurements
# 
# Results:
#   - Significantly higher Testosterone in males than females

myData = gse41080GEM$pheno

ggplot(myData, aes(x=sex, y = testosterone)) + geom_violin(fill='grey', trim=F)+
  geom_jitter(width=0.1, size=2) + theme_cowplot() + 
  draw_label('p<0.0001', x=1.5, y=20) + ggtitle("GSE41080 Testosterone Levels") +
  xlab('') +ylab('Serum Free Testosterone (nmol/L)')
t.test(myData$testosterone[myData$sex=='male'], myData$testosterone[myData$sex == 'female'])


a = ggplot(myData, aes(x=age.years, y=testosterone, col=sex)) + geom_smooth()+
  geom_point(size=2) + scale_color_manual(values = c(femColor, malColor))+
  theme_cowplot() + ggtitle('GSE41080 - All Individuals') +
  xlab('Age (years)') +ylab('Serum Free Testosterone (nmol/L)') +
  theme(legend.position = c(0.8, 0.8))

#....................................................#
##### Does Testosterone have an effect on iSEXS? #####
#....................................................#
# Purpose: 
#   - Testosterone levels have been related to vaccine response
#   - I'm trying to see if A-iSEXS is connected to vaccine response
#   - Could be through testosterone levels!
#
# Results: 
#    - No correlation between A-iSEXS and testosterone levels
load('../../../../1_metaAnaly/sexMetaObj.RData')

myData$autoScore = calculateScore(sexMetaObj$filterResults$autosomeOnly,gse41080GEM)
myData$xyScore = calculateScore(sexMetaObj$filterResults$xy, gse41080GEM)

myDataYoung = subset(myData, age< 60)
hist(myDataYoung$age)

myDataYoung_fem = subset(myDataYoung, sex == 'female')
myDataYoung_mal = subset(myDataYoung, sex == 'male')

### No significant relationship between T and auto-iSEXS
b = ggplot(myDataYoung, aes(x=testosterone, y=autoScore, col=sex)) +
  geom_point(size=2)+ geom_smooth(method = 'lm')+
  scale_color_manual(values = c(femColor, malColor)) +
  theme_cowplot() +
  xlab('Serum Free Testosterone (nmol/L)') + ylab('Autosomal-iSEXS Score') + 
  ggtitle('GSE41080 - Individuals < 60 years old') + 
  draw_label('N.S.', x = 2.5, y=2.5, colour=femColor) +
  draw_label('N.S.', x=20, y=1, colour = malColor) +
  theme(legend.position = c(0.8, 0.8))
cor.test(myDataYoung$autoScore, myDataYoung$testosterone) # p=0.12
cor.test(myDataYoung_fem$autoScore, myDataYoung_fem$testosterone) # p=0.9
cor.test(myDataYoung_mal$autoScore, myDataYoung_mal$testosterone) # p=0.4


### No significant relationship between T and XY-iSEXS within males or females
ggplot(myDataYoung, aes(x=testosterone, y=xyScore, col=sex)) +
  geom_point(size=2)+ geom_smooth(method = 'lm') +
  scale_color_manual(values = c(femColor, malColor)) +
  theme_cowplot()
cor.test(myDataYoung$xyScore, myDataYoung$testosterone) # p=0.0007, probs by definition, not surprising
cor.test(myDataYoung_fem$xyScore, myDataYoung_fem$testosterone) # p=0.3
cor.test(myDataYoung_mal$xyScore, myDataYoung_mal$testosterone) # p=0.8


### Look at all ages combined
ggplot(myData, aes(x=testosterone, y=autoScore, col=sex)) +
  geom_point(size=2.5, aes(shape=age.years < 60))+
  scale_color_manual(values = c(femColor, malColor))+
  scale_shape_manual(values = c(1,19)) + 
  theme_cowplot()

### Look at old people
myDataOld = subset(myData, age > 60)
myDataOld_fem = subset(myDataOld, sex=='female')
myDataOld_mal = subset(myDataOld, sex == 'male')

ggplot(myDataOld, aes(x=testosterone, y=autoScore, col=sex)) +
  geom_point(size=2)+
  scale_color_manual(values = c(femColor, malColor)) +
  theme_cowplot()
cor.test(myDataOld_fem$testosterone, myDataOld_fem$autoScore) # p=0.3
cor.test(myDataOld_mal$testosterone, myDataOld_mal$autoScore) # p=0.6

### Save the plot
Tplot = plot_grid(a,b, nrow=1, ncol=2)
save_plot('gse41080_testosterone.pdf', Tplot, ncol=2, nrow=1)

#..............#
##### Save #####
#..............#
save(gse41080GEM, file = 'gse41080.RData')
