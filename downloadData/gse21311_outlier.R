# October 2018 - What's up with that gse21311 outlier?

# Purpose: 
#    - Get an idea of what's happening with the outlier in gse21311?
#    - Is it Klinefelter's syndrome?

# Results: 
#    - Outlier has female-like levels of XIST and male-like levels of RPS4Y1
#         - Suggesting that either XXY sort of phenotype, or technical failure
#    - Y-Chr genes alone perfectly separate males and females
#    - Hard to find non-XIST genes that separate males from females

# Conclusion:
#    - Since Y-chr genes still separate males from females, it suggests that
#      the arrays did not fail
#    - Since outlier as female-like XIST and male-like Y-Chr expression, suggests XXY
#    - But, I'm not willing to definitively say


rocPlot(xySig, valMetaObj$originalData$gse21311) 

violinPlot(xySig, valMetaObj$originalData$gse21311, 'sex') +
  theme_cowplot(font_size = 18) + scale_color_manual(values = c(femColor, malColor))

# Nope, not there
'KDM5D' %in% valMetaObj$originalData$gse21311$keys

quickScatter(valMetaObj$originalData$gse21311, 'RPS4Y1', 'XIST', 'sex') +
  scale_color_manual(values = c(femColor, malColor)) + 
  theme_cowplot(font_size = 18)

# one male  has female-like XIST
quickViolin(valMetaObj$originalData$gse21311, 'XIST', 'sex')+ 
  theme_cowplot(font_size = 18) + theme(legend.position = 'none') + xlab('')

### What if we make a score from Y only?
ySig = xySig
ySig$posGeneNames = c('POOPSMGEE')
ySig$negGeneNames = ySig$negGeneNames[ISEXSLoc[ySig$negGeneNames] == 'Y']


violinPlot(ySig, valMetaObj$originalData$gse21311, 'sex') +
  theme_cowplot(font_size = 18) + scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Y-Chromosome Gene Score')




violinPlot(autoSig, valMetaObj$originalData$gse21311, 'sex') +
  theme_cowplot(font_size = 18) + scale_color_manual(values = c(femColor, malColor)) +
  xlab('') + ylab('Autosomal-iSEXS Score')
