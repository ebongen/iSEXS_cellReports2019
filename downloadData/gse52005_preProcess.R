
GSE52005 = getGEOData("GSE52005")
GSE52005 = GSE52005$originalData$GSE52005

GSE52005$rawPheno = GSE52005$pheno
GSE52005$pheno = cleanUpPheno(GSE52005$rawPheno, T)

GSE52005$pheno$age = as.numeric(as.character(GSE52005$pheno$`age_(years)`))


phenoBase = subset(GSE52005$pheno, time_point == "D0")
hist(phenoBase$age)
table(phenoBase$age, phenoBase$sex)


# Fix expr
GSE52005$expr = GSE52005$expr - min(GSE52005$expr) + 1
GSE52005$expr = log2(GSE52005$expr)

gse52005_base = subsetGEMFromPheno(GSE52005, phenoBase)

boxplot(log2(gse52005_base$expr), main = "Baseline")

##### Check Sex Labels #####

"XIST" %in% gse52005_base$keys
"RPS4Y1" %in% gse52005_base$keys

quickViolin(gse52005_base, "RPS4Y1", "sex")
quickViolin(gse52005_base, "XIST", "sex")

timecourse_iSEXS(gse52005_base, sexMetaObj$filterResults$FDR0.05_es0.4_nStudies2_looaFALSE_hetero0, groupColumn = "sex", timeColumn = "age", timeUnits = "years", main = "Full iSEXS")
timecourse_iSEXS(gse52005_base, sexMetaObj$filterResults$xy, groupColumn = "sex", timeColumn = "age", timeUnits = "years", main = "XY-iSEXS")
timecourse_iSEXS(gse52005_base, sexMetaObj$filterResults$autosomeOnly, groupColumn = "sex", timeColumn = "age", timeUnits = "years", main = "Autosomal-iSEXS")


save(GSE52005, gse52005_base, file = "/labs/khatrilab/ebongen/sexDifferences2.0/0_datasets/2_age/prePuberty/gse52005_vaccination/gse52005.RData")
