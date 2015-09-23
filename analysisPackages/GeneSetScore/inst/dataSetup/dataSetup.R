##########################
#GeneSetBinomialMethods  #
##########################
setwd("/shared/cs_researcher/He_C/QH_HOPP/Wade_C/oncoscape/code/GeneSetBinomialMethods2_package/")
install.packages(pkgs = "GeneSetBinomialMethods_0.99.1.tar.gz", repos = NULL, type="source")
library(GeneSetBinomialMethods)

###################
#STTR Data        #
###################
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

install.packages(pkgs = system.file(package="GeneSetBinomialMethods", "STTRDataPackages", "PatientHistory_0.99.4.tar.gz"), repos = NULL, type="source")
install.packages(pkgs = system.file(package="GeneSetBinomialMethods", "STTRDataPackages", "SttrDataPackage_0.99.8.tar.gz"), repos = NULL, type="source")
install.packages(pkgs = system.file(package="GeneSetBinomialMethods", "STTRDataPackages", "TCGAgbm_0.99.9.tar.gz"), repos = NULL, type="source")
install.packages(pkgs = system.file(package="GeneSetBinomialMethods", "STTRDataPackages", "TCGAlgg_0.99.5.tar.gz"), repos = NULL, type="source")
install.packages(pkgs = system.file(package="GeneSetBinomialMethods", "STTRDataPackages", "DEMOdz_0.99.9.tar.gz"), repos = NULL, type="source")

library(PatientHistory)
library(SttrDataPackage)
library(TCGAgbm)
library(TCGAlgg)
library(DEMOdz)

TCGAgbm_data <- TCGAgbm()
TCGAmrna <- matrices(TCGAgbm_data)$mtx.mrna
TCGAhistory <- history(TCGAgbm_data)@eventtable

rownames(TCGAmrna) <- substr(rownames(TCGAmrna), 1, nchar(rownames(TCGAmrna))-3)

randomGeneset <- list()
set.seed(123)
for(i in 1:10) {
  randomGeneset[[i]] <- sample(colnames(TCGAmrna), sample.int(10, 1), replace = TRUE)
}
names(randomGeneset)<- paste("geneset", 1:10, sep = "")

STTRdata <- GeneSetBinomialMethods(expressionData = TCGAmrna, clinicalData = TCGAhistory, geneSets = randomGeneset)

GeneSetBinomialMethods::show(STTRdata)

save(STTRdata, file = "GeneSetBinomialMethods/data/STTRdata.RData")

######################
#Original data       #
######################
load(system.file(package="GeneSetBinomialMethods", "dataSetup/rawData", "tbl.mrnaUnified.TCGA.GBM.RData"))
load(system.file(package="GeneSetBinomialMethods", "dataSetup/rawData", "tbl.ptHistory.RData"))
load(system.file(package="GeneSetBinomialMethods", "dataSetup/rawData", "msigdb.RData"))

rownames(tbl.clinical) <- tbl.clinical$ID
tbl.clinical <- tbl.clinical[, -1]

ORIGdata <- GeneSetBinomialMethods(expressionData = tbl.mrna, clinicalData = tbl.clinical, geneSets = genesets)

GeneSetBinomialMethods::show(ORIGdata)

save(ORIGdata, file = "GeneSetBinomialMethods/data/ORIGdata.RData")
