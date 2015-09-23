###################
#Build commands   #
###################
# cd /shared/cs_researcher/He_C/QH_HOPP/Wade_C/oncoscape/GeneSetBinomialMethods_package/GeneSetBinomialMethods/inst/STTRDataPackages/
# R CMD build SttrDataPackage

cd /shared/cs_researcher/He_C/QH_HOPP/Wade_C/oncoscape/GeneSetBinomialMethods_package/
#R CMD Rd2pdf GeneSetBinomialMethods --force --title='GeneSetBinomialMethods'
R CMD build GeneSetBinomialMethods
R CMD check GeneSetBinomialMethods

install.packages(pkgs = "/shared/cs_researcher/He_C/QH_HOPP/Wade_C/oncoscape/GeneSetBinomialMethods_package/GeneSetBinomialMethods_0.99.1.tar.gz", repos = NULL, type="source")
library(GeneSetBinomialMethods)
BiocGenerics:::testPackage("GeneSetBinomialMethods")

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