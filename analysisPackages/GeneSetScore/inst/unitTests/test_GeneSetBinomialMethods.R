########################
#Required Libraries    #
########################
library(org.Hs.eg.db)
library(PatientHistory)
library(SttrDataPackage)
library(TCGAgbm)
library(TCGAlgg)
library(DEMOdz)
library(RUnit)
library(foreach)

#------------------------------------------------------------------------------------------------------------------------#
printf = function (...) print (noquote (sprintf (...)))

#------------------------------------------------------------------------------------------------------------------------#
runTests <- function() {
  test_GeneSetBinomialMethods()
  test_show.GeneSetBinomialMethods()
  test_accessors()
  test_randomSample()
  test_analysisDataSetup()
  test_geneSetScoreTest()
} # runTests

#------------------------------------------------------------------------------------------------------------------------#
test_GeneSetBinomialMethods = function () {
  printf("--- test_GeneSetBinomialMethods")
  
  checkEquals(slotNames("GeneSetBinomialMethods"), c("expressionData", "clinicalData", "geneSets"))
  
  #######################
  #TCGA External Data   #
  #######################
  #No shared observations in the gene expression data
  TCGAgbm_data <- TCGAgbm()
  TCGAmrna <- matrices(TCGAgbm_data)$mtx.mrna
  TCGAhistory <- history(TCGAgbm_data)@eventtable
  
  #No shared samples in expression and clinical data
  test1 <- try(GeneSetBinomialMethods(expressionData = TCGAmrna, clinicalData = TCGAhistory, geneSets = list()), silent = TRUE)
  test1 <- make.names(test1[1])
  checkEquals(test1, "Error.in.validObject..Object.......invalid.class..GeneSetBinomialMethods..object..There.are.no.shared.samples.between.the.expression.data.and.clinical.data..")
  
  #Make the sample names the same
  rownames(TCGAmrna) <- substr(rownames(TCGAmrna), 1, nchar(rownames(TCGAmrna))-3)
  
  #Duplicated observations in the gene expression data (since the clinical data is a data frame we don't need to do the same check)
  TCGAmrna.dupExpr = TCGAmrna
  rownames(TCGAmrna.dupExpr)[c(1,4, 5, 6)] <- "duplicated"
  
  test2 <- try(GeneSetBinomialMethods(expressionData = TCGAmrna.dupExpr, clinicalData = TCGAhistory, geneSets = list()), silent = TRUE)
  test2 <- make.names(test2[1])
  checkEquals(test2, "Error.in.validObject..Object.......invalid.class..GeneSetBinomialMethods..object..None.of.the.methods.for.objects.of.class.GeneSetBinomialMethods.can.handle.duplicated.observations...Please.input.expression.data.that.does.not.have.duplicated.observations..")
  
  #Warning if none of the genes in the geneset are in the expression data
  badGenes <- paste("gene", 1:50, sep = "")
  
  badGeneset <- list()
  set.seed(123)
  for(i in 1:10) {
    badGeneset[[i]] <- sample(badGenes, sample.int(10, 1), replace = TRUE)
  }
  
  names(badGeneset)<- paste("geneset", 1:10, sep = "")
  
  badGenesetObj <- .withWarnings(GeneSetBinomialMethods(expressionData = TCGAmrna, clinicalData = TCGAhistory, geneSets = badGeneset))
  
  checkEquals(make.names(badGenesetObj[[2]][[1]]), "simpleWarning.in.validityMethod.object...None.of.the.genes.in.the.specified.geneset.are.in.the.supplied.expression.data..")
  
  #Create a random gene set that works.
  externalDataGeneset <- list()
  set.seed(123)
  for(i in 1:10) {
    externalDataGeneset[[i]] <- sample(colnames(TCGAmrna), sample.int(10, 1), replace = TRUE)
  }
  names(externalDataGeneset)<- paste("geneset", 1:10, sep = "")
  
  #Create a valid object
  externalData <- GeneSetBinomialMethods(expressionData = TCGAmrna, clinicalData = TCGAhistory, geneSets = externalDataGeneset)
  checkTrue(validObject(externalData))
  
  checkEquals(dim(externalData@expressionData), c(154, 20457))
  checkEquals(dim(externalData@clinicalData), c(592, 409)) 
  checkTrue(length(externalData@geneSets)== 10)
  
  checkEquals(externalData@geneSets[[1]], c("SLC25A17", "IL36G", "TMEM100"))
  checkEquals(externalData@geneSets[[4]], c("MRPL54", "NUFIP1", "FAM149B1", "CBLL1", "ZBTB39", "TP73", "PRSS30P", "SLC39A7" ))
  checkEquals(externalData@geneSets[[8]], c("TAF13", "ARHGAP42", "KIFC1"))
  
  data(STTRdata)
  checkIdentical(externalData, STTRdata)
  
  #######################
  #Package default data #
  #######################
  load(system.file(package="GeneSetBinomialMethods", "dataSetup/rawData", "tbl.mrnaUnified.TCGA.GBM.RData"))
  load(system.file(package="GeneSetBinomialMethods", "dataSetup/rawData", "tbl.ptHistory.RData"))
  load(system.file(package="GeneSetBinomialMethods", "dataSetup/rawData", "msigdb.RData"))
  
  rownames(tbl.clinical) <- tbl.clinical$ID
  tbl.clinical <- tbl.clinical[, -1]
  
  internalData <- GeneSetBinomialMethods(expressionData = tbl.mrna, clinicalData = tbl.clinical, geneSets = genesets)
  
  checkEquals(dim(internalData@expressionData), c(315, 11864))
  checkEquals(dim(internalData@clinicalData), c(583, 12)) 
  checkTrue(length(internalData@geneSets) == 10295)
  checkEquals(length(which(lapply(internalData@geneSets, length) <= 10)), 575)
  
  checkEquals(names(head(internalData@geneSets)),   c("NUCLEOPLASM", "EXTRINSIC_TO_PLASMA_MEMBRANE", "ORGANELLE_PART", "CELL_PROJECTION_PART", "CYTOPLASMIC_VESICLE_MEMBRANE", "GOLGI_MEMBRANE"))
  checkEquals(names(head(which(lapply(internalData@geneSets, length) < 10))), c("INTERCALATED_DISC", "V$SEF1_C", "TCGATGG,MIR-213", "CGGTGTG,MIR-220", "ACGCACA,MIR-210", "GTCGATC,MIR-369-5P"))
  
  data(ORIGdata)
  checkIdentical(internalData, ORIGdata)
  
  ###########################################################################
  #Silly data! (checks the warning that the IDs can be coerced to numeric)  #
  ###########################################################################
  simpleMat <- matrix(c(1,2,3,4,5,6,7,8,9), ncol = 3)
  rownames(simpleMat) <- 1:3
  colnames(simpleMat) <- c("Moe", "Larry", "Curly")
  
  simpleDF <- data.frame(matrix(c(1,2,3,4,5,6,7,8,9), ncol = 3))
  rownames(simpleDF) <- as.numeric(rownames(simpleDF))
  
  checkTrue(!suppressWarnings(any(is.na(as.numeric(rownames(simpleMat))))))
  checkTrue(!suppressWarnings(any(is.na(as.numeric(rownames(simpleDF))))))
  
  sillyData <- .withWarnings(GeneSetBinomialMethods(expressionData = simpleMat, clinicalData = simpleDF, geneSets = list()))
  
  checkEquals(make.names(sillyData[[2]][[1]]), "simpleWarning.in.validityMethod.object...The.rownames.of.both.the.expression.data.and.clinical.data.are.numeric...Please.verify.these.are.the.unique.ID.values.you.want.to.use..")
} # test_GeneSetBinomialMethods

#------------------------------------------------------------------------------------------------------------------------
test_show.GeneSetBinomialMethods = function() 
{
  printf("--- test_show.GeneSetBinomialMethods")
  
  data(ORIGdata)
  ORIGdata.chk <- show(ORIGdata)
  checkEquals(ORIGdata.chk, "This is an object of class GeneSetBinomialMethods.\n The expression data has 315 observations and 11864 genes (getExpressionData()).\n The clinical data has 583 observations and 12 variables (getClinicalData()).\n There are 10295 gene sets (getGeneSets()).")
  
  data(STTRdata)
  STTRdata.chk <- show(STTRdata)
  checkEquals(STTRdata.chk, "This is an object of class GeneSetBinomialMethods.\n The expression data has 154 observations and 20457 genes (getExpressionData()).\n The clinical data has 592 observations and 409 variables (getClinicalData()).\n There are 10 gene sets (getGeneSets()).")
} # test_show.GeneSetBinomialMethods

#------------------------------------------------------------------------------------------------------------------------
test_accessors = function () {
  printf("--- test_accessors")
  
  #Original Data
  data(ORIGdata)
  
  checkEquals(dim(getExpressionData(ORIGdata)), c(315, 11864))
  checkEquals(dim(getClinicalData(ORIGdata)), c(583, 12)) 
  checkTrue(length(getGeneSets(ORIGdata)) == 10295)
  checkEquals(length(which(lapply(getGeneSets(ORIGdata), length) <= 10)), 575)
  
  checkEquals(names(head(getGeneSets(ORIGdata))),   c("NUCLEOPLASM", "EXTRINSIC_TO_PLASMA_MEMBRANE", "ORGANELLE_PART", "CELL_PROJECTION_PART", "CYTOPLASMIC_VESICLE_MEMBRANE", "GOLGI_MEMBRANE"))
  checkEquals(names(head(which(lapply(getGeneSets(ORIGdata), length) < 10))), c("INTERCALATED_DISC", "V$SEF1_C", "TCGATGG,MIR-213", "CGGTGTG,MIR-220", "ACGCACA,MIR-210", "GTCGATC,MIR-369-5P"))
  
  #STTR Data
  data(STTRdata)
  
  checkEquals(dim(getExpressionData(STTRdata)), c(154, 20457))
  checkEquals(dim(getClinicalData(STTRdata)), c(592, 409)) 
  checkTrue(length(getGeneSets(STTRdata))== 10)
} # test_accessors

#------------------------------------------------------------------------------------------------------------------------
test_randomSample = function()
{
  printf("--- test_randomSample")
  
  ########################################################
  #Original Data: For a fixed seed get the same samples  #
  ########################################################
  set40.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, nG1 = 40, nG2 = 40, all = FALSE, seed = 123)
  
  
  checkEquals(unlist(set40.orig$shortSurvivors), c("TCGA.16.1047", "TCGA.08.0359", "TCGA.06.2557", "TCGA.08.0345", "TCGA.28.5218", "TCGA.14.0787", "TCGA.06.A5U1", "TCGA.06.0127", "TCGA.06.A6S1",
                                                   "TCGA.06.2567", "TCGA.08.0373", "TCGA.06.2559", "TCGA.06.0173", "TCGA.06.A6S0", "TCGA.06.5412", "TCGA.08.0392", "TCGA.76.6656", "TCGA.02.0060",
                                                   "TCGA.14.1453", "TCGA.19.1790", "TCGA.08.0510", "TCGA.41.3392", "TCGA.OX.A56R", "TCGA.32.4210", "TCGA.41.3393", "TCGA.41.2571", "TCGA.27.2526",
                                                   "TCGA.06.0140", "TCGA.76.6657", "TCGA.06.6391", "TCGA.06.0413", "TCGA.06.0162", "TCGA.76.4928", "TCGA.06.0157", "TCGA.02.0037", "TCGA.14.1795",
                                                   "TCGA.19.2621", "TCGA.76.6193", "TCGA.12.0657", "TCGA.81.5910"))
  
  checkEquals(unlist(set40.orig$longSurvivors), c("TCGA.02.0114", "TCGA.08.0245", "TCGA.27.2528", "TCGA.41.2572", "TCGA.28.5207", "TCGA.02.0046", "TCGA.06.0138", "TCGA.12.0776", "TCGA.06.0160",
                                                  "TCGA.06.0133", "TCGA.28.5204", "TCGA.06.0129", "TCGA.02.0337", "TCGA.06.0168", "TCGA.14.0789", "TCGA.32.2491", "TCGA.16.1045", "TCGA.02.0038",
                                                  "TCGA.12.1600", "TCGA.08.0358", "TCGA.08.0389", "TCGA.02.0332", "TCGA.02.0106", "TCGA.26.5135", "TCGA.02.0260", "TCGA.02.0338", "TCGA.19.2629",
                                                  "TCGA.06.0165", "TCGA.02.0075", "TCGA.19.5952", "TCGA.32.4719", "TCGA.08.0355", "TCGA.02.0289", "TCGA.08.0511", "TCGA.02.0014", "TCGA.06.0124",
                                                  "TCGA.02.0446", "TCGA.12.1094", "TCGA.02.0089", "TCGA.14.0736"))
  
  checkTrue(all(getClinicalData(ORIGdata)[which(rownames(getClinicalData(ORIGdata)) %in% unlist(set40.orig$shortSurvivors)), "survival"] <= 0.5))
  checkTrue(all(getClinicalData(ORIGdata)[which(rownames(getClinicalData(ORIGdata)) %in% unlist(set40.orig$longSurvivors)), "survival"] > 0.5))
  
  ########################################################
  #STTR Data: For a fixed seed get the same samples  #
  ########################################################
  data(STTRdata)
  
  set40_50.sttr <-  randomSample(object = STTRdata, survivalVar = "Survival", cut = 100, nG1 = 40, nG2 = 50, all = FALSE, seed = 123)
  
  checkTrue(all(getClinicalData(STTRdata)[unlist(set40_50.sttr$shortSurvivors), "Survival"] <= 100))
  checkTrue(all(getClinicalData(STTRdata)[unlist(set40_50.sttr$longSurvivors), "Survival"] > 100))
  
  checkTrue(length(unlist(set40_50.sttr$shortSurvivors)) == 40)
  checkTrue(length(unlist(set40_50.sttr$longSurvivors)) == 50)
  
  checkEquals(unlist(set40_50.sttr$shortSurvivors), c("TCGA.14.1794", "TCGA.14.3476", "TCGA.06.2569", "TCGA.06.5418", "TCGA.26.5139", "TCGA.14.0813", "TCGA.06.0119",
                                                      "TCGA.08.0345", "TCGA.76.4928", "TCGA.27.2526", "TCGA.08.0352", "TCGA.06.0749", "TCGA.06.A6S1", "TCGA.32.1976",
                                                      "TCGA.14.0781", "TCGA.06.0413", "TCGA.14.1395", "TCGA.14.0787", "TCGA.19.0962", "TCGA.41.4097", "TCGA.06.0201",
                                                      "TCGA.06.0140", "TCGA.28.2502", "TCGA.02.0439", "TCGA.28.2499", "TCGA.12.5301", "TCGA.28.1757", "TCGA.12.5299",
                                                      "TCGA.32.1977", "TCGA.19.5955", "TCGA.19.2624", "TCGA.06.A5U1", "TCGA.28.1747", "TCGA.06.0402", "TCGA.02.0055",
                                                      "TCGA.14.0812", "TCGA.19.2621", "TCGA.76.6192", "TCGA.14.0862", "TCGA.06.0145"))
  
  checkEquals(unlist(set40_50.sttr$longSurvivors), c("TCGA.19.1392", "TCGA.08.0244", "TCGA.02.2483", "TCGA.06.0644", "TCGA.28.5216", "TCGA.06.0875", "TCGA.12.3646",
                                                     "TCGA.06.0645", "TCGA.19.1787", "TCGA.27.2518", "TCGA.28.5220", "TCGA.27.1837", "TCGA.02.0281", "TCGA.32.2632",
                                                     "TCGA.06.5410", "TCGA.06.0152", "TCGA.12.1099", "TCGA.02.0057", "TCGA.02.0089", "TCGA.26.5135", "TCGA.32.4210",
                                                     "TCGA.02.0289", "TCGA.06.0195", "TCGA.76.4931", "TCGA.06.0211", "TCGA.02.0325", "TCGA.RR.A6KC", "TCGA.41.2573",
                                                     "TCGA.16.1056", "TCGA.28.6450", "TCGA.06.5416", "TCGA.32.2495", "TCGA.02.0106", "TCGA.08.0517", "TCGA.02.0028",
                                                     "TCGA.12.0691", "TCGA.06.0149", "TCGA.76.6657", "TCGA.26.1443", "TCGA.06.1086", "TCGA.19.5958", "TCGA.26.1799",
                                                     "TCGA.19.1788", "TCGA.06.1805", "TCGA.76.4925", "TCGA.19.5952", "TCGA.06.1084", "TCGA.06.2564", "TCGA.14.0736",
                                                     "TCGA.02.0051"))
  
  checkTrue(all(getClinicalData(STTRdata)[which(rownames(getClinicalData(STTRdata)) %in% unlist(set40_50.sttr$shortSurvivors)), "survival"] <= 100))
  checkTrue(all(getClinicalData(STTRdata)[which(rownames(getClinicalData(STTRdata)) %in% unlist(set40_50.sttr$longSurvivors)), "survival"] > 100))
  
  ##########################
  #nG1 or nG2 is truncated #
  ##########################  
  nG1_truncated_orig <- suppressWarnings(randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.1, nG1 = 30, nG2 = 20, all = FALSE, seed = 123))
  checkEquals(length(nG1_truncated_orig$shortSurvivors), 25)
  checkEquals(length(nG1_truncated_orig$longSurvivors), 20)
  
  nG2_truncated_orig <- suppressWarnings(randomSample(object = ORIGdata, survivalVar = "survival", cut = 10, nG1 = 30, nG2 = 20, all = FALSE, seed = 123))
  checkEquals(length(nG2_truncated_orig$shortSurvivors), 30)
  checkEquals(length(nG2_truncated_orig$longSurvivors), 2)
  
  nG1_nG2_truncated_orig <- suppressWarnings(randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, nG1 = 500, nG2 = 500, all = FALSE, seed = 123))
  checkEquals(length(nG1_nG2_truncated_orig$shortSurvivors), 130)
  checkEquals(length(nG1_nG2_truncated_orig$longSurvivors), 313)
  
  nG1_nG2_truncated_sttr <- suppressWarnings(randomSample(object = STTRdata, survivalVar = "Survival", cut = 100, nG1 = 1000, nG2 = 1000, all = FALSE, seed = 123))
  checkEquals(length(nG1_nG2_truncated_sttr$shortSurvivors), 80)
  checkEquals(length(nG1_nG2_truncated_sttr$longSurvivors), 512)
  
  ########################################################
  #Get all of the subjects when nG1 and nG2 are the same #
  ########################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.75, all = TRUE)
  checkEquals(length(unlist(setAll.orig$shortSurvivors)), 171)
  checkEquals(length(unlist(setAll.orig$longSurvivors)), 272)
  
  checkEquals(length(getClinicalData(ORIGdata)[which(rownames(getClinicalData(ORIGdata)) %in% unlist(setAll.orig$shortSurvivors)), "survival"]), 171)
  checkEquals(length(getClinicalData(ORIGdata)[which(rownames(getClinicalData(ORIGdata)) %in% unlist(setAll.orig$longSurvivors)), "survival"]), 272)
  
  setAll.sttr <- randomSample(object = STTRdata, survivalVar = "Survival", cut = 200, all = TRUE)
  checkEquals(length(unlist(setAll.sttr$shortSurvivors)), 174)
  checkEquals(length(unlist(setAll.sttr$longSurvivors)), 418)
  
  checkEquals(length(getClinicalData(STTRdata)[which(rownames(getClinicalData(STTRdata)) %in% unlist(setAll.sttr$shortSurvivors)), "Survival"]), 174)
  checkEquals(length(getClinicalData(STTRdata)[which(rownames(getClinicalData(STTRdata)) %in% unlist(setAll.sttr$longSurvivors)), "Survival"]), 418)
  
  ##############################################
  #Verify that random samples are not the same #
  ##############################################
  checkSeed <- vector(mode = "logical")
  checkDiffShort <- vector(mode = "logical")
  checkDiffLong <- vector(mode = "logical")
  
  for(i in 1:5) {
    setRand1 <- randomSample(object = ORIGdata, survivalVar = "survival", nG1 = 20, nG2 = 30, cut = 0.5, all = FALSE)
    setRand2 <- randomSample(object = ORIGdata, survivalVar = "survival", nG1 = 20, nG2 = 20, cut = 0.5, all = FALSE)
    
    checkSeed[[i]] = (setRand1$seed != setRand2$seed)
    checkDiffShort[[i]] = any(setdiff(unlist(setRand1$shortSurvivors), unlist(setRand2$shortSurvivors)) > 0)
    checkDiffLong[[i]] = any(setdiff(unlist(setRand1$longSurvivors), unlist(setRand2$longSurvivors)) > 0)
  }
  
  checkEquals(length(unlist(setRand1$shortSurvivors)), 20)
  checkEquals(length(unlist(setRand1$longSurvivors)), 30)
  
  checkTrue(any(checkSeed == TRUE))
  checkTrue(any(checkDiffShort == TRUE))
  checkTrue(any(checkDiffLong == TRUE))
  
  ###############
  #Null values  #
  ###############
  setErr1 <- try(randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = FALSE), silent = TRUE)
  checkEquals(setErr1[1], "Error in .local(object, ...) : \n  If either of the arguments nG1 or nG2 is null, then all must equal TRUE\n")
  
  setErr2 <- try(randomSample(object = ORIGdata, nG1 = NULL, nG2 = 30, cut = 0.5, all = FALSE), silent = TRUE)
  checkEquals(setErr2[1], "Error in .local(object, ...) : \n  If either of the arguments nG1 or nG2 is null, then all must equal TRUE\n")
  
  setErr3 <- try(randomSample(object = STTRdata, nG1 = NULL, nG2 = 30, cut = 100, all = FALSE), silent = TRUE)
  checkEquals(setErr2[1], "Error in .local(object, ...) : \n  If either of the arguments nG1 or nG2 is null, then all must equal TRUE\n")
  
  ###############
  #No Samples   #
  ###############
  checkIdentical(range(getClinicalData(STTRdata)[, "Survival"]), c(3, 3881))
  checkIdentical(range(na.omit(getClinicalData(ORIGdata)[, "survival"])), c(0.01, 10.63))
  
  noSamples1 <- try(randomSample(object = ORIGdata, survivalVar = "survival", nG1 = 30, nG2 = 20, cut = -2, all = FALSE, seed = 123), silent = TRUE)
  checkEquals(noSamples1[1],  "Error in .local(object, ...) : \n  Either one or both groups contain no subjects for the given cut value.  The survival range is [0.01, 10.63].\n")
  
  noSamples2 <- try(randomSample(object = ORIGdata, survivalVar = "survival", nG1 = 10, nG2 = 50, cut = 11, all = FALSE, seed = 123), silent = TRUE)
  checkEquals(noSamples2[1], "Error in .local(object, ...) : \n  Either one or both groups contain no subjects for the given cut value.  The survival range is [0.01, 10.63].\n")
  
  noSamples1 <- try(randomSample(object = STTRdata, survivalVar = "Survival", nG1 = 30, nG2 = 20, cut = -2, all = FALSE, seed = 123), silent = TRUE)
  checkEquals(noSamples1[1],  "Error in .local(object, ...) : \n  Either one or both groups contain no subjects for the given cut value.  The survival range is [3, 3881].\n")
} #test_randomSample

#------------------------------------------------------------------------------------------------------------------------
test_analysisDataSetup = function()
{
  printf("--- test_analysisDataSetup")
  
  data(ORIGdata)
  data(STTRdata)
  
  ################################
  #Original Data: No Covariates  #
  ################################
  set40.orig <- randomSample(object = ORIGdata,  survivalVar = "survival", nG1 = 40, nG2 = 40, cut = 0.5, all = FALSE, seed = 12345)
  
  analysisDataSetup_nocov_orig <- analysisDataSetup(
    object = ORIGdata,
    sampleIDsG1 = set40.orig$shortSurvivors,
    sampleIDsG2 = set40.orig$longSurvivors,
    covariates = NULL,
    geneSet = "KANG_CISPLATIN_RESISTANCE_DN",
    sampleDescription ="TCGA GBM long vs. short survivors (no covariates)",
    geneSetDescription ="msgidb:KANG_CISPLATIN_RESISTANCE_DN")
  
  checkTrue(analysisDataSetup_nocov_orig$sampleDescription == "TCGA GBM long vs. short survivors (no covariates)")
  checkTrue(analysisDataSetup_nocov_orig$geneSetDescription == "msgidb:KANG_CISPLATIN_RESISTANCE_DN")
  checkEquals(dim(analysisDataSetup_nocov_orig$analysisData), c(38, 7))
  
  checkTrue(length(intersect(analysisDataSetup_nocov_orig$unmatchedSamples, analysisDataSetup_nocov_orig$matchedSamples)) == 0)
  
  checkEquals(analysisDataSetup_nocov_orig$matchedGenes, c("ASB8", "COL9A2", "GPR172A", "NBN", "RPS6KA4", "SEMA3C"))
  
  checkEquals(analysisDataSetup_nocov_orig$unmatchedGenes, c("TMX4", "UBE2K"))
  
  checkEquals(analysisDataSetup_nocov_orig$matchedSamples, c("TCGA.02.0001", "TCGA.02.0003", "TCGA.02.0006", "TCGA.02.0023", "TCGA.02.0028", "TCGA.02.0051", "TCGA.02.0060",
                                                             "TCGA.02.0079", "TCGA.02.0106", "TCGA.02.0266", "TCGA.02.0326", "TCGA.02.0333", "TCGA.02.0339", "TCGA.02.0456",
                                                             "TCGA.06.0124", "TCGA.06.0130", "TCGA.06.0141", "TCGA.06.0142", "TCGA.06.0155", "TCGA.06.0173", "TCGA.06.0194",
                                                             "TCGA.06.0197", "TCGA.06.0219", "TCGA.06.0649", "TCGA.08.0345", "TCGA.08.0348", "TCGA.08.0349", "TCGA.08.0358",
                                                             "TCGA.08.0385", "TCGA.08.0509", "TCGA.08.0510", "TCGA.08.0514", "TCGA.08.0516", "TCGA.12.0657", "TCGA.12.0707",
                                                             "TCGA.14.0783", "TCGA.14.1396", "TCGA.14.1401"))
  
  checkEquals(analysisDataSetup_nocov_orig$unmatchedSamples,  c("TCGA.28.1747", "TCGA.76.6193", "TCGA.02.0071", "TCGA.76.4925", "TCGA.32.1980", "TCGA.14.1795", "TCGA.06.2559", "TCGA.14.1794", "TCGA.06.5418",
                                                                "TCGA.19.1790", "TCGA.19.1789", "TCGA.14.3476", "TCGA.41.3392", "TCGA.19.A6J4", "TCGA.76.4934", "TCGA.19.2624", "TCGA.27.2526", "TCGA.14.0862",
                                                                "TCGA.19.1788", "TCGA.12.5299", "TCGA.28.5218", "TCGA.06.A6S0", "TCGA.41.4097", "TCGA.06.5856", "TCGA.41.2572", "TCGA.76.6286", "TCGA.14.0866",
                                                                "TCGA.19.1787", "TCGA.12.0662", "TCGA.06.1806", "TCGA.14.0790", "TCGA.06.2558", "TCGA.27.1835", "TCGA.12.1599", "TCGA.27.1834", "TCGA.12.0769",
                                                                "TCGA.12.0653", "TCGA.76.4932", "TCGA.14.1823", "TCGA.14.1037", "TCGA.06.0165", "TCGA.06.5408"))
  
  checkIdentical(analysisDataSetup_nocov_orig$analysisData[, 1], c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1))
  
  checkIdentical(round(analysisDataSetup_nocov_orig$analysisData[, 4], 2), c(-0.99,  0.06, -0.04, -0.48, -1.20, -1.08,  0.30,  0.24, -0.02,  2.14,  1.19,  1.51, -0.06,  2.08,  0.94, -0.63, -0.02,
                                                                             -1.37,  0.24, -0.39,  0.39,  0.34, -1.01, -0.97,  0.38,  0.32,  0.64,  1.76,  0.63,  0.09,  0.09,  0.07, -0.58,  0.15,
                                                                             -0.62, -0.42,  1.32, -0.95))
  
  ################################
  #STTR Data: No Covariates      #
  ################################
  set60.sttr <- randomSample(object = STTRdata,  survivalVar = "Survival", nG1 = 60, nG2 = 60, cut = 100, all = FALSE, seed = 12345)
  
  analysisDataSetup_nocov_sttr <- analysisDataSetup(
    object = STTRdata,
    sampleIDsG1 = set60.sttr$shortSurvivors,
    sampleIDsG2 = set60.sttr$longSurvivors,
    covariates = NULL,
    geneSet = "geneset8",
    sampleDescription ="STTR data (no covariates)",
    geneSetDescription ="geneset8")
  
  checkTrue(analysisDataSetup_nocov_sttr$sampleDescription == "STTR data (no covariates)")
  checkTrue(analysisDataSetup_nocov_sttr$geneSetDescription == "geneset8")
  checkEquals(dim(analysisDataSetup_nocov_sttr$analysisData), c(28, 4))
  
  checkTrue(length(intersect(analysisDataSetup_nocov_sttr$unmatchedSamples, analysisDataSetup_nocov_sttr$matchedSamples)) == 0)
  
  checkEquals(analysisDataSetup_nocov_sttr$matchedGenes, c("ARHGAP42", "KIFC1", "TAF13"))
  
  checkEquals(length(analysisDataSetup_nocov_sttr$unmatchedGenes), 0)
  
  checkEquals(analysisDataSetup_nocov_sttr$matchedSamples, c("TCGA.02.0047", "TCGA.02.0055", "TCGA.06.0157", "TCGA.06.0174", "TCGA.06.0219", "TCGA.06.0646", "TCGA.06.0750",
                                                             "TCGA.06.2557", "TCGA.06.2569", "TCGA.12.5299", "TCGA.14.0787", "TCGA.14.2554", "TCGA.15.0742", "TCGA.26.5134",
                                                             "TCGA.26.5136", "TCGA.26.5139", "TCGA.27.1832", "TCGA.27.1834", "TCGA.27.2521", "TCGA.27.2526", "TCGA.28.1753",
                                                             "TCGA.28.2499", "TCGA.28.5208", "TCGA.28.5209", "TCGA.32.1980", "TCGA.41.2571", "TCGA.41.4097", "TCGA.76.4928"))
  
  checkEquals(analysisDataSetup_nocov_sttr$unmatchedSamples, c("TCGA.06.0402", "TCGA.12.5301", "TCGA.02.0048", "TCGA.06.6699", "TCGA.14.1795", "TCGA.28.2502", "TCGA.76.6193",
                                                               "TCGA.41.3392", "TCGA.02.0033", "TCGA.19.1789", "TCGA.08.0345", "TCGA.81.5910", "TCGA.06.0413", "TCGA.28.1746",
                                                               "TCGA.14.1453", "TCGA.06.0140", "TCGA.14.0812", "TCGA.28.1755", "TCGA.08.0392", "TCGA.32.1977", "TCGA.06.0119",
                                                               "TCGA.08.0352", "TCGA.76.6192", "TCGA.76.6661", "TCGA.26.6174", "TCGA.32.1976", "TCGA.06.0169", "TCGA.14.3476",
                                                               "TCGA.06.6391", "TCGA.06.A5U1", "TCGA.28.1756", "TCGA.19.2621", "TCGA.06.A6S1", "TCGA.08.0385", "TCGA.12.0657",
                                                               "TCGA.14.1794", "TCGA.02.0439", "TCGA.76.4934", "TCGA.06.A6S0", "TCGA.06.0201", "TCGA.14.0867", "TCGA.06.0394",
                                                               "TCGA.06.0152", "TCGA.08.0529", "TCGA.06.6698", "TCGA.02.0107", "TCGA.06.0148", "TCGA.06.0179", "TCGA.12.0656",
                                                               "TCGA.76.6285", "TCGA.06.0182", "TCGA.02.0001", "TCGA.19.1788", "TCGA.27.2518", "TCGA.19.1389", "TCGA.14.1450",
                                                               "TCGA.02.0086", "TCGA.02.0325", "TCGA.06.0197", "TCGA.19.1386", "TCGA.02.0285", "TCGA.RR.A6KA", "TCGA.06.1086",
                                                               "TCGA.27.2527", "TCGA.08.0518", "TCGA.02.0006", "TCGA.26.6173", "TCGA.06.A7TL", "TCGA.12.1602", "TCGA.06.1805",
                                                               "TCGA.08.0355", "TCGA.32.1978", "TCGA.06.0150", "TCGA.19.5951", "TCGA.28.1750", "TCGA.32.4208", "TCGA.06.0162",
                                                               "TCGA.12.1098", "TCGA.26.1443", "TCGA.06.0881", "TCGA.02.0057", "TCGA.06.0877", "TCGA.41.2575", "TCGA.02.0084",
                                                               "TCGA.08.0517", "TCGA.19.0960", "TCGA.74.6581", "TCGA.02.0266", "TCGA.27.1838", "TCGA.02.0111", "TCGA.12.0826",
                                                               "TCGA.06.0155"))
  
  checkIdentical(analysisDataSetup_nocov_sttr$analysisData[, 1], 
                 c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1))
  
  checkIdentical(round(analysisDataSetup_nocov_sttr$analysisData[, 4], 2), c(0.20,  2.53, -0.50, -0.83,  2.32, -0.59, -0.42,  0.55,  0.94,  1.26, -0.89,  0.78, -0.53, -2.09, -1.94, -0.50, -0.49,
                                                                             -0.04, -0.02,  0.50, -1.44,  0.30,  0.85,  0.38,  0.16, -0.39, -1.06, -0.28))
  
  ################################
  #Original Data: One Covariate  #
  ################################
  set40.orig <- randomSample(object = ORIGdata,  survivalVar = "survival", nG1 = 40, nG2 = 40, cut = 0.5, all = FALSE, seed = 12345)
  
  analysisDataSetup_singlecov_orig <- analysisDataSetup(
    object = ORIGdata,
    sampleIDsG1 = set40.orig$shortSurvivors,
    sampleIDsG2 = set40.orig$longSurvivors,
    covariates = "ageAtDx",
    geneSet = "KANG_CISPLATIN_RESISTANCE_DN",
    sampleDescription ="TCGA GBM long vs. short survivors (single covariate)",
    geneSetDescription ="msgidb:KANG_CISPLATIN_RESISTANCE_DN")
  
  checkTrue(analysisDataSetup_singlecov_orig$sampleDescription == "TCGA GBM long vs. short survivors (single covariate)")
  checkTrue(analysisDataSetup_singlecov_orig$geneSetDescription == "msgidb:KANG_CISPLATIN_RESISTANCE_DN")
  checkEquals(dim(analysisDataSetup_singlecov_orig$analysisData), c(38, 8))
  
  checkTrue(length(intersect(analysisDataSetup_singlecov_orig$unmatchedSamples, analysisDataSetup_singlecov_orig$matchedSamples)) == 0)
  
  checkEquals(analysisDataSetup_singlecov_orig$matchedGenes, c("ASB8", "COL9A2", "GPR172A", "NBN", "RPS6KA4", "SEMA3C"))
  
  checkEquals(analysisDataSetup_singlecov_orig$unmatchedGenes, c("TMX4", "UBE2K"))
  
  checkEquals(analysisDataSetup_singlecov_orig$matchedSamples, c("TCGA.02.0001", "TCGA.02.0003", "TCGA.02.0006", "TCGA.02.0023", "TCGA.02.0028", "TCGA.02.0051", "TCGA.02.0060",
                                                                 "TCGA.02.0079", "TCGA.02.0106", "TCGA.02.0266", "TCGA.02.0326", "TCGA.02.0333", "TCGA.02.0339", "TCGA.02.0456",
                                                                 "TCGA.06.0124", "TCGA.06.0130", "TCGA.06.0141", "TCGA.06.0142", "TCGA.06.0155", "TCGA.06.0173", "TCGA.06.0194",
                                                                 "TCGA.06.0197", "TCGA.06.0219", "TCGA.06.0649", "TCGA.08.0345", "TCGA.08.0348", "TCGA.08.0349", "TCGA.08.0358",
                                                                 "TCGA.08.0385", "TCGA.08.0509", "TCGA.08.0510", "TCGA.08.0514", "TCGA.08.0516", "TCGA.12.0657", "TCGA.12.0707",
                                                                 "TCGA.14.0783", "TCGA.14.1396", "TCGA.14.1401"))
  
  checkEquals(analysisDataSetup_singlecov_orig$unmatchedSamples, c("TCGA.28.1747", "TCGA.76.6193", "TCGA.02.0071", "TCGA.76.4925", "TCGA.32.1980", "TCGA.14.1795", "TCGA.06.2559", "TCGA.14.1794", "TCGA.06.5418",
                                                                   "TCGA.19.1790", "TCGA.19.1789", "TCGA.14.3476", "TCGA.41.3392", "TCGA.19.A6J4", "TCGA.76.4934", "TCGA.19.2624", "TCGA.27.2526", "TCGA.14.0862",
                                                                   "TCGA.19.1788", "TCGA.12.5299", "TCGA.28.5218", "TCGA.06.A6S0", "TCGA.41.4097", "TCGA.06.5856", "TCGA.41.2572", "TCGA.76.6286", "TCGA.14.0866",
                                                                   "TCGA.19.1787", "TCGA.12.0662", "TCGA.06.1806", "TCGA.14.0790", "TCGA.06.2558", "TCGA.27.1835", "TCGA.12.1599", "TCGA.27.1834", "TCGA.12.0769",
                                                                   "TCGA.12.0653", "TCGA.76.4932", "TCGA.14.1823", "TCGA.14.1037", "TCGA.06.0165", "TCGA.06.5408"))
  
  checkIdentical(analysisDataSetup_singlecov_orig$analysisData[, 2], c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1))
  
  checkIdentical(round(analysisDataSetup_singlecov_orig$analysisData[, 5], 2), c(-0.99,  0.06, -0.04, -0.48, -1.20, -1.08,  0.30,  0.24, -0.02,  2.14,  1.19,  1.51, -0.06,  2.08,  0.94, -0.63, -0.02,
                                                                                 -1.37,  0.24, -0.39,  0.39,  0.34, -1.01, -0.97,  0.38,  0.32,  0.64,  1.76,  0.63,  0.09,  0.09,  0.07, -0.58,  0.15,
                                                                                 -0.62, -0.42,  1.32, -0.95))
  
  #################################
  #Original Data: Two Covariates  #
  #################################
  set40.orig <- randomSample(object = ORIGdata,  survivalVar = "survival", nG1 = 40, nG2 = 40, cut = 0.5, all = FALSE, seed = 12345)
  
  analysisDataSetup_twocov_orig <- analysisDataSetup(
    object = ORIGdata,
    sampleIDsG1 = set40.orig$shortSurvivors,
    sampleIDsG2 = set40.orig$longSurvivors,
    covariates = c("ageAtDx", "FirstProgression"),
    geneSet = "KANG_CISPLATIN_RESISTANCE_DN",
    sampleDescription ="TCGA GBM long vs. short survivors (two covariates)",
    geneSetDescription ="msgidb:KANG_CISPLATIN_RESISTANCE_DN")
  
  checkTrue(analysisDataSetup_twocov_orig$sampleDescription == "TCGA GBM long vs. short survivors (two covariates)")
  checkTrue(analysisDataSetup_twocov_orig$geneSetDescription == "msgidb:KANG_CISPLATIN_RESISTANCE_DN")
  checkEquals(dim(analysisDataSetup_twocov_orig$analysisData), c(38, 9))
  
  checkTrue(length(intersect(analysisDataSetup_twocov_orig$unmatchedSamples, analysisDataSetup_twocov_orig$matchedSamples)) == 0)
  
  checkEquals(analysisDataSetup_twocov_orig$matchedGenes, c("ASB8", "COL9A2", "GPR172A", "NBN", "RPS6KA4", "SEMA3C"))
  
  checkEquals(analysisDataSetup_twocov_orig$unmatchedGenes, c("TMX4", "UBE2K"))
  
  checkEquals(analysisDataSetup_twocov_orig$matchedSamples, c("TCGA.02.0001", "TCGA.02.0003", "TCGA.02.0006", "TCGA.02.0023", "TCGA.02.0028", "TCGA.02.0051", "TCGA.02.0060",
                                                              "TCGA.02.0079", "TCGA.02.0106", "TCGA.02.0266", "TCGA.02.0326", "TCGA.02.0333", "TCGA.02.0339", "TCGA.02.0456",
                                                              "TCGA.06.0124", "TCGA.06.0130", "TCGA.06.0141", "TCGA.06.0142", "TCGA.06.0155", "TCGA.06.0173", "TCGA.06.0194",
                                                              "TCGA.06.0197", "TCGA.06.0219", "TCGA.06.0649", "TCGA.08.0345", "TCGA.08.0348", "TCGA.08.0349", "TCGA.08.0358",
                                                              "TCGA.08.0385", "TCGA.08.0509", "TCGA.08.0510", "TCGA.08.0514", "TCGA.08.0516", "TCGA.12.0657", "TCGA.12.0707",
                                                              "TCGA.14.0783", "TCGA.14.1396", "TCGA.14.1401"))
  
  checkEquals(analysisDataSetup_twocov_orig$unmatchedSamples, c("TCGA.28.1747", "TCGA.76.6193", "TCGA.02.0071", "TCGA.76.4925", "TCGA.32.1980", "TCGA.14.1795", "TCGA.06.2559", "TCGA.14.1794", "TCGA.06.5418",
                                                                "TCGA.19.1790", "TCGA.19.1789", "TCGA.14.3476", "TCGA.41.3392", "TCGA.19.A6J4", "TCGA.76.4934", "TCGA.19.2624", "TCGA.27.2526", "TCGA.14.0862",
                                                                "TCGA.19.1788", "TCGA.12.5299", "TCGA.28.5218", "TCGA.06.A6S0", "TCGA.41.4097", "TCGA.06.5856", "TCGA.41.2572", "TCGA.76.6286", "TCGA.14.0866",
                                                                "TCGA.19.1787", "TCGA.12.0662", "TCGA.06.1806", "TCGA.14.0790", "TCGA.06.2558", "TCGA.27.1835", "TCGA.12.1599", "TCGA.27.1834", "TCGA.12.0769",
                                                                "TCGA.12.0653", "TCGA.76.4932", "TCGA.14.1823", "TCGA.14.1037", "TCGA.06.0165", "TCGA.06.5408"))
  
  
  checkIdentical(analysisDataSetup_twocov_orig$analysisData[, 3], c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1))
  
  checkIdentical(round(analysisDataSetup_twocov_orig$analysisData[, 4], 2), c(0.27,  1.11,  0.21, -0.31,  0.06, -0.01,  5.38, -1.14, -0.09, -1.68, -1.53,  0.29,  0.01,  0.11,  0.25, -0.80,  0.60,
                                                                              0.36, -1.13,  0.09, -0.05, -1.30, -0.95, -0.05,  0.61, -1.50,  0.05, -0.99,  1.44,  0.51, -0.45,  0.55,  0.78, -0.18,
                                                                              -0.18, -0.47,  0.65,  0.06))
  
  ###########################################################################################################
  #Original Data:  Only a single gene is shared between the specified gene set and the TCGA expression data #
  ###########################################################################################################
  set_20_30.orig <- randomSample(object = ORIGdata, survivalVar = "survival", nG1 = 20, nG2 = 30, cut = 0.5, all = FALSE, seed = -1684486789)
  
  analysisDataSetup_singleGene_orig <- analysisDataSetup(
    object = ORIGdata,
    sampleIDsG1 = set_20_30.orig$shortSurvivors,
    sampleIDsG2 = set_20_30.orig$longSurvivors,
    covariates = NULL,
    geneSet = "CAFFAREL_RESPONSE_TO_THC_8HR_3_UP",
    sampleDescription ="TCGA GBM long vs. short survivors (single gene)",
    geneSetDescription ="msgidb:CAFFAREL_RESPONSE_TO_THC_8HR_3_UP")

  checkEquals(analysisDataSetup_singleGene_orig$matchedGenes, c("RHOB"))
  
  checkEquals(analysisDataSetup_singleGene_orig$unmatchedGenes, c("PATL1", "AIMP2", "EPCAM", "CPEB4"))
  
  checkEquals(analysisDataSetup_singleGene_orig$matchedSamples, c("TCGA.02.0009", "TCGA.02.0028", "TCGA.02.0033", "TCGA.02.0060", "TCGA.02.0338", "TCGA.02.0439", "TCGA.06.0140",
                                                                  "TCGA.06.0158", "TCGA.06.0173", "TCGA.06.0189", "TCGA.06.0195", "TCGA.06.0402", "TCGA.06.0648", "TCGA.08.0350",
                                                                  "TCGA.08.0359", "TCGA.08.0373", "TCGA.08.0380", "TCGA.08.0390", "TCGA.08.0392", "TCGA.08.0524", "TCGA.08.0529",
                                                                  "TCGA.12.0707", "TCGA.12.0773", "TCGA.12.0778", "TCGA.12.0822", "TCGA.12.1090", "TCGA.12.1096", "TCGA.14.1454",
                                                                  "TCGA.16.1062", "TCGA.16.1063"))
  
  checkEquals(analysisDataSetup_singleGene_orig$unmatchedSamples, c("TCGA.06.5410", "TCGA.76.6656", "TCGA.19.1790", "TCGA.06.0119", "TCGA.76.4929", "TCGA.06.2557", "TCGA.14.3476",
                                                                    "TCGA.14.1043", "TCGA.06.A6S1", "TCGA.28.6450", "TCGA.14.0812", "TCGA.12.3650", "TCGA.32.1986", "TCGA.32.1973",
                                                                    "TCGA.32.2615", "TCGA.06.1806", "TCGA.27.1834", "TCGA.14.2554", "TCGA.14.0740", "TCGA.19.0957"))
  
  checkIdentical(analysisDataSetup_singleGene_orig$analysisData[, 1], c(0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  
  
  checkIdentical(round(analysisDataSetup_singleGene_orig$analysisData[, 2], 2), c(0.81,  0.14,  0.16,  1.41, -3.03,  0.59,  0.23,  0.74,  0.18, 0.07, -1.60,  0.17, -0.05, -0.33,  1.20,  0.11, -0.96,
                                                                                  -1.80, -3.04, -1.41, -0.22,  0.94, -0.10,  1.06,  0.69, -3.20,  0.41, -1.33, -0.13, -0.65))
  
  ###########################################################################################
  #Original Data:  No shared genes between specified gene set and the TCGA expression data  #
  ###########################################################################################
  set_20_30.orig <- randomSample(object = ORIGdata, survivalVar = "survival", nG1 = 20, nG2 = 30, cut = 0.5, all = FALSE, seed = -1684486789)
  
  obj_chk.orig = ORIGdata
  
  checkEquals(which(names(getGeneSets(obj_chk.orig)) == "chr4p11"), 3315)
  checkEquals(length(which(colnames(getExpressionData(obj_chk.orig)) %in% getGeneSets(obj_chk.orig)["chr4p11"][[1]])), 0)
  
  analysisDataSetup_noGene_orig <- try(analysisDataSetup(
    obj = ORIGdata,
    sampleIDsG1 = set_20_30.orig$shortSurvivors,
    sampleIDsG2 = set_20_30.orig$longSurvivors,
    covariates = NULL,
    geneSet = "chr4p11",
    sampleDescription ="TCGA GBM long vs. short survivors (no shared genes)",
    geneSetDescription ="msgidb:chr4p11"), silent = TRUE)
  
  checkEquals(analysisDataSetup_noGene_orig[[1]], "Error in .local(object, ...) : \n  None of the genes in the specified gene sets are contained in the gene expression data.\n")
} #test_analysisDataSetup

#------------------------------------------------------------------------------------------------------------------------
test_geneSetScoreTest = function () 
{  
  printf("--- test_geneSetScoreTest")
  
  data(ORIGdata)
  data(STTRdata)
  
  #############################################################################
  #Score test with a subset of the samples and no covariates (original data)  #
  #############################################################################
  set40.orig <- randomSample(object = ORIGdata,  survivalVar = "survival", nG1 = 40, nG2 = 40, cut = 0.5, all = FALSE, seed = 12345)
  
  skat_nocov_orig <- geneSetScoreTest(
    analysisDataSetup(object = ORIGdata,
    sampleIDsG1 = set40.orig$shortSurvivors,
    sampleIDsG2 = set40.orig$longSurvivors,
    geneSet = "KANG_CISPLATIN_RESISTANCE_DN",
    covariates = NULL,
    sampleDescription ="TCGA GBM long vs. short survivors",
    geneSetDescription ="msgidb:KANG_CISPLATIN_RESISTANCE_DN"))
  
  checkEquals(skat_nocov_orig$null.model, "group ~ 1")
  checkTrue(is.null(skat_nocov_orig$null.warnings))
  checkEquals(skat_nocov_orig$alternative.model, "group ~ 1 + ASB8 + COL9A2 + GPR172A + NBN + RPS6KA4 + SEMA3C")
  checkEquals(round(skat_nocov_orig$p.value, 4), 0.3837)
  
  chkSummary <- summary(skat_nocov_orig)
  checkEquals(chkSummary, "Null Model:  group ~ 1 \nAlternative Model:  group ~ 1 + ASB8 + COL9A2 + GPR172A + NBN + RPS6KA4 + SEMA3C \n\nP-value: 0.3837")

  #############################################################################
  #Score test with a subset of the samples and no covariates (STTR data)      #
  #############################################################################
  set60.sttr <- randomSample(object = STTRdata,  survivalVar = "Survival", nG1 = 60, nG2 = 60, cut = 100, all = FALSE, seed = 12345)
  
  skat_nocov_sttr <- geneSetScoreTest(
    analysisDataSetup(object = STTRdata,
    sampleIDsG1 = set60.sttr$shortSurvivors,
    sampleIDsG2 = set60.sttr$longSurvivors,
    geneSet = "geneset9",
    covariates = NULL,
    sampleDescription ="long vs. short survivors",
    geneSetDescription ="geneset9"))
  
  checkEquals(skat_nocov_sttr$null.model, "group ~ 1")
  checkTrue(is.null(skat_nocov_sttr$null.warnings))
  checkEquals(skat_nocov_sttr$alternative.model, "group ~ 1 + LINC00575 + TMEM244 + CRISP1 + HCP5 + NAIP + PMPCB + SALL4 + TMPRSS3")
  checkEquals(round(skat_nocov_sttr$p.value, 4), 0.1907)

  chkSummary <- summary(skat_nocov_sttr)
  checkEquals(chkSummary, "Null Model:  group ~ 1 \nAlternative Model:  group ~ 1 + LINC00575 + TMEM244 + CRISP1 + HCP5 + NAIP + PMPCB + SALL4 + TMPRSS3 \n\nP-value: 0.1907")

  ########################################################################
  #Score test with all of the samples and no covariates (original data)  #
  ########################################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = TRUE)
  
  skatAll_orig <- geneSetScoreTest(
    analysisDataSetup(object = ORIGdata,
    sampleIDsG1 = setAll.orig$shortSurvivors,
    sampleIDsG2 = setAll.orig$longSurvivors,
    covariates = NULL,
    geneSet = "KANG_CISPLATIN_RESISTANCE_DN",
    sampleDescription ="TCGA GBM long vs. short survivors",
    geneSetDescription ="msgidb:KANG_CISPLATIN_RESISTANCE_DN"))
  
  checkEquals(skatAll_orig$null.model, "group ~ 1")
  checkTrue(is.null(skatAll_orig$null.warnings))
  checkEquals(skatAll_orig$alternative.model, "group ~ 1 + ASB8 + COL9A2 + GPR172A + NBN + RPS6KA4 + SEMA3C")
  checkEquals(round(skatAll_orig$p.value, 4), 0.3711)
  
  chkSummary <- summary(skatAll_orig)
  checkEquals(chkSummary, "Null Model:  group ~ 1 \nAlternative Model:  group ~ 1 + ASB8 + COL9A2 + GPR172A + NBN + RPS6KA4 + SEMA3C \n\nP-value: 0.3711")

  ##########################################################################################
  #Score test with all of the samples and age at diagnosis as a covariate (original data)  #
  ##########################################################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = TRUE)
  
  skatAll_onecov_orig <- geneSetScoreTest(
    analysisDataSetup(object = ORIGdata,
    sampleIDsG1 = setAll.orig$shortSurvivors,
    sampleIDsG2 = setAll.orig$longSurvivors,
    covariates = "ageAtDx",
    geneSet = "KANG_CISPLATIN_RESISTANCE_DN",
    sampleDescription ="TCGA GBM long vs. short survivors",
    geneSetDescription ="msgidb:KANG_CISPLATIN_RESISTANCE_DN"))
  
  checkEquals(skatAll_onecov_orig$null.model, "group ~ 1 + ageAtDx")
  checkTrue(is.null(skatAll_onecov_orig$null.warnings))
  checkEquals(skatAll_onecov_orig$alternative.model, "group ~ 1 + ageAtDx + ASB8 + COL9A2 + GPR172A + NBN + RPS6KA4 + SEMA3C")
  checkEquals(round(skatAll_onecov_orig$p.value, 4), 0.5523)
  
  chkSummary <- summary(skatAll_onecov_orig)
  checkEquals(chkSummary, "Null Model:  group ~ 1 + ageAtDx \nAlternative Model:  group ~ 1 + ageAtDx + ASB8 + COL9A2 + GPR172A + NBN + RPS6KA4 + SEMA3C \n\nP-value: 0.5523")
  
  #########################################################################################
  #Score test with all of the samples, three covariates and two warnings (original data)  #
  #########################################################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = TRUE)
  
  skatAll_threecov_warn_orig <- geneSetScoreTest(
    analysisDataSetup(object = ORIGdata,
    sampleIDsG1 = setAll.orig$shortSurvivors,
    sampleIDsG2 = setAll.orig$longSurvivors,
    covariates = c("ageAtDx", "FirstProgression", "ChemoAgent"),
    geneSet = "ORGANELLAR_SMALL_RIBOSOMAL_SUBUNIT",
    sampleDescription ="TCGA GBM long vs. short survivors",
    geneSetDescription ="msgidb:ORGANELLAR_SMALL_RIBOSOMAL_SUBUNIT"))
  
  checkEquals(as.character(skatAll_threecov_warn_orig$null.warnings), 
              c("120  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",
                "glm.fit: fitted probabilities numerically 0 or 1 occurred"))
  
  checkEquals(skatAll_threecov_warn_orig$null.model, "group ~ 1 + ageAtDx + FirstProgression + ChemoAgent")
  checkEquals(skatAll_threecov_warn_orig$alternative.model, "group ~ 1 + ageAtDx + FirstProgression + ChemoAgent + MRPL42 + MRPS11 + MRPS15 + MRPS16 + MRPS18A + MRPS22 + MRPS35")
  checkEquals(round(skatAll_threecov_warn_orig$p.value, 4), 0.2252)
  
  chkSummary <- summary(skatAll_threecov_warn_orig)
  checkEquals(chkSummary, "Null Model:  group ~ 1 + ageAtDx + FirstProgression + ChemoAgent \nAlternative Model:  group ~ 1 + ageAtDx + FirstProgression + ChemoAgent + MRPL42 + MRPS11 + MRPS15 + MRPS16 + MRPS18A + MRPS22 + MRPS35 \n\nP-value: 0.2252 \n\nNull Model Warnings: 120  samples have either missing phenotype or missing covariates. They are excluded from the analysis!; glm.fit: fitted probabilities numerically 0 or 1 occurred")
  
  #########################################################################################
  #Score test with all of the samples, three covariates and two warnings (STTR data)      #
  #########################################################################################
  setAll.sttr <- randomSample(object = STTRdata, survivalVar = "Survival", cut = 100, all = TRUE)

  skatAll_threecov_warn_sttr <- geneSetScoreTest(
    analysisDataSetup(object = STTRdata,
    sampleIDsG1 = setAll.sttr$shortSurvivors,
    sampleIDsG2 = setAll.sttr$longSurvivors,
    covariates = c("AgeDx", "TimeFirstProgression", "Birth.gender"),
    geneSet = "geneset8",
    sampleDescription ="long vs. short survivors",
    geneSetDescription ="geneset8"))
  
  checkEquals(as.character(skatAll_threecov_warn_sttr$null.warnings), 
              c("68  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",
                "glm.fit: fitted probabilities numerically 0 or 1 occurred"))
  
  checkEquals(skatAll_threecov_warn_sttr$null.model, "group ~ 1 + AgeDx + TimeFirstProgression + Birth.gender")
  checkEquals(skatAll_threecov_warn_sttr$alternative.model, "group ~ 1 + AgeDx + TimeFirstProgression + Birth.gender + ARHGAP42 + KIFC1 + TAF13")
  checkEquals(round(skatAll_threecov_warn_sttr$p.value, 4), 0.0904)
  
  chkSummary <- summary(skatAll_threecov_warn_sttr)
  checkEquals(chkSummary, "Null Model:  group ~ 1 + AgeDx + TimeFirstProgression + Birth.gender \nAlternative Model:  group ~ 1 + AgeDx + TimeFirstProgression + Birth.gender + ARHGAP42 + KIFC1 + TAF13 \n\nP-value: 0.0904 \n\nNull Model Warnings: 68  samples have either missing phenotype or missing covariates. They are excluded from the analysis!; glm.fit: fitted probabilities numerically 0 or 1 occurred")
    
  ###################################################################################
  #Score test with a subset of the samples and only a single gene (original data)   #
  ###################################################################################
  set40.orig <- randomSample(object = ORIGdata,  survivalVar = "survival", nG1 = 40, nG2 = 40, cut = 0.5, all = FALSE, seed = 12345)
  
  skat_singleGene_orig <- geneSetScoreTest(
    analysisDataSetup(object = ORIGdata,
    sampleIDsG1 = set40.orig$shortSurvivors,
    sampleIDsG2 = set40.orig$longSurvivors,
    covariates = NULL,
    geneSet = "CAFFAREL_RESPONSE_TO_THC_8HR_3_UP",
    sampleDescription ="TCGA GBM long vs. short survivors",
    geneSetDescription ="msgidb:CAFFAREL_RESPONSE_TO_THC_8HR_3_UP"))
  
  checkEquals(skat_singleGene_orig$null.model, "group ~ 1")
  checkTrue(is.null(skat_singleGene_orig$null.warnings))
  checkEquals(skat_singleGene_orig$alternative.model, "group ~ 1 + RHOB")
  checkEquals(round(skat_singleGene_orig$p.value, 4), 0.2347)
  
  chkSummary <- summary(skat_singleGene_orig)
  checkEquals(chkSummary, "Null Model:  group ~ 1 \nAlternative Model:  group ~ 1 + RHOB \n\nP-value: 0.2347")
  
  ###########################################################################
  #Score test with missing values in the expression data  (original data)   #
  ###########################################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = TRUE)
  
  genes <- getGeneSets(ORIGdata)[["GSE1460_INTRATHYMIC_T_PROGENITOR_VS_DP_THYMOCYTE_UP"]][
    which(getGeneSets(ORIGdata)[["GSE1460_INTRATHYMIC_T_PROGENITOR_VS_DP_THYMOCYTE_UP"]] %in% colnames(getExpressionData(ORIGdata)))]
  
  missings <- apply(getExpressionData(ORIGdata)[, genes], 2, function(x) all(!is.na(x)))
  
  checkEquals(names(which(missings == FALSE)), "S100G")
  
  skat_missingGene_orig <- geneSetScoreTest(
    analysisDataSetup(object = ORIGdata,
    sampleIDsG1 = setAll.orig$shortSurvivors,
    sampleIDsG2 = setAll.orig$longSurvivors,
    covariates = c("ageAtDx", "FirstProgression", "ChemoAgent"),
    geneSet = "GSE1460_INTRATHYMIC_T_PROGENITOR_VS_DP_THYMOCYTE_UP",
    sampleDescription ="TCGA GBM long vs. short survivors",
    geneSetDescription ="msgidb:GSE1460_INTRATHYMIC_T_PROGENITOR_VS_DP_THYMOCYTE_UP"))
  
  checkEquals(as.character(skat_missingGene_orig$null.warnings), 
              c("118  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",
                "glm.fit: fitted probabilities numerically 0 or 1 occurred"))
  
  checkEquals(skat_missingGene_orig$null.model, "group ~ 1 + ageAtDx + FirstProgression + ChemoAgent")
  checkEquals(skat_missingGene_orig$alternative.model, "group ~ 1 + ageAtDx + FirstProgression + ChemoAgent + ABCC9 + ABCE1 + ABHD14A + ADAMTS1 + ADRB3 + ALCAM + ANXA4 + APEX1 + AQP9 + ARMCX1 + ATIC + ATP10A + BAZ2B + BMP15 + BMX + BRIP1 + CA6 + CASP1 + CCDC102B + CD55 + CD74 + CDH8 + CDKAL1 + CFH + CIDEB + CKAP4 + CLCA4 + CLIP3 + CREB3L2 + CRTAP + CTNNA3 + CUL1 + DHODH + DSC3 + DTX4 + DUSP22 + EEF1G + EIF2B3 + EIF4EBP1 + ELF2 + ELN + EMD + ERAL1 + ERLIN1 + ETF1 + ETFB + EXOSC5 + FAIM + FAM49A + FXYD2 + FZD6 + GAS6 + GNG11 + GNL2 + GPD1L + GPRC5A + GPX1 + GPX7 + GRB14 + GSTM3 + GUCY1A3 + GUCY1B3 + HDGFRP3 + HES1 + HLA-DMA + IFI44 + IFI44L + IGLL1 + IL19 + IMPACT + IMPDH1 + INTS7 + KCNK10 + LAMC1 + LARS2 + LDB1 + LILRA4 + LUZP2 + MALL + MAP9 + MCTP2 + MDFIC + MEP1B + MEST + METTL8 + MFAP4 + MINA + MRE11A + MRPL11 + MRTO4 + MS4A5 + MSRB2 + MTCP1 + MYOZ2 + NBN + NEK11 + NOS1 + NOTCH1 + NUBPL + OBFC1 + PAICS + PBX1 + PCCB + PCNA + PDE3B + PDZRN3 + PFKM + PIM1 + PKD2 + PKIG + PLA2G4C + PLAC8 + POFUT1 + POLR3B + PON2 + PPAP2B + PPEF2 + PRDX4 + PRKCA + PRMT5 + PTGER3 + PTPLA + RAB40B + RNF130 + RNF144A + RPS8 + RRAS2 + RRP1 + RUVBL1 + S100G + SAMM50 + SAMSN1 + SCN3A + SCN9A + SEC31B + SELL + SH3BGR + SLC15A3 + SND1 + SOCS2 + STK3 + STOM + STOML2 + TBCE + TCF7L2 + TFAM + TFRC + TIMP1 + TIMP4 + TLR7 + TMED3 + TNFSF4 + TRIM29 + TUBA1A + TUFT1 + TYRP1 + UCHL5 + UCK2 + ULK4 + VIM + WDR3 + WDR41 + XKR8 + XPOT + XYLT1 + ZCCHC14 + ZMYND10 + ZNF132 + ZNF282 + ZNF550 + ZNF643")
  checkEquals(round(skat_missingGene_orig$p.value, 4), 0.1886)
  
  chkSummary <- summary(skat_missingGene_orig)
  checkEquals(chkSummary, "Null Model:  group ~ 1 + ageAtDx + FirstProgression + ChemoAgent \nAlternative Model:  group ~ 1 + ageAtDx + FirstProgression + ChemoAgent + ABCC9 + ABCE1 + ABHD14A + ADAMTS1 + ADRB3 + ALCAM + ANXA4 + APEX1 + AQP9 + ARMCX1 + ATIC + ATP10A + BAZ2B + BMP15 + BMX + BRIP1 + CA6 + CASP1 + CCDC102B + CD55 + CD74 + CDH8 + CDKAL1 + CFH + CIDEB + CKAP4 + CLCA4 + CLIP3 + CREB3L2 + CRTAP + CTNNA3 + CUL1 + DHODH + DSC3 + DTX4 + DUSP22 + EEF1G + EIF2B3 + EIF4EBP1 + ELF2 + ELN + EMD + ERAL1 + ERLIN1 + ETF1 + ETFB + EXOSC5 + FAIM + FAM49A + FXYD2 + FZD6 + GAS6 + GNG11 + GNL2 + GPD1L + GPRC5A + GPX1 + GPX7 + GRB14 + GSTM3 + GUCY1A3 + GUCY1B3 + HDGFRP3 + HES1 + HLA-DMA + IFI44 + IFI44L + IGLL1 + IL19 + IMPACT + IMPDH1 + INTS7 + KCNK10 + LAMC1 + LARS2 + LDB1 + LILRA4 + LUZP2 + MALL + MAP9 + MCTP2 + MDFIC + MEP1B + MEST + METTL8 + MFAP4 + MINA + MRE11A + MRPL11 + MRTO4 + MS4A5 + MSRB2 + MTCP1 + MYOZ2 + NBN + NEK11 + NOS1 + NOTCH1 + NUBPL + OBFC1 + PAICS + PBX1 + PCCB + PCNA + PDE3B + PDZRN3 + PFKM + PIM1 + PKD2 + PKIG + PLA2G4C + PLAC8 + POFUT1 + POLR3B + PON2 + PPAP2B + PPEF2 + PRDX4 + PRKCA + PRMT5 + PTGER3 + PTPLA + RAB40B + RNF130 + RNF144A + RPS8 + RRAS2 + RRP1 + RUVBL1 + S100G + SAMM50 + SAMSN1 + SCN3A + SCN9A + SEC31B + SELL + SH3BGR + SLC15A3 + SND1 + SOCS2 + STK3 + STOM + STOML2 + TBCE + TCF7L2 + TFAM + TFRC + TIMP1 + TIMP4 + TLR7 + TMED3 + TNFSF4 + TRIM29 + TUBA1A + TUFT1 + TYRP1 + UCHL5 + UCK2 + ULK4 + VIM + WDR3 + WDR41 + XKR8 + XPOT + XYLT1 + ZCCHC14 + ZMYND10 + ZNF132 + ZNF282 + ZNF550 + ZNF643 \n\nP-value: 0.1886 \n\nNull Model Warnings: 118  samples have either missing phenotype or missing covariates. They are excluded from the analysis!; glm.fit: fitted probabilities numerically 0 or 1 occurred")
  
  ################################################################################
  #Random gene sets with all of the samples and no covariates  (original data)   #
  ################################################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = TRUE)
  all.geneSets <- getGeneSets(ORIGdata)
  
  set.seed(-1327681546)
  samp.geneSets <- sample(all.geneSets, 10)
  
  run.start.genesets <- proc.time()
  skatRes <- foreach(i = 1:length(samp.geneSets), .errorhandling="pass") %do% {  
    res <- geneSetScoreTest(
      analysisDataSetup(object = ORIGdata,
      sampleIDsG1 = setAll.orig$shortSurvivors,
      sampleIDsG2 = setAll.orig$longSurvivors,
      covariates = NULL,
      geneSet = names(samp.geneSets)[i],
      sampleDescription = "TCGA GBM long vs. short survivors (random gene set, all samples)",
      geneSetDescription = names(samp.geneSets)[i]))
  }
  names(skatRes) <- names(samp.geneSets)
  run.time.genesets <- proc.time() - run.start.genesets  
  
  checkEquals(names(skatRes), names(samp.geneSets))
  
  checkEquals(length(skatRes), 10)
  checkEquals(any(!unlist(lapply(skatRes, class)) == "geneSetScoreTest"), FALSE)
  
  checkEquals(skatRes[[1]]$null.model, "group ~ 1")
  checkTrue(is.null(skatRes[[1]]$null.warnings))
  checkEquals(skatRes[[1]]$alternative.model, "group ~ 1 + CASP10 + CDH2 + IGFBP2 + SEPT7 + TFG")
  checkEquals(round(skatRes[[1]]$p.value, 4), 0.2439)
  
  checkEquals(skatRes[[2]]$null.model, "group ~ 1")
  checkTrue(is.null(skatRes[[2]]$null.warnings))
  checkEquals(skatRes[[2]]$alternative.model, "group ~ 1 + AMFR + ATP6V0A2 + ATRX + BAHD1 + BPHL + BRCA1 + BTD + CAMK2G + CCNF + CETN3 + CHD3 + CLP1 + CLPX + CPSF4 + CSTF3 + CYP4F12 + EIF5B + EP400 + ERAL1 + ERCC2 + EXTL3 + FANCG + FANCI + FDXR + FNTB + FOXD1 + FRYL + GALNT2 + GRIK5 + GRIP2 + GSK3B + HNRNPL + HTR7 + ILVBL + IMPA1 + INPP5E + JRK + KIAA0195 + KIAA0586 + KLHL18 + KRT33A + MC2R + MEA1 + MFN1 + MSH3 + MTX1 + MUTYH + NFRKB + NFYB + NKRF + NMT1 + NR2C1 + PAXIP1 + PCGF1 + PDXDC1 + PFDN6 + PHF10 + PIGB + PIGF + PIK3R2 + PLEKHB1 + POP4 + PPP5C + RAP1A + RBBP8 + REV3L + RFC5 + RPS6KB2 + SCAMP1 + SH2B1 + SLC24A1 + SLC25A11 + SLC30A3 + SPAST + SPRED2 + SSR1 + SYNJ2 + TAF2 + TMEM11 + TRIM27 + UBE4B + WDR62 + ZNF500 + ZNF592")
  checkEquals(round(skatRes[[2]]$p.value, 4), 0.745)
  
  checkEquals(skatRes[[10]]$null.model, "group ~ 1")
  checkTrue(is.null(skatRes[[10]]$null.warnings))
  checkEquals(skatRes[[10]]$alternative.model, "group ~ 1 + ABCC9 + ABCE1 + ABHD14A + ADAMTS1 + ADRB3 + ALCAM + ANXA4 + APEX1 + AQP9 + ARMCX1 + ATIC + ATP10A + BAZ2B + BMP15 + BMX + BRIP1 + CA6 + CASP1 + CCDC102B + CD55 + CD74 + CDH8 + CDKAL1 + CFH + CIDEB + CKAP4 + CLCA4 + CLIP3 + CREB3L2 + CRTAP + CTNNA3 + CUL1 + DHODH + DSC3 + DTX4 + DUSP22 + EEF1G + EIF2B3 + EIF4EBP1 + ELF2 + ELN + EMD + ERAL1 + ERLIN1 + ETF1 + ETFB + EXOSC5 + FAIM + FAM49A + FXYD2 + FZD6 + GAS6 + GNG11 + GNL2 + GPD1L + GPRC5A + GPX1 + GPX7 + GRB14 + GSTM3 + GUCY1A3 + GUCY1B3 + HDGFRP3 + HES1 + HLA-DMA + IFI44 + IFI44L + IGLL1 + IL19 + IMPACT + IMPDH1 + INTS7 + KCNK10 + LAMC1 + LARS2 + LDB1 + LILRA4 + LUZP2 + MALL + MAP9 + MCTP2 + MDFIC + MEP1B + MEST + METTL8 + MFAP4 + MINA + MRE11A + MRPL11 + MRTO4 + MS4A5 + MSRB2 + MTCP1 + MYOZ2 + NBN + NEK11 + NOS1 + NOTCH1 + NUBPL + OBFC1 + PAICS + PBX1 + PCCB + PCNA + PDE3B + PDZRN3 + PFKM + PIM1 + PKD2 + PKIG + PLA2G4C + PLAC8 + POFUT1 + POLR3B + PON2 + PPAP2B + PPEF2 + PRDX4 + PRKCA + PRMT5 + PTGER3 + PTPLA + RAB40B + RNF130 + RNF144A + RPS8 + RRAS2 + RRP1 + RUVBL1 + S100G + SAMM50 + SAMSN1 + SCN3A + SCN9A + SEC31B + SELL + SH3BGR + SLC15A3 + SND1 + SOCS2 + STK3 + STOM + STOML2 + TBCE + TCF7L2 + TFAM + TFRC + TIMP1 + TIMP4 + TLR7 + TMED3 + TNFSF4 + TRIM29 + TUBA1A + TUFT1 + TYRP1 + UCHL5 + UCK2 + ULK4 + VIM + WDR3 + WDR41 + XKR8 + XPOT + XYLT1 + ZCCHC14 + ZMYND10 + ZNF132 + ZNF282 + ZNF550 + ZNF643")
  checkEquals(round(skatRes[[10]]$p.value, 4), 0.6569) 
  
  ###########################################################################
  #Random gene sets with all of the samples and covariate (original data)   #
  ###########################################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = TRUE)
  all.geneSets <- getGeneSets(ORIGdata)
  
  set.seed(-1327681546)
  samp.geneSets <- sample(all.geneSets, 10)
  
  res <- list()
  run.start.genesets <- proc.time()
  skatRes <- foreach(i = 1:length(samp.geneSets), .errorhandling="pass") %do% {  
    res[[i]] <- geneSetScoreTest(
      analysisDataSetup(object = ORIGdata,
      sampleIDsG1 = setAll.orig$shortSurvivors,
      sampleIDsG2 = setAll.orig$longSurvivors,
      covariates = "ageAtDx",
      geneSet = names(samp.geneSets)[i],
      sampleDescription = "TCGA GBM long vs. short survivors (random gene set, all samples)",
      geneSetDescription = names(samp.geneSets)[i]) )
  }
  names(skatRes) <- names(samp.geneSets)
  run.time.genesets <- proc.time() - run.start.genesets  
  
  checkEquals(names(skatRes), names(samp.geneSets))
  
  checkEquals(length(skatRes), 10)
  checkEquals(any(!unlist(lapply(skatRes, class)) == "geneSetScoreTest"), FALSE)
  
  checkEquals(skatRes[[1]]$null.model, "group ~ 1 + ageAtDx")
  checkTrue(is.null(skatRes[[1]]$null.warnings))
  checkEquals(skatRes[[1]]$alternative.model, "group ~ 1 + ageAtDx + CASP10 + CDH2 + IGFBP2 + SEPT7 + TFG")
  checkEquals(round(skatRes[[1]]$p.value, 4), 0.3061)
  
  checkEquals(skatRes[[2]]$null.model, "group ~ 1 + ageAtDx")
  checkTrue(is.null(skatRes[[2]]$null.warnings))
  checkEquals(skatRes[[2]]$alternative.model, "group ~ 1 + ageAtDx + AMFR + ATP6V0A2 + ATRX + BAHD1 + BPHL + BRCA1 + BTD + CAMK2G + CCNF + CETN3 + CHD3 + CLP1 + CLPX + CPSF4 + CSTF3 + CYP4F12 + EIF5B + EP400 + ERAL1 + ERCC2 + EXTL3 + FANCG + FANCI + FDXR + FNTB + FOXD1 + FRYL + GALNT2 + GRIK5 + GRIP2 + GSK3B + HNRNPL + HTR7 + ILVBL + IMPA1 + INPP5E + JRK + KIAA0195 + KIAA0586 + KLHL18 + KRT33A + MC2R + MEA1 + MFN1 + MSH3 + MTX1 + MUTYH + NFRKB + NFYB + NKRF + NMT1 + NR2C1 + PAXIP1 + PCGF1 + PDXDC1 + PFDN6 + PHF10 + PIGB + PIGF + PIK3R2 + PLEKHB1 + POP4 + PPP5C + RAP1A + RBBP8 + REV3L + RFC5 + RPS6KB2 + SCAMP1 + SH2B1 + SLC24A1 + SLC25A11 + SLC30A3 + SPAST + SPRED2 + SSR1 + SYNJ2 + TAF2 + TMEM11 + TRIM27 + UBE4B + WDR62 + ZNF500 + ZNF592")
  checkEquals(round(skatRes[[2]]$p.value, 4), 0.5352)
  
  checkEquals(skatRes[[10]]$null.model, "group ~ 1 + ageAtDx")
  checkTrue(is.null(skatRes[[10]]$null.warnings))
  checkEquals(skatRes[[10]]$alternative.model, "group ~ 1 + ageAtDx + ABCC9 + ABCE1 + ABHD14A + ADAMTS1 + ADRB3 + ALCAM + ANXA4 + APEX1 + AQP9 + ARMCX1 + ATIC + ATP10A + BAZ2B + BMP15 + BMX + BRIP1 + CA6 + CASP1 + CCDC102B + CD55 + CD74 + CDH8 + CDKAL1 + CFH + CIDEB + CKAP4 + CLCA4 + CLIP3 + CREB3L2 + CRTAP + CTNNA3 + CUL1 + DHODH + DSC3 + DTX4 + DUSP22 + EEF1G + EIF2B3 + EIF4EBP1 + ELF2 + ELN + EMD + ERAL1 + ERLIN1 + ETF1 + ETFB + EXOSC5 + FAIM + FAM49A + FXYD2 + FZD6 + GAS6 + GNG11 + GNL2 + GPD1L + GPRC5A + GPX1 + GPX7 + GRB14 + GSTM3 + GUCY1A3 + GUCY1B3 + HDGFRP3 + HES1 + HLA-DMA + IFI44 + IFI44L + IGLL1 + IL19 + IMPACT + IMPDH1 + INTS7 + KCNK10 + LAMC1 + LARS2 + LDB1 + LILRA4 + LUZP2 + MALL + MAP9 + MCTP2 + MDFIC + MEP1B + MEST + METTL8 + MFAP4 + MINA + MRE11A + MRPL11 + MRTO4 + MS4A5 + MSRB2 + MTCP1 + MYOZ2 + NBN + NEK11 + NOS1 + NOTCH1 + NUBPL + OBFC1 + PAICS + PBX1 + PCCB + PCNA + PDE3B + PDZRN3 + PFKM + PIM1 + PKD2 + PKIG + PLA2G4C + PLAC8 + POFUT1 + POLR3B + PON2 + PPAP2B + PPEF2 + PRDX4 + PRKCA + PRMT5 + PTGER3 + PTPLA + RAB40B + RNF130 + RNF144A + RPS8 + RRAS2 + RRP1 + RUVBL1 + S100G + SAMM50 + SAMSN1 + SCN3A + SCN9A + SEC31B + SELL + SH3BGR + SLC15A3 + SND1 + SOCS2 + STK3 + STOM + STOML2 + TBCE + TCF7L2 + TFAM + TFRC + TIMP1 + TIMP4 + TLR7 + TMED3 + TNFSF4 + TRIM29 + TUBA1A + TUFT1 + TYRP1 + UCHL5 + UCK2 + ULK4 + VIM + WDR3 + WDR41 + XKR8 + XPOT + XYLT1 + ZCCHC14 + ZMYND10 + ZNF132 + ZNF282 + ZNF550 + ZNF643")
  checkEquals(round(skatRes[[10]]$p.value, 4), 0.8078)
  
  ################################################################################
  #Random gene sets with all of the samples and two covariates (original data)   #
  ################################################################################
  setAll.orig <- randomSample(object = ORIGdata, survivalVar = "survival", cut = 0.5, all = TRUE)
  all.geneSets <- getGeneSets(ORIGdata)
  
  set.seed(-1327681546)
  samp.geneSets <- sample(all.geneSets, 10)
  
  res <- list()
  run.start.genesets <- proc.time()
  skatRes <- foreach(i = 1:length(samp.geneSets), .errorhandling="pass") %do% {  
    res[[i]] <- geneSetScoreTest(
      analysisDataSetup(obj = ORIGdata,
      sampleIDsG1 = setAll.orig$shortSurvivors,
      sampleIDsG2 = setAll.orig$longSurvivors,
      covariates = c("ageAtDx", "ChemoAgent"),
      geneSet = names(samp.geneSets)[i],
      sampleDescription = "TCGA GBM long vs. short survivors (random gene set, all samples)",
      geneSetDescription = names(samp.geneSets)[i]))
  }
  names(skatRes) <- names(samp.geneSets)
  run.time.genesets <- proc.time() - run.start.genesets  
  
  checkEquals(names(skatRes), names(samp.geneSets))
  
  checkEquals(length(skatRes), 10)
  checkEquals(any(!unlist(lapply(skatRes, class)) == "geneSetScoreTest"), FALSE)
  
  checkEquals(skatRes[[1]]$null.model, "group ~ 1 + ageAtDx + ChemoAgent")
  checkTrue(!is.null(skatRes[[1]]$null.warnings))
  checkEquals(skatRes[[1]]$alternative.model, "group ~ 1 + ageAtDx + ChemoAgent + CASP10 + CDH2 + IGFBP2 + SEPT7 + TFG")
  checkEquals(round(skatRes[[1]]$p.value, 4), 0.6578)
  
  checkEquals(skatRes[[2]]$null.model, "group ~ 1 + ageAtDx + ChemoAgent")
  checkTrue(!is.null(skatRes[[2]]$null.warnings))
  checkEquals(skatRes[[2]]$alternative.model, "group ~ 1 + ageAtDx + ChemoAgent + AMFR + ATP6V0A2 + ATRX + BAHD1 + BPHL + BRCA1 + BTD + CAMK2G + CCNF + CETN3 + CHD3 + CLP1 + CLPX + CPSF4 + CSTF3 + CYP4F12 + EIF5B + EP400 + ERAL1 + ERCC2 + EXTL3 + FANCG + FANCI + FDXR + FNTB + FOXD1 + FRYL + GALNT2 + GRIK5 + GRIP2 + GSK3B + HNRNPL + HTR7 + ILVBL + IMPA1 + INPP5E + JRK + KIAA0195 + KIAA0586 + KLHL18 + KRT33A + MC2R + MEA1 + MFN1 + MSH3 + MTX1 + MUTYH + NFRKB + NFYB + NKRF + NMT1 + NR2C1 + PAXIP1 + PCGF1 + PDXDC1 + PFDN6 + PHF10 + PIGB + PIGF + PIK3R2 + PLEKHB1 + POP4 + PPP5C + RAP1A + RBBP8 + REV3L + RFC5 + RPS6KB2 + SCAMP1 + SH2B1 + SLC24A1 + SLC25A11 + SLC30A3 + SPAST + SPRED2 + SSR1 + SYNJ2 + TAF2 + TMEM11 + TRIM27 + UBE4B + WDR62 + ZNF500 + ZNF592")
  checkEquals(round(skatRes[[2]]$p.value, 4), 0.2776)
  
  checkEquals(skatRes[[10]]$null.model, "group ~ 1 + ageAtDx + ChemoAgent")
  checkTrue(!is.null(skatRes[[10]]$null.warnings))
  checkEquals(skatRes[[10]]$alternative.model, "group ~ 1 + ageAtDx + ChemoAgent + ABCC9 + ABCE1 + ABHD14A + ADAMTS1 + ADRB3 + ALCAM + ANXA4 + APEX1 + AQP9 + ARMCX1 + ATIC + ATP10A + BAZ2B + BMP15 + BMX + BRIP1 + CA6 + CASP1 + CCDC102B + CD55 + CD74 + CDH8 + CDKAL1 + CFH + CIDEB + CKAP4 + CLCA4 + CLIP3 + CREB3L2 + CRTAP + CTNNA3 + CUL1 + DHODH + DSC3 + DTX4 + DUSP22 + EEF1G + EIF2B3 + EIF4EBP1 + ELF2 + ELN + EMD + ERAL1 + ERLIN1 + ETF1 + ETFB + EXOSC5 + FAIM + FAM49A + FXYD2 + FZD6 + GAS6 + GNG11 + GNL2 + GPD1L + GPRC5A + GPX1 + GPX7 + GRB14 + GSTM3 + GUCY1A3 + GUCY1B3 + HDGFRP3 + HES1 + HLA-DMA + IFI44 + IFI44L + IGLL1 + IL19 + IMPACT + IMPDH1 + INTS7 + KCNK10 + LAMC1 + LARS2 + LDB1 + LILRA4 + LUZP2 + MALL + MAP9 + MCTP2 + MDFIC + MEP1B + MEST + METTL8 + MFAP4 + MINA + MRE11A + MRPL11 + MRTO4 + MS4A5 + MSRB2 + MTCP1 + MYOZ2 + NBN + NEK11 + NOS1 + NOTCH1 + NUBPL + OBFC1 + PAICS + PBX1 + PCCB + PCNA + PDE3B + PDZRN3 + PFKM + PIM1 + PKD2 + PKIG + PLA2G4C + PLAC8 + POFUT1 + POLR3B + PON2 + PPAP2B + PPEF2 + PRDX4 + PRKCA + PRMT5 + PTGER3 + PTPLA + RAB40B + RNF130 + RNF144A + RPS8 + RRAS2 + RRP1 + RUVBL1 + S100G + SAMM50 + SAMSN1 + SCN3A + SCN9A + SEC31B + SELL + SH3BGR + SLC15A3 + SND1 + SOCS2 + STK3 + STOM + STOML2 + TBCE + TCF7L2 + TFAM + TFRC + TIMP1 + TIMP4 + TLR7 + TMED3 + TNFSF4 + TRIM29 + TUBA1A + TUFT1 + TYRP1 + UCHL5 + UCK2 + ULK4 + VIM + WDR3 + WDR41 + XKR8 + XPOT + XYLT1 + ZCCHC14 + ZMYND10 + ZNF132 + ZNF282 + ZNF550 + ZNF643")
  checkEquals(round(skatRes[[10]]$p.value, 4), 0.3497)
} # test_geneSetScoreTest