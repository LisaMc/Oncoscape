#############################################
#Setup analysis data for a given gene set   #
#Author:  Wade K. Copeland                  #
#############################################
setClass("analysisData") #Virtual class, class that exists for which no actual objects can be created by a call to new

setMethod("analysisDataSetup", signature(object = "GeneSetBinomialMethods"),
          function (object, sampleIDsG1, sampleIDsG2, covariates = NULL, geneSet,sampleDescription="", geneSetDescription="") {
            myGenes <- getGeneSets(object)[geneSet][[1]]
            
            if(length(which(colnames(getExpressionData(object)) %in% myGenes)) == 0) {
              stop("None of the genes in the specified gene sets are contained in the gene expression data.")
            }
            
            locG1 <- which(rownames(getClinicalData(object)) %in% sampleIDsG1)
            subjG1 <- rownames(getClinicalData(object)[locG1, ])
            clinicalG1 <- data.frame(getClinicalData(object)[locG1, covariates], row.names = subjG1, stringsAsFactors = FALSE)
            colnames(clinicalG1) <- covariates
            clinicalG1$group <- 1
            
            locG2 <- which(rownames(getClinicalData(object)) %in% sampleIDsG2)
            subjG2 <- rownames(getClinicalData(object)[locG2, ])
            clinicalG2 <- data.frame(getClinicalData(object)[locG2, covariates], row.names = subjG2, stringsAsFactors = FALSE)
            colnames(clinicalG2) <- covariates
            clinicalG2$group <- 0
            
            clinical <- rbind(clinicalG1, clinicalG2)
            #which(rownames(clinical) %in% sampleIDsG1)
            #which(rownames(clinical) %in% sampleIDsG2)
            
            if(length(which(colnames(getExpressionData(object)) %in% myGenes)) > 1) {
              geneExpression <- getExpressionData(object)[, which(colnames(getExpressionData(object)) %in% myGenes)]
              #dim(geneExpression)
              
              unmatchedGenes = setdiff(myGenes, colnames(geneExpression))
              
              matchedGenes <- colnames(geneExpression)
            } else {
              geneExpression <- getExpressionData(object)[, which(colnames(getExpressionData(object)) %in% myGenes)]
              #length(geneExpression)
              
              geneExpression <- matrix(geneExpression, 
                                       dimnames = list(rownames(getExpressionData(object)), 
                                                       colnames(getExpressionData(object))[which(colnames(getExpressionData(object)) %in% myGenes)]))
              
              unmatchedGenes = setdiff(myGenes, colnames(geneExpression))
              
              matchedGenes <- colnames(geneExpression)
            }
            
            analysisData <- merge(clinical, geneExpression, by = "row.names")
            rownames(analysisData) <- analysisData$Row.names
            analysisData <- analysisData[, -1]
            
            #Of the samples provided, which are contained in both the clinical and gene expression data?
            unmatchedSamples <- setdiff(unlist(c(sampleIDsG1, sampleIDsG2)), rownames(analysisData))
            
            matchedSamples <- rownames(analysisData)
            
            output <- list(sampleDescription = sampleDescription,
                           geneSetDescription = geneSetDescription,
                           matchedSamples = matchedSamples,
                           unmatchedSamples = unmatchedSamples,
                           matchedGenes = matchedGenes,
                           unmatchedGenes = unmatchedGenes,
                           analysisData = analysisData)
            
            
            output <- structure(output, class = "analysisData")
            
            return(output)
          })