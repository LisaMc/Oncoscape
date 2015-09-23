#####################################
#Function to take random samples    #
#Author:  Wade K. Copeland          #
#####################################
setClass("geneSetScoreTest", contains = "analysisData")
setGeneric("geneSetScoreTest", function(data, ...) standardGeneric ('geneSetScoreTest'))

setMethod("geneSetScoreTest", signature(data = "analysisData"),
          function(data) {
            skatData <- data
            
            covariates <- setdiff(colnames(skatData$analysisData), c("group", skatData$matchedGenes))

            #Replace the factor level "" with a missing value.
            for(i in length(covariates)) {
              skatData$analysisData[which(skatData$analysisData[, covariates[i]] == ""), covariates[i]] <- NA
            }
            
            #Data set that only includes the genes (note that the genes matrix cannot have any missing values)
            myGenes <- data.matrix(skatData$analysisData[, (-1)*which(colnames(skatData$analysisData) %in% c(covariates, "group"))])
            
            colnames(myGenes) <- colnames(skatData$analysisData)[(-1)*which(colnames(skatData$analysisData) %in% c(covariates, "group"))]
            
            rownames(myGenes) <- rownames(skatData$analysisData)
            
            myGenes <- na.omit(myGenes)
            
            #Analysis data that does not include missing values
            skatData$analysisData <- skatData$analysisData[rownames(myGenes), ]
            
            
            if(length(covariates) == 0) {
              myformula <- as.formula(paste("skatData$analysisData[, 'group'] ~ 1"))
              
              null <- .withWarnings(SKAT_Null_Model(myformula, out_type = "D", Adjustment = FALSE))
              
              skatRes <- SKAT(myGenes, null$value, kernel = "linear", is_check_genotype = FALSE)
              
              mNull <- paste("group ~ 1", sep = "")
              
              mAlt <- paste("group ~ 1 + ", paste(colnames(myGenes), collapse = " + "), sep = "")
            } else {
              myformula <- as.formula(paste("skatData$analysisData[, 'group'] ~ 1 +",  paste("skatData$analysisData[ ,'", covariates, "']", sep = "", collapse = " + ")))
              
              null <- .withWarnings(SKAT_Null_Model(myformula, out_type = "D", Adjustment = FALSE))
              
              skatRes <- SKAT(myGenes, null$value, kernel = "linear", is_check_genotype = FALSE)
              
              mNull <- paste("group ~ 1 + ", paste(covariates, collapse = " + "), sep = "")
              
              mAlt <- paste("group ~ 1 + ", paste(covariates, collapse = " + "), " + ", paste(colnames(myGenes), collapse = " + "), sep = "")
            }
            
            res <- list(
              sampleDescription = skatData$sampleDescription, 
              geneSetDescription = skatData$geneSetDescription, 
              analysisData = skatData$analysisData, 
              matchedSamples = skatData$matchedSamples, 
              unmatchedSamples = skatData$unmatchedSamples,
              matchedGenes = skatData$matchedGenes,
              unmatchedGenes = skatData$unmatchedGenes, 
              null.model = mNull, 
              null.warnings = unlist(null$warnings),
              alternative.model = mAlt, 
              p.value = skatRes$p.value
            )
            
            res <- structure(res, class = "geneSetScoreTest")
            return(res)
          })

setMethod("summary", signature(object = "geneSetScoreTest"), 
          function(object) {   
            
            if(is.null(object$null.warnings)) {
              summary.object <- paste("Null Model: ", object$null.model, "\nAlternative Model: ", object$alternative.model, "\n\nP-value:", round(object$p.value, 4))
            } else {
              summary.object <- paste("Null Model: ", object$null.model, "\nAlternative Model: ", object$alternative.model, "\n\nP-value:", round(object$p.value, 4), 
                                      "\n\nNull Model Warnings:", paste(unlist(object$null.warnings), collapse = "; ")) 
            }
            
            cat(summary.object)
            
            invisible(summary.object)
          }         
)