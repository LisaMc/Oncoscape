############################################################################
#Class generator function for the class to contain the data                #
#Author:  Wade K. Copeland                                                 #
############################################################################
GeneSetBinomialMethods <- setClass("GeneSetBinomialMethods", 
                                   slots = c(expressionData = "matrix",
                                             clinicalData = "data.frame",
                                             geneSets = "list"),
                                   validity  = function(object) {
                                     rowExpr <- rownames(object@expressionData)
                                     rowClinical <- rownames(object@clinicalData)
                                     
                                     if(length(intersect(rowExpr, rowClinical)) == 0) { 
                                       return("There are no shared samples between the expression data and clinical data.")
                                     }
                                     
                                     if(any(duplicated(rowExpr)) == TRUE) {
                                       return("None of the methods for objects of class GeneSetBinomialMethods can handle duplicated observations.  Please input expression data that does not have duplicated observations.")
                                     }
                                     
                                     if(!suppressWarnings(any(is.na(as.numeric(rowExpr)))) & !suppressWarnings(any(is.na(as.numeric(rowClinical))))) {
                                       warning("The rownames of both the expression data and clinical data are numeric.  Please verify these are the unique ID values you want to use.")
                                     }
                                     
                                     if(!any(unique(unlist(object@geneSets)) %in% unique(colnames(object@expressionData)))) {
                                       warning("None of the genes in the specified geneset are in the supplied expression data.")
                                     }
                                     
                                     return(TRUE)
                                   }
)

#################################################
#Show method for the class generator function   #
#Author:  Wade K. Copeland                      #
#################################################
setMethod("show", signature(object = "GeneSetBinomialMethods"), 
          function(object) {          
            objectInfo <- paste("This is an object of class GeneSetBinomialMethods.\n",
                                "The expression data has", dim(object@expressionData)[1], "observations and", dim(object@expressionData)[2], "genes (getExpressionData()).\n",
                                "The clinical data has", dim(object@clinicalData)[1], "observations and", dim(object@clinicalData)[2], "variables (getClinicalData()).\n",
                                "There are", length(object@geneSets), "gene sets (getGeneSets()).")
            
            cat(objectInfo)
            
            invisible(objectInfo)
          }         
)

################################
#Non-s4 helper functions       #
#Author:  Wade K. Copeland     #
################################
#Function to retrieve warnings without sending them to the standard output.  This is an internal!
.withWarnings <- function(expr) {
  myWarnings <- NULL
  
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  
  val <- withCallingHandlers(expr, warning = wHandler)
  
  return(list(value = val, warnings = myWarnings))
}

################################
#Create some generic functions #
#Author:  Wade K. Copeland     #
################################
setGeneric('getExpressionData', function(object, ...) standardGeneric ('getExpressionData'))
setGeneric('getClinicalData', function(object, ...) standardGeneric ('getClinicalData'))
setGeneric('getGeneSets', function(object, ...) standardGeneric ('getGeneSets'))
setGeneric("randomSample", function(object, ...) standardGeneric ('randomSample'))
setGeneric("analysisDataSetup", function(object, ...) standardGeneric ('analysisDataSetup'))

##############################################
#GeneSetBinomialMethods accessor functions   #
#Author:  Wade K. Copeland                   #
##############################################
setMethod('getExpressionData', signature(object = "GeneSetBinomialMethods"),
          function(object) { 
            return(object@expressionData)
          })

setMethod('getClinicalData', signature(object = "GeneSetBinomialMethods"),
          function(object) { 
            return(object@clinicalData)
          })

setMethod('getGeneSets', signature(object = "GeneSetBinomialMethods"),
          function(object) { 
            return(object@geneSets)
          })