#####################################
#Function to take random samples    #
#Author:  Wade K. Copeland          #
#####################################
setClass("randomSample") #Virtual class, class that exists for which no actual objects can be created by a call to new

setMethod("randomSample", signature(object = "GeneSetBinomialMethods"),
          function(object, survivalVar, cut, nG1 = NULL, nG2 = NULL, all = FALSE, seed = sample(.Random.seed, 1)) {
            
            if((is.null(nG1) | is.null(nG2)) & all == FALSE) {
              stop("If either of the arguments nG1 or nG2 is null, then all must equal TRUE")
            }
            
            if((length(which(getClinicalData(object)[, survivalVar] <= cut)) == 0 | length(which(getClinicalData(object)[, survivalVar] > cut)) == 0) == TRUE) {
              stop(paste("Either one or both groups contain no subjects for the given cut value.  The survival range is [", 
                         range(na.omit(getClinicalData(object)[, survivalVar]))[1], ", ", range(na.omit(getClinicalData(object)[, survivalVar]))[2], "].", sep = ""))
            }
            
            if(all == FALSE) {
              if(length(which(getClinicalData(object)[, survivalVar] <= cut)) < nG1) {
                nG1 <- length(which(getClinicalData(object)[, survivalVar] <= cut))
                warning(paste("nG1 is truncated to have length", nG1)) 
              }
              
              if(length(which(getClinicalData(object)[, survivalVar] > cut)) < nG2) {
                nG2 <- length(which(getClinicalData(object)[, survivalVar] > cut))
                warning(paste("nG2 is truncated to have length", nG2)) 
              }
              
              set.seed(seed)
              shortSurvivors <- as.list(row.names(getClinicalData(object)[sample(which(getClinicalData(object)[, survivalVar] <= cut), nG1), ]))
              
              set.seed(seed)
              longSurvivors <- as.list(row.names(getClinicalData(object)[sample(which(getClinicalData(object)[, survivalVar] > cut), nG2), ]))
              
              res <- list(shortSurvivors = shortSurvivors, longSurvivors = longSurvivors, seed = seed)
              res <- structure(res, class = "randomSample")
              return(res)
            } else {
              shortSurvivors.all <- as.list(rownames(getClinicalData(object)[which(getClinicalData(object)[, survivalVar] <= cut), ]))
              longSurvivors.all <- as.list(rownames(getClinicalData(object)[which(getClinicalData(object)[, survivalVar] > cut), ]))
              
              res <- list(shortSurvivors = shortSurvivors.all, longSurvivors = longSurvivors.all)
              res <- structure(res, class = "randomSample")
              return(res)
            }
          })