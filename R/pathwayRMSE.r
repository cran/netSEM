##' Summarize root mean square error (RMSE) for direct and indirect pathway from netSEMp1 result.
##'
##' pathwayRMSE gives a summary about RMSE.
##'
##' @title Summary of pathway RMSE
##' @param x An object of class "netSEM", the returned list from netSEMm().
##' @param maxlen The maximum length of chosen mechanism.
##' @param ... A S3 generic/method consistency.
##' @return A data frame of result. A summary of RMSE results for all of pathways.
##' @importFrom gtools permutations 
##' @importFrom stats predict
##' @export
##'
##' @seealso \link[netSEM]{netSEMp1}, \link[netSEM]{pathwayPredict}
##'
##' @examples
##' \dontrun{
##' data(acrylic)
##' ans <- netSEMp1(acrylic)
##' pathwayRMSE(ans,maxlen=2)
##' }


pathwayRMSE <- function(x, maxlen=2, ... ) { 
  
  ## Define stressor, mechanism, and response variable
  S <- x$exogenous[1]
  M <- x$exogenous[2:length(x$exogenous)]
  R <- x$endogenous # Response
  
  if (maxlen > length(M)) {
    stop("The value of maxlen must be less than the number of mechanism")
  }
  
  result <- data.frame()
  
  ## Direct pathway
  relName <- x$bestModels[S,R] # Type of best model between S and R
  theModel <- x$allModels[[S,R,relName]]
  pred <- predict(theModel)
  RMSE_dir <- round(sqrt(mean((pred - x$data[, x$endogenous])^2)), 4)
  result[1,1] <- 0 # 1st column is length
  result[1,2] <- paste0(S, "-->", R) # 2nd column is path
  result[1,3] <- RMSE_dir # 3rd column is RMSE
  colnames(result) <- c("length", "path", "RMSE")
  
  for (k in 1:maxlen) {
    # Pick k mechanism without replacement
    # Get all permutations
    nM <- length(M)
    Pr <- permutations(n = nM, r = k, v = M) # Pr is a matrix where each row contains a vector of length r
    temp <- result # May contain NA's
    res <- data.frame()
    for (j in 1:nrow(Pr)) {
      path <- paste0(S, "-->", paste(Pr[j, ], collapse = "-->"), "-->", R)
      suppressWarnings(RMSE <- pathwayPredict(x, path)$RMSE)
      res[j, 1] <- k
      res[j, 2] <- paste0(S, "-->", paste(Pr[j, ], collapse = "-->"), "-->", R)
      res[j, 3] <- RMSE
      colnames(res) <- c("length", "path", "RMSE")
    }
    result <- rbind(temp, res)
  }
  result_final <- result[which(!is.na(result$RMSE)), ] # Remove rows with NA
  rownames(result_final) <- 1:nrow(result_final)
  return(data.frame(result_final))
}