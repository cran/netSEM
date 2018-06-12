##' summarize root mean square error (RMSE) for direct and indrect pathway from netSEM result
##'
##' pathwayRMSE gives a summary about RMSE
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
##' @seealso \link[netSEM]{netSEMm}, \link[netSEM]{pathwayPredict}
##'
##' @examples
##' data(acrylic)
##' ans <- netSEMm(acrylic)
##' pathwayRMSE(ans,maxlen=2)
##' 


pathwayRMSE <- function(x,maxlen=2, ... ){
  
  ## define stressor, mechanism, and response variable
  S <- x$exogenous[1]
  M <- x$exogenous[2:length(x$exogenous)]
  R <- x$endogenous
  
  if(maxlen > length(M)){
    stop("The value of maxlen must be less than the number of mechanism")
  }
  
  result <- data.frame()
  
  ## direct pathway
  relName <- x$bestModels[S,R]
  theModel <- x$allModels[[S,R,relName]]
  pred <- predict(theModel)
  RMSE_dir <- round(sqrt(mean((pred-x$data[,x$endogenous])^2)),4)
  result[1,1] <- 0
  result[1,2] <- paste0(S,"-->",R)
  result[1,3] <- RMSE_dir
  colnames(result) <- c("length","path","RMSE")
  
  for(k in 1:maxlen){
    #pick k mechanism without replacement
    #get all permutations
    nM <- length(M)
    Pr <- permutations(n=nM,r=k,v=M)
    temp <- result
    res <- data.frame()
    for(j in 1:nrow(Pr)){
      path <- paste0(S,"-->",paste(Pr[j,],collapse="-->"),"-->",R)
      suppressWarnings(RMSE <- pathwayPredict(x, path)$RMSE)
      res[j,1] <- k
      res[j,2] <- paste0(S,"-->",paste(Pr[j,],collapse="-->"),"-->",R)
      res[j,3] <- RMSE
      colnames(res) <- c("length","path","RMSE")
    }
    result <- rbind(temp,res)
  }
  result_final <- result[which(!is.na(result$RMSE)),]
  rownames(result_final) <- 1:nrow(result_final)
  return(data.frame(result_final))
}


