#' Predict response variable values along a pathway
#'
#' @title Summary of predicted and observed response values along a pathway
#' 
#' @param x An object of class "netSEM", the returned list from \code{netSEMm}.
#' @param path A string form for a pathway, the default output format from \link[netSEM]{pathwayRMSE}.
#' @param newdata A data frame of the stress variable. The default is NULL.
#'
#' @importFrom stats predict
#' @seealso \link[netSEM]{pathwayRMSE}, \link[netSEM]{netSEMm}
#' @return An object of class pathway, which is a list of the following items:
##'
##' \itemize{
##' \item "pathway": A string form that shows the pathway.
##' \item "RMSE": A value of the root mean squared error.
##' \item "Resp": A matrix. The first column is the observed response values and 
##' the second is the predicted response values.
##' }
#' @export
#' 
#' @examples
#' # Load the sample acrylic data set
#' data(acrylic)
#' ans <- netSEMm(acrylic)
#' paths <- pathwayRMSE(ans,maxlen=3)
#' response <- pathwayPredict(ans, paths[10,2])
#' response
#' 

pathwayPredict <- function(x, path, newdata = NULL) {
  
  if (!"netSEM" %in% class(x)) stop("object must be of class 'netSEM'")
  
  vars <- unlist(strsplit(path, "-->"))
  nVars <- length(vars)
  
  if (is.null(newdata)) {
    stressVals <- data.frame(unlist(x$data[vars[1]]))
  } else {
    stressVals <- data.frame(newdata)
  }
  
  colnames(stressVals) <- vars[1]
  resp <- data.frame(rep(0, nrow(stressVals)))
  for (i in 1:(nVars - 1)) {
    type <- try(x$bestModels[vars[i], vars[i+1]], silent = TRUE)
    if (class(type) == "try-error") {
      warning(sprintf("Path between %s and %s not found in netSEM object, returning NA", vars[i], vars[i+1]))
      return(list(pathway=path,RMSE=NA,Resp=NA)) 
    } else if (type %in% c(NA,-1)) {
      warning(sprintf("No statistically significant path between %s and %s, returning NA", vars[i], vars[i+1]))
      return(list(pathway=path,RMSE=NA,Resp=NA)) 
    }
    if (type == "CP") {
      colnames(resp) <- "Var.v"
      if (i == 1) {
        colnames(stressVals) <- "Var.v"
      }
    }
    model <- x$allModels[[vars[i], vars[i+1], type]]
    if (i != 1) {
      resp <- data.frame(predict(model, resp))
      colnames(resp) <- vars[i+1]
    } else {
      resp <- data.frame(predict(model, stressVals))
      colnames(resp) <- vars[i+1]
    }
  }
  pred.resp <- unname(unlist(resp))
  RMSE <- round(sqrt(mean((pred.resp-x$data[,x$endogenous])^2)),4)
  resp.final <- data.frame(obs.resp=x$data[,x$endogenous],pred.resp=round(pred.resp,3))
  return(list(pathway=path,RMSE=RMSE,Resp=resp.final))
  
}
