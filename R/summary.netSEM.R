##' summarize netSEM result
##'
##' summary.netSEM gives a summary about the net-SEM analysis.
##'
##' @title Summary of net-SEM
##' @param object An object of class "netSEM", the returned list from netSEMm().
##' @param ... A S3 generic/method consistency.
##' @return NULL. A summary of data and fitting result is printed on screen.
##'
##' @export
##'
##' @examples
##' data(acrylic)
##' ans <- netSEMm(acrylic)
##' summary(ans)

summary.netSEM <- function(object, ...){
  cat("netSEM: It consideres all exagenous variables and finds the best regression model \n")
  cat("\n")
  cat("Exogenous (Multiple exogenous):", colnames(object$bestModels)[-1], "\n")
  cat("Endegenous:", colnames(object$bestModels)[1], "\n")

  cat("\n")
  cat("Chosen models:\n")
  
  sapply(1:nrow(object$table), function(i) {
    paste(object$table[i,2], " ---> ", object$table[i,1],
          " | ", object$table[i,3], " (Adjusted R-square = ",
          round(as.numeric(object$table[i,5]), 5), ")", collapse = "\n",sep = "")
    }
  )
}
