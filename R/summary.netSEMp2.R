##' summarize netSEM result
##'
##' summary.netSEMp2 gives a breif summary about the netSEM Principle 2 analysis.
##'
##' @title Summary of netSEMp2
##' @param object An object of class "netSEMp2", the returned list from netSEMp2().
##' @param ... A S3 generic/method consistency.
##' @return NULL. A summary of data and fitting result is printed on screen.
##'
##' @export
##'
##' @examples
##' data(acrylic)
##' ans <- netSEMp2(acrylic)
##' summary(ans)
##' 
##' 

summary.netSEMp2 <- function(object, ...) {
  cat("netSEMp2: It consideres all exagenous variables and finds the multiple regression model \n")
  cat("\n")
  #cat("Main predictor:", colnames(object$res.best)[1], "\n")
  #cat("Response:", colnames(object$res.best)[2], "\n")
  #cat("Intermediate variables:", colnames(object$res.best)[-c(1,2)], "\n")
  #cat("\n")
  cat("Chosen models:\n")
  sapply(1:nrow(object$res.print), function(i) {
    paste(object$res.print[i,2], "--->", object$res.print[i,1], "|", "GModel:", object$res.print[i,3], 
          "(Adjusted R-square = ", round(as.numeric(object$res.print[i,4]), 3), ")", collapse="\n")
  }
  )
}