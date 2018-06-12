##' This function creat different Subset of the result table based on the cutoffs.
##' 
##' 
##' @title Subset of the result table 
##' @param x An object of class "netSEM", the returned list from \code{netSEMm}. Plotting uses the first element of this list (res.print) in which the first column of it is endogenous variable, second column is variable and other columns are corresponding best functional form, r-squared, adj-r-squared, P-value1, P-value2 and P-value3.
##' @param cutoff A threshold value for adjusted R-squared. The maximum number of cutoff is 3.
##' @return The dataframe of different subset of the result table.
##' @importFrom DiagrammeR grViz
##' @export
##'
##' @examples
##' # Load acrylic data set
##' data(acrylic)
##' # Build a semi-gSEM model
##' ans <- netSEMm(acrylic)
##' # Subset dataset with different cutoff
##' res <- subsetData(ans,cutoff=c(0.2))
##' res <- subsetData(ans,cutoff=c(0.2,0.5))
##' res <- subsetData(ans,cutoff=c(0.2,0.5,0.8))
##'

subsetData <- function(x,cutoff=c(0.2,0.5,0.8)){
  
  rtp1 <- x$table[,1:5]
  rtp1[,-c(1:3)] <- round(rtp1[,-c(1:3)],3)
  
  rtp1.a <- rtp1[-c(1:nrow(rtp1)),]
  rtp1.b <- rtp1[-c(1:nrow(rtp1)),]
  rtp1.c <- rtp1[-c(1:nrow(rtp1)),]
  rtp1.d <- rtp1[-c(1:nrow(rtp1)),]
  
  ### set up cutoffs ############
  count <- length(cutoff)
  
  if (count == 1) {
    c1 <- cutoff[1]  
    rtp1.a <- rtp1[rtp1[, "adj-R-Sqr"] <  c1 , ]
    rtp1.b <- rtp1[rtp1[, "adj-R-Sqr"] >= c1 , ]
  }
  if (count == 2) {
    c1 <- cutoff[1] ; c2 <- cutoff[2] 
    rtp1.a <- rtp1[rtp1[, "adj-R-Sqr"] < c1, ]
    rtp1.b <- rtp1[rtp1[, "adj-R-Sqr"] >= c1 & rtp1[, "adj-R-Sqr"] < c2, ]
    rtp1.c <- rtp1[rtp1[, "adj-R-Sqr"] >= c2 , ]
  }
  if (count == 3) {
    c1 <- cutoff[1] ; c2 <- cutoff[2] ; c3 <- cutoff[3] 
    rtp1.a <- rtp1[rtp1[, "adj-R-Sqr"] < c1, ]
    rtp1.b <- rtp1[rtp1[, "adj-R-Sqr"] >= c1 & rtp1[, "adj-R-Sqr"] < c2, ]
    rtp1.c <- rtp1[rtp1[, "adj-R-Sqr"] >= c2 & rtp1[, "adj-R-Sqr"] < c3, ]
    rtp1.d <- rtp1[rtp1[, "adj-R-Sqr"] >= c3 , ]
  }
  
  return(list(rtp1.a,rtp1.b,rtp1.c,rtp1.d))
  
}