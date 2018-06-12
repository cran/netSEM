# This function Calculate R^2 and adj R^2 from nls result

nlsR2 <- function(nlsAns, # return from a nls() call
                  y, # endogenouse vector
                  p) # number of parameters
{
  n <- length(y)
  rsds <- residuals(nlsAns)
  rr <- sum(rsds^2) / var(y) / (length(y) - 1)
  R2 <- 1 - rr
  adjR2 <- 1 - rr * (n - 1) / (n - p)
  return(list(R2 = R2, adjR2 = adjR2))
}
