##' Extract and display an equation of a pairwise path between two variables.
##'
##' Extract the "best" model between any two variables. The model name and the model equation are printed on screen. The model coefficients, as well as the model R object are also returned.
##' @title Extract Path Model Equation
##' @param x object of class "netSEM", which is the return value of function \code{netSEMm}.
##' @param from character string. Name of the predictor.
##' @param to character string. Name of the endogenous variable.
##' @param round a positive integer. The coefficients are rounded to this decimal place.
##' @return a list of the following items: 
##' \itemize{
##'   \item "model": the best fitted model.
##'   \item "model.print": a character string of the model equation.
##' }
##' 
##' @export
##' @seealso \link[netSEM]{netSEMm}
##' @examples
##' ##' ## Load the sample acrylic data set
##' data(acrylic)
##'
##' ## Run netSEM principle three
##' ans <- netSEMm(acrylic)
##'
##' ## Extract relations between IrradTot and IAD2
##' cf <- path(ans,from ="IAD2",to="IrradTot")
##' print(cf)
##' 

path <- function(x, from, to, round = 3){
  ###------------------
  ### Model Names
  ###------------------
  ## IMPORTANT: These model names are taken from the main 'netSEM' function.
  ## It needs to be updated if they are changed in 'netSEM'
  modelNames <- c("SL", "Quad", "SQuad", "Exp", "Log", "SQRoot", "ISQRoot", "nls", "CPSLSL", "CP")

  ###------------------
  ### Checking arguments
  ###------------------
  if(!inherits(x, "netSEM"))
    stop("x is not of class 'netSEM'.")
  vars <- colnames(x$bestModels)
  if(from %in% modelNames)
    stop("'from' variable is not found.")
  if(to %in% modelNames)
    stop("'to' variable is not found.")

  ###------------------
  ### Extract coefficients
  ###------------------
  if(is.na(x$bestModels[from, to])){
    cat("No relation was found from", from, "to", to, "\n")
    return(NA)
  }else{
    relName <- x$bestModels[from, to]
    theModel <- x$allModels[[from, to, relName]]
    cat("Model type:", relName, "\n")
    if(relName == "SL"){
      coefs <- round(theModel$coef, round)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from)
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else if(relName == "Quad"){
      coefs <- round(theModel$coef, round)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from, cf(coefs[3]), " * ", from, "^2")
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else if(relName == "SQuad"){
      coefs <- round(theModel$coef, round)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from, "^2")
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else if(relName == "Exp"){
      coefs <- round(theModel$coef, round)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * exp(", from, ")")
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else if(relName == "Log"){
      coefs <- round(theModel$coef, round)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * log(", from, ")")
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else if(relName == "SQRoot"){
      coefs <- round(theModel$coef, round)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from, "^0.5")
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else if(relName == "ISQRoot"){
      coefs <- round(theModel$coef, round)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from, "^(-0.5)")
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else if(relName == "nls"){
       coefs <- round(summary(theModel)$coef[,1], round)
       ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * exp(", cf(coefs[3], TRUE), " * ", from, ")")
       cat("Model equation (round = ", round, "):\n", sep = "")
       print(ans.print)
    } else if(relName == "CPSLSL"){
      ans.print <- NA
    } else if(relName == "CP"){
      coefs <- round(theModel$coef, round)
      breakpoint <- round(theModel$psi[2], 1)
      ans.print <- paste0(to, " = ", cf(coefs[1], TRUE), cf(coefs[2]), " * ", from, cf(coefs[3]), " * (", from,"-",breakpoint, ")_+")
      cat("Model equation (round = ", round, "):\n", sep = "")
      print(ans.print)
    } else {
      ans.print <- NA
    }
    invisible(list(model = theModel, model.print = ans.print))
  }
}
