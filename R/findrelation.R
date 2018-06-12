#' Find relationship between any two variables.
#'
#'
#' @title Find relationship between any two variables
#' @import segmented 
#' @importFrom stats predict
#' @importFrom stats lm
#' @importFrom stats residuals
#' @importFrom stats runif
#' @importFrom stats var
#' 
#' @param iVar The first parameter.
#' @param iResp The second parameter.
#' @param criterion What model fit parameter is used for ranking? (Default: adj.R2)
#' @param str TRUE/FALSE Is this a 'strength' type model fitting? (whether to use ^2 and sqrt functions)
#' @param modelNames A vector of model names chosen from "SL", "Quad", "SQuad", "Exp", "Log", "nls","CP"
#' @param cRes.all A dataframe of all results, passed in by the netSEMm function
#' @param cRes.best A dataframe of best results, passed in by the netSEMm function
#' @param cRes.print A dataframe of results to print, passed in by the netSEMm function
#' @param x A dataframe of data values, passed in by the netSEMm function
#' @param nlsInits A vector of nls initialization coefficients, passed in by the netSEMm function
#' @param nRes number of cells in the print variable, value passed in by the netSEmm function
#' @return a list of the following items: 
#' \itemize{
#'   \item "cRes.all": A dataframe of all results.
#'   \item "cRes.best": A dataframe of best results.
#'   \item "cRes.print": A dataframe of results to print
#' }
#' 
#'
#' 
#' @export
#' 



###############################
### The following function fits paired relations
###############################
findrelation <- function(iVar, iResp, criterion = "adj.R2", modelNames,str = FALSE,cRes.all,
                         cRes.best,cRes.print,x,nlsInits,nRes){
  #cat("\nRunning findrelation...\n")
  #cat("\niVar is:",iVar)
  #cat("\niResp is:",iResp,"\n")
  ## saves all adjusted R2 for all 8 functional forms
  adj.R2 <- rep(0, length(modelNames))
  ## Save "lm" object for all 8 functional forms
  model.res <- vector("list", length(modelNames))
  Resp.v <- x[[iResp]]          # data of main endogenous variable
  Var.v  <- x[[iVar]]           # data of predictor
  Resp   <- colnames(x)[iResp]  # name of endogenous variable
  Var    <- colnames(x)[iVar]   # name of predictor
  #cat("Regress Dep ~ Indep ", Resp, "~", Var, "\n")
  lm4.flag <- TRUE
  lm5.flag <- TRUE
  lm6.flag <- TRUE
  lm7.flag <- TRUE
  lm9.flag <- TRUE

  ##--------------------------------------------------
  ## Functional Form 1: Linear
  ##--------------------------------------------------
  Rel1 <- paste0(Resp, "~", Var)
  lm1 <- do.call("lm", list(Rel1, data = as.name("x")))
  adj.R2[1] <- summary(lm1)$adj.r.squared
  model.res[[1]] <- lm1
  cRes.all[[iVar, iResp, "SL"]] <- lm1

  ##--------------------------------------------------
  ## Functional Form 2: Quadratic
  ##--------------------------------------------------
  Rel2 <- paste0(Resp,"~",Var,"+","I(", Var,"^2)")
  lm2 <- do.call("lm", list(Rel2, data = as.name("x")))
  adj.R2[2] <- summary(lm2)$adj.r.squared
  model.res[[2]] <- lm2
  cRes.all[[iVar, iResp, "Quad"]] <- lm2

  ##--------------------------------------------------
  ## Functional Form 3: Quadratic (no linear term)
  ##--------------------------------------------------
  Rel3 <- paste0(Resp, "~", "I(", Var,"^2)")
  lm3 <- do.call("lm", list(Rel3, data = as.name("x")))
  adj.R2[3] <- summary(lm3)$adj.r.squared
  model.res[[3]] <- lm3
  cRes.all[[iVar, iResp, "SQuad"]] <- lm3

  ##--------------------------------------------------
  ## Functional Form 4: Exponential
  ##--------------------------------------------------
  lm4.flag <- all(exp(Var.v[!is.na(Var.v)]) != Inf)
  if (lm4.flag) {
    Rel4 <- paste0(Resp, "~", "exp(", Var,")")
    lm4 <- do.call("lm", list(Rel4, data = as.name("x")))
    adj.R2[4] <- summary(lm4)$adj.r.squared
    model.res[[4]] <- lm4
    cRes.all[[iVar, iResp, "Exp"]] <- lm4
  }else{
    adj.R2[4] <- -Inf
    model.res[[4]] <- NA
    cRes.all[[iVar, iResp, "Exp"]] <- NA
  }

  ##---------------------------------------------------
  ## Functional Form 5: Log
  ##---------------------------------------------------
  lm5.flag <- all(Var.v > 0 & Var.v < Inf, na.rm = TRUE)
  if(lm5.flag){
    Rel5 <- paste0(Resp, "~", "log(", Var,")")
    lm5 <- do.call("lm", list(Rel5, data = as.name("x")))
    adj.R2[5] <- summary(lm5)$adj.r.squared
    model.res[[5]] <- lm5
    cRes.all[[iVar, iResp, "Log"]] <- lm5
  }else{
    adj.R2[5] <- -Inf
    model.res[[5]] <- NA
    cRes.all[[iVar, iResp, "Log"]] <- NA
  }

  # This is to be able to turn off these two functional forms
  if(str){ # Begining of if(str)

    ##----------------------------------------------------
    ## Functional Form 6: SquareRoot (no linear term)
    ##----------------------------------------------------
    lm6.flag <- all(Var.v > 0 & Var.v < Inf, na.rm = TRUE)
    if(lm6.flag){
      Rel6 <- paste0(Resp, "~", "I(", Var,"^0.5)")
      lm6 <- do.call("lm", list(Rel6, data=as.name("x")))
      adj.R2[6] <- summary(lm6)$adj.r.squared
      model.res[[6]] <- lm6
      cRes.all[[iVar, iResp, "SQRoot"]] <- lm6
    }
    else {
      adj.R2[6] <- -Inf
      model.res[[6]] <- NA
      cRes.all[[iVar, iResp, "SQRoot"]] <- NA
    }

    ##-------------------------------------------------------
    ## Functional Form 7: Inverse SquareRoot (no linear term)
    ##-------------------------------------------------------
    lm7.flag <- all(Var.v > 0 & Var.v < Inf, na.rm = TRUE)
    if(lm7.flag){
      Rel7 <- paste0(Resp, "~", "I(1/(", Var,"^0.5))")
      lm7 <- do.call("lm", list(Rel7, data=as.name("x")))
      adj.R2[7] <- summary(lm7)$adj.r.squared
      model.res[[7]] <- lm7
      cRes.all[[iVar, iResp, "ISQRoot"]] <- lm7
    }
    else {
      adj.R2[7] <- -Inf
      model.res[[7]] <- NA
      cRes.all[[iVar, iResp, "ISQRoot"]] <- NA
    }

  } # End of if(str)

  ##--------------------------------------------------------
  ## Functional Form 8: nls
  ##--------------------------------------------------------
  ## if(Resp == "G" & Var == "time") browser()

  Rel8 <- paste0(Resp, " ~ ", "a1 + a2 * exp(a3 * ",Var, ")")
  model8 <- NA
  SSE <- NA
  ## Find the best (if any) nls models from all initial values
  for (i in 1:nrow(nlsInits)) {
    nls1 <- tryCatch({
      do.call("nls", list(Rel8, data = as.name("x"),
			         start = list(a1 = nlsInits[i, 1],a2 = nlsInits[i,2], a3 = nlsInits[i,3])))
    }, error = function(e) e)

    if (inherits(nls1, "nls")) {
      ## if model8 is still NA, initialize it
      if (!inherits(model8, "nls")) {
	       model8 <- nls1
	       SSE <- sum(residuals(nls1) ^ 2)
	       init <- nlsInits[i,]
      }else{
	       newSSE <- sum(residuals(nls1) ^ 2)
	       if (newSSE < SSE) {
	         model8 <- nls1
	         SSE <- newSSE
	         init <- nlsInits[i,]
	       }else{}
      }
    }
  }

  # This is needed to be able to drop out the 'str' funcforms
  #  and still have nls save into the cRes.all object correctly.
  # It will be '8' if the 'str' offs occupy '6' and '7',
  #  but otherwise it will be '6'.
  if (str) {
    nlsNum <- 8
  }else{
    nlsNum <- 6
  }

  ## Assign result
  if (!inherits(model8, "nls")) {
    adj.R2[nlsNum] <- -Inf
    model.res[[nlsNum]] <- NA
    cRes.all[[iVar, iResp, "nls"]] <- NA
  }else{
    R2s <- nlsR2(model8, y = x[[Resp]], p = 3)  ## apply 'nlsR2' function
    adj.R2[nlsNum] <- R2s$adjR2
    model.res[[nlsNum]] <- model8
    cRes.all[[iVar, iResp, "nls"]] <- model8
  }

  ##--------------------------------------------------
  ## Functional Form 9: Change Point
  ##--------------------------------------------------
  cpNum <- nlsNum + 1
  lm9 <- segmented.lm(lm(Resp.v~Var.v,data = x),seg.Z = ~Var.v, psi = NA, 
                      control = seg.control(K = 1,stop.if.error = FALSE,
                                            n.boot = 0,it.max = 20))
  lm9.flag <- all(summary(lm9)$adj.r.squared == summary(lm1)$adj.r.squared)
  if (lm9.flag) {
    adj.R2[cpNum] <- -Inf
    model.res[[cpNum]] <- NA
    cRes.all[[iVar, iResp, "CP"]] <- NA
  }else{
    adj.R2[cpNum] <- summary(lm9)$adj.r.squared
    model.res[[cpNum]] <- lm9
    cRes.all[[iVar, iResp, "CP"]] <- lm9
  }
  
  
  ## Choose the best model, if no model works, write -1.
  ## Populate the 'table' item in function return.
  ## Populate the 'bestModels' item in function return.
  if (max(adj.R2, na.rm = TRUE) >= 0.001) {
    
    step <- attributes(cRes.best)$Step[iResp] + 1
    attributes(cRes.best)$Step[iVar] <- min(step, attributes(cRes.best)$Step[iVar])
    attributes(cRes.best)$diag.Step[iVar, iResp] <- step
    nbest <- which.max(adj.R2) # location of best model
    cRes.best[iVar, iResp] <- modelNames[nbest]

    ## A row for the print table
    table1r <- matrix(NA, nrow = 1, ncol = nRes)
    colnames(table1r) <- c("endogenous", "Variable", "Model", "R-Sqr", "adj-R-Sqr",
			                     "Pval1", "Pval2", "Pval3")
    table1r <- as.data.frame(table1r)
    table1r[1, c("endogenous", "Variable")] = c(Resp, Var)
    table1r[1, c("Model")] <- modelNames[nbest]

    ## if(modelNames[nbest] == "nls") browser()
    ## Best model from Resp ~ Var
    finalModel <- model.res[[nbest]]
    if (inherits(finalModel, "lm")) {
      model <- summary(finalModel)
      table1r[1,c("R-Sqr", "adj-R-Sqr", "Pval1", "Pval2")] <-
	           c(model$r.sq, model$adj.r.sq, model$coef[1,4], model$coef[2,4])
      if (modelNames[nbest] == "Quad" | modelNames[nbest] == "CP" ) {
	       table1r[1, "Pval3"] <- model$coeff[3, 4]
      }else{
	       table1r[1, "Pval3"] <- NA
      }
    }else if (inherits(finalModel, "nls")) {
      table1r[1,c("R-Sqr", "adj-R-Sqr", "Pval1", "Pval2", "Pval3")] <-
	           c(R2s$R2, R2s$adjR2, NA, NA, NA)
    }
    cRes.print <- rbind(cRes.print, table1r)
  }else{
    cRes.best[iVar, iResp] <- -1
  }
  
  
  return(list(cRes.all,cRes.best,cRes.print))

  
}


