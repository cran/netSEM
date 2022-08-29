##' This function builds an netSEM model using principle 2. 
##' 
##' Principle 2 resembles the multiple regression principle in the way multiple predictors are considered simultaneously. Specifically, the first-level predictors to the system level variable, such as, Time and unit level variables, acted on the system level variable collectively by an additive model. This collective additive model can be found with a generalized stepwise variable selection (using the stepAIC() function in R, which performs variable selection on the basis of AIC or BIC) and this proceeds iteratively.
##' 
##' Data is analysed first using Principle 1 to find the best models. If needed, transformations based on the best models are applied to the predictors. Starting from the system response variable, each variable is regressed on all other variables except for the system response in an additive multiple regression model, which is reduced by a stepwise selection using stepAIC(). Then, for each selected variable, fitted regression for 6 selected functional forms (8 if strength type problem) and pick the best.
##'
##' @title network Structural Equation Modelling (netSEM) - Principle 2


##' @param x A dataframe. By default it considers all columns as exogenous variables, except the first column which stores the system endogenous variable.
##' @param exogenous, by default it consideres all columns as exogenous variables except column number 1, which is the main endogenous response.
##' @param endogenous A character string of the column name of the main endogenous OR a numeric number indexing the column of the main endogenous.
##' @param str A boolean, whether or not this is a 'strength' type problem.
##' @param criterion By default uses AIC to identify best model. Alternative choice is BIC.
##' @return A list of the following items:
##'
##' \itemize{
##' \item "res.print": A matrix. For each row, first column is the endogenous variable, second column is the exogenous variable, 
##' the other columns show corresponding summary information.
##'}

##' @seealso \link[netSEM]{netSEMp1} 
##'
##' @export
##' 
##' @import MASS
##' @importFrom stats coef
##' @examples
##' # Using AIC criterion
##' data(acrylic)
##' ans <- netSEMp2(acrylic, criterion = "AIC")
##' 
##' # Using AIC criterion
##' ans <- netSEMp2(acrylic, criterion = "BIC")
##' 
##' 
##' \dontrun{
##' # Using simulated data
##'  s <- runif(100,0,2)
##' m3 <- 1+2.5*s+rnorm(100,0,0.5)
##' m1 <- runif(100,1,4)
##' m2 <- -1-m1+m3+rnorm(100,0,0.3)
##'  y <- 2+2*exp(m1/3)+(m2-1)^2-m3+rnorm(100,0,0.5)
##' # Check the pairwise plot 
##' sim <- data.frame(cbind(y,s,m1,m2,m3))
##' pairs(sim)
##' ans <- netSEMp2(sim)
##' 
##' }

 
netSEMp2 <- function(x, exogenous = NULL, endogenous = NULL, str = FALSE, criterion = "AIC"){ 
  # Default is when str = FALSE (not considering strength type problem)
  
  # Assign 1st column to exogenous if exogenous is not specified in function
  # Define range of exogenous variable
  if(!missing(exogenous)){
    if(!is.character(exogenous)){
      if(is.wholenumber(exogenous)){
        if(exogenous < 1 | exogenous > length(colnames(x)))
          stop("exogenous location out of range!")
        exogenous.loc <- exogenous
      }
    }else{ 
      if(!(exogenous %in% colnames(x))){
        stop(paste0("exogenous '", exogenous, "' does not exist!"))
      }
      exogenous.loc <- which(colnames(x) == exogenous)
    }
    neworder <- 1:length(colnames(x))
    neworder[1] <- exogenous.loc
    neworder[exogenous.loc] <- 1
    x <- x[neworder]
  }
  
  # Define range of endogenous variable
  if(!missing(endogenous)){
    endogenous.loc <- which(colnames(x) == endogenous)
    if(!is.character(endogenous)){
      if(is.wholenumber(endogenous)){
        if(endogenous < 1 | endogenous > length(colnames(x)))
          stop("endogenous location out of range!")
        endogenous.loc <- endogenous
      }
    }else{
      if(!(endogenous %in% colnames(x))){
        stop(paste0("endogenous '", endogenous, "' does not exist!"))
      }
      endogenous.loc <- which(colnames(x) == endogenous)
    }
    neworder <- 1:length(colnames(x))
    neworder[2] <- endogenous.loc
    neworder[endogenous.loc] <- 2
    x <- x[neworder]
  }
  
  ###############################
  ## The following function does multiple selection
  ############################### 
  Multiple.relation <- function(Mx){
    # status_matrix (matrix) appears later
    # i_status starts from 1
    status_ind <- status_matrix[i_status,]
    Resp <- colnames(status_matrix)[which(status_ind==1)]
    Var <- colnames(status_matrix)[which(status_ind==0)]
    ## Get rid of NA's
    # Mx is a data frame from Multiple.relation function
    x1 <- Mx[c(Resp,Var)]
    x1 <- x1[apply(x1, 1, FUN = function(x) sum(is.na(x)) == 0),]
    
    ################## START: APPLY TRANSFORMATION TO PREDICTORS ###################
    # Define the equation for each of the linear and nonlinear functional forms
    trans <- rep(NA,length=length(Var)) # Save equations with coefficients
    trans1 <- rep(NA,length=length(Var)) # Save correlations without coefficients
    names(trans) <- Var
    for ( i in 1:length(Var) ) {
      # Use best models from netSEMp1 (p1.bestModels)
      # If "SL", x ---> x
      if ( p1.bestModels[Var[i],Resp] == "SL" ) {
        trans1[i] <- NA
        names(trans1)[i] <- Var[i]
        Quad_value <- coef(p1.allModels[Var[i],Resp,'SL'][[1]]) # Quad_value[1] is intercept
        if ( Quad_value[2]>0 ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),"*",Var[i],sep="")
        } else {  # When Quad_value[2]<0
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),"*",Var[i],sep="")
        }
      }
      # If "SQuad", x ---> x^2
      # Simple Quadratic relation.
      if ( p1.bestModels[Var[i],Resp] == "SQuad" ) {
        x1[,Var[i]] <- (x1[,Var[i]])^2
        trans1[i] <- paste(Var[i],"_t=", Var[i],"^2",sep="")
        names(trans1)[i] <- Var[i]
        names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
        Quad_value <- coef(p1.allModels[Var[i],Resp,'SQuad'][[1]])
        if ( Quad_value[2]>0 ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),"*",Var[i],"^2",sep="")
        } else {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),"*",Var[i],"^2",sep="")
        }
      }
      # If "Exp", x ---> exp(x)
      if ( p1.bestModels[Var[i],Resp] == "Exp" ) {
        x1[,Var[i]] <- exp(x1[,Var[i]])
        trans1[i] <- paste(Var[i],"_t=exp(", Var[i],")",sep="")
        names(trans1)[i] <- Var[i]
        names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
        Quad_value <- coef(p1.allModels[Var[i],Resp,'Exp'][[1]])
        if ( Quad_value[2]>0 ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),"*e^",Var[i],sep="")
        } else {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),"*e^",Var[i],sep="")
        }
      }
      # If "SQRoot"", x ---> sqrt(x)
      # For strength type problem
      if ( p1.bestModels[Var[i],Resp] == "SQRoot") {
        x1[,Var[i]] <- sqrt(x1[,Var[i]])
        trans1[i] <- paste(Var[i],"_t=sqroot(",Var[i],")",sep="")
        names(trans1)[i] <- Var[i]
        names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="") 
        Quad_value <- coef(p1.allModels[Var[i],Resp,'SQRoot'][[1]])
        if ( Quad_value[2]>0 ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),"*sqroot_",Var[i],sep="") 
          
        } 
        
        else {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),"*sqroot_",Var[i],sep="")
          
        }
      }
      # If "ISQRoot", x ---> 1/sqrt(x)
      # For strength type problem
      if ( p1.bestModels[Var[i],Resp] == "ISQRoot") {
        x1[,Var[i]] <- (1/sqrt(x1[,Var[i]]))
        trans1[i] <- paste(Var[i],"_t=isqroot(",Var[i],")",sep="")
        names(trans1)[i] <- Var[i]
        names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
        Quad_value <- coef(p1.allModels[Var[i],Resp,'ISQRoot'][[1]])
        if ( Quad_value[2]>0 ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),"*isqroot_",Var[i],sep="") 
          
        } 
        else {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),"*isqroot_",Var[i],sep="")
        }
      }
      # If "Log", x ---> log(x)
      if ( p1.bestModels[Var[i],Resp] == "Log" ) {
        x1[,Var[i]] <- log(x1[,Var[i]])
        trans1[i] <- paste(Var[i],"_t=log(", Var[i],")",sep="")
        names(trans1)[i] <- Var[i]
        names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
        Quad_value <- coef(p1.allModels[Var[i],Resp,'Log'][[1]])
        if ( Quad_value[2]>0 ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),"*log_",Var[i],sep="")
        } else {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),"*log_",Var[i],sep="")
        }
      }
      # If "Quad", x ---> (x+b/2c)^2
      if ( p1.bestModels[Var[i],Resp] == "Quad" ) {
        temp_column <- x1[,Var[i]]^2
        x1 <- cbind(x1,temp_column) # Adds additional column called temp_column to x1
        names(x1)[ncol(x1)] <- paste(Var[i],"__2",sep="")
        
        Quad_value <- coef(p1.allModels[Var[i],Resp,'Quad'][[1]])
        names(trans1)[i] <- Var[i]
        # Quad_value[1] is intercept, Quad_value[2] is coefficient, Quad_value[3] is coefficient only for Quad (x^2)
        if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),'*',Var[i],'+',
                            format(Quad_value[3],scitific=T,digits=4),'*',Var[i],"^2",sep="")
        }
        if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),'*',Var[i],
                            format(Quad_value[3],scitific=T,digits=4),'*',Var[i],"^2",sep="")
        }
        if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),'*',Var[i],'+',
                            format(Quad_value[3],scitific=T,digits=4),'*',Var[i],"^2",sep="")
        }
        if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),'*',Var[i],
                            format(Quad_value[3],scitific=T,digits=4),'*',Var[i],"^2",sep="")
        }
      }
      # If "CP", x ---> x + pmax(x-tau,0)
      if ( p1.bestModels[Var[i],Resp] == "CP" ) {
        tau <- p1.allModels[Var[i],Resp,'CP'][[1]]$psi[2]
        temp_column <- pmax(x1[,Var[i]]-tau,0)
        x1 <- cbind(x1,temp_column)
        names(x1)[ncol(x1)] <- paste(Var[i],"_c",sep="")
        
        Quad_value <- coef(p1.allModels[Var[i],Resp,'CP'][[1]])
        names(trans1)[i] <- Var[i]
        if(tau > 0){
          trans1[i] <- paste(Var[i],"_c=(", Var[i],"-",format(tau,scitific=T,digits=4),")_c",sep=""  )
          if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],'+',
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],"-",format(tau,scitific=T,digits=4),")_c",sep="")
          }
          if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],"-",format(tau,scitific=T,digits=4),")_c",sep="")
          }
          if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],'+',
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],"-",format(tau,scitific=T,digits=4),")_c",sep="")
          }
          if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],"-",format(tau,scitific=T,digits=4),")_c",sep="")
          }
        }else{
          trans1[i] <- paste(Var[i],"_c=(", Var[i],format(tau,scitific=T,digits=4),")_c",sep=""  )
          if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],'+',
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],format(tau,scitific=T,digits=4),")_c",sep="")
          }
          if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],format(tau,scitific=T,digits=4),")_c",sep="")
          }
          if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],'+',
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],format(tau,scitific=T,digits=4),")_c",sep="")
          }
          if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
            trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                              format(Quad_value[2],scitific=T,digits=4),'*',Var[i],
                              format(Quad_value[3],scitific=T,digits=4),'*(',Var[i],format(tau,scitific=T,digits=4),")_c",sep="")
          }
        }
        
      }
      # If "nls", x ---> exp(cx)
      # Nonlinear least squares (nls)
      if ( p1.bestModels[Var[i],Resp] == "nls" ) {
        Quad_value <- coef(p1.allModels[Var[i],Resp,'nls'][[1]])
        x1[,Var[i]] <- exp(Quad_value[3]*x1[,Var[i]])
        trans1[i] <- paste(Var[i],"_t=exp(", format(Quad_value[3],scitific=T,digits=4),"*",Var[i],")",sep="")
        names(trans1)[i] <- Var[i]
        names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
        if ( Quad_value[2]>0 ) {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),"+",
                            format(Quad_value[2],scitific=T,digits=4),"*exp(",
                            format(Quad_value[3],scitific=T,digits=4),Var[i],")",sep="")
        } else {
          trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=4),
                            format(Quad_value[2],scitific=T,digits=4),"*exp(",
                            format(Quad_value[3],scitific=T,digits=4),Var[i],")",sep="")
        }
      }
      # If "-1", No model
      if ( p1.bestModels[Var[i],Resp] == "-1" ) {
        trans1[i] <- NA
        names(trans1)[i] <- Var[i]
      }
      
    }
    ################## END: APPLY TRANSFORMATION TO PREDICTORS ###################
    
    Var_char <- paste(names(x1)[-1], collapse = "+") # Create equation to put into stepAIC
    Rel <- paste0(Resp, "~", Var_char) # Make complete equation with a response variable
    lm.full.model <- do.call("lm", list(Rel, data=as.name("x1"))) # Create a linear regression model with Rel equation
    # stepAIC() checks for all possible variable combinations and selects the model with fewer terms and high adjusted R^2
    if (criterion == "AIC") {
      lm.best <-
        stepAIC(lm.full.model, direction = "backward", trace = FALSE)
    } else  {
      lm.best <-
        stepAIC(
          lm.full.model,
          direction = "backward",
          trace = FALSE,
          k = log(nrow(x))
        )
    }
    
    R2 <- summary(lm.best)$r.squared
    adj.R2 <- summary(lm.best)$adj.r.squared
    best.vars <- names(coef(lm.best))[-1] # Best variables from stepAIC()
    
    coeff_number <- round(lm.best$coeff,digits=4)
    names(coeff_number) <- NULL
    coeff_names <- names(lm.best$coeff)
    lm.best.text <- paste(Resp,"=",sep="") # lm.best.text generates the same equation with corresponding coefficient numbers
    for ( ii in 1:length(coeff_names) ) {
      if (as.numeric(gsub(" ","",coeff_number[ii]))<0) {
        if (coeff_names[ii]=='(Intercept)') {
          lm.best.text <- paste(lm.best.text,coeff_number[ii],sep="")
        } else {
          if (length(grep("_t",coeff_names[ii]))+length(grep("__2",coeff_names[ii]))+length(grep("_c",coeff_names[ii]))==0) {
            lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',coeff_names[ii],sep="")
          }
          if (length(grep("_t",coeff_names[ii]))!=0) {
            names_temp <- substr(coeff_names[ii],1,regexpr("_t",coeff_names[ii])[1]-1) # Removes _t
            names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
            lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',names_temp2,sep="")
          }
          if (length(grep("_c",coeff_names[ii]))!=0) {
            names_temp <- substr(coeff_names[ii],1,regexpr("_c",coeff_names[ii])[1]-1) # Removes _c
            names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
            lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',names_temp2,sep="")
          }
          if (length(grep("__2",coeff_names[ii]))!=0) { # Remove __2 and substitute with ^2
            names_temp <- gsub("__","^",coeff_names[ii])
            lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',names_temp,sep="")
          }
        }
      } else {
        if (coeff_names[ii]=='(Intercept)') {
          lm.best.text <- paste(lm.best.text,gsub(" ","",coeff_number[ii]),sep="")
        } else {
          if (length(grep("_t",coeff_names[ii]))+length(grep("__2",coeff_names[ii]))+length(grep("_c",coeff_names[ii]))==0) {
            lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',coeff_names[ii],sep="")
          }
          if (length(grep("_t",coeff_names[ii]))!=0) {
            names_temp <- substr(coeff_names[ii],1,regexpr("_t",coeff_names[ii])[1]-1)
            names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
            lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',names_temp2,sep="")
          }
          if (length(grep("_c",coeff_names[ii]))!=0) {
            names_temp <- substr(coeff_names[ii],1,regexpr("_c",coeff_names[ii])[1]-1)
            names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
            lm.best.text <- paste(lm.best.text,'+',coeff_number[ii],'*',names_temp2,sep="")
          }
          if (length(grep("__2",coeff_names[ii]))!=0) {
            names_temp <- gsub("__","^",coeff_names[ii])
            lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',names_temp,sep="")
          }
        }
      }
    }
    
    # Make a new row in the Res.print table for each of the predictor variables in the additive model
    Res.print.local <- matrix(NA, nrow = 1, ncol = nRes) # Make annother matrix to save Principl 2 results (from stepAIC best models)
    colnames(Res.print.local) <- c("endogenous",  "Variable",  "GModel", 
                                   "GAdj-R2",       "IModel",  "IAdj-R2",
                                   "Gcoeff",       "GModelB",  "Model" )
    for (n in 1:length(best.vars)){
      Res.print.newrow <- matrix(NA, nrow = 1, ncol = nRes) # Temporary row that will add rows to Res.print and Res.print.local each loop
      Res.print.newrow <- as.data.frame(Res.print.newrow)
      colnames(Res.print.newrow) <- c("endogenous",  "Variable", "GModel", 
                                      "GAdj-R2",     "IModel",   "IAdj-R2",
                                      "Gcoeff",       "GModelB", "Model" )
      Res.print.newrow[1,1] <- Resp # 1st column is Endogenous
      #Res.print.newrow[1,3] <- format(adj.R2,scitific=T,digits=2)		
      #Res.print.newrow[1,4] <- paste(round(summary(lm.best)$coeff[,4],digits=3),collapse="  ")		
      #Res.print.newrow[1,4] <- paste(round(R2,digits=3),round(adj.R2,digits=3),sep="  ")			
      Res.print.newrow[1,4] <- format(adj.R2,scitific=T,digits=4) # 4th column is GAdj-R2
      Res.print.newrow[1,7] <- paste(format(lm.best$coeff,digits=4),collapse=" ") # 7th column is Gcoeff 
      best.vars.rep <- best.vars # Principle 2 results from stepAIC()
      for (iitemp in 1:length(best.vars)) {
        if (regexpr("_t",best.vars[iitemp])!=-1) {
          temp_str <- trans1[substr(best.vars[iitemp],1,regexpr("_t",best.vars[iitemp])-1)]
          best.vars.rep[iitemp] <- substr(temp_str,regexpr("=",temp_str)+1,nchar(temp_str))
        }
      }
      Res.print.newrow[1,8] <- paste(Resp,"=",gsub("__2","^2",paste(best.vars.rep,collapse="+")),sep="") # 8th column is GModelB (no coefficients)
      best.vars.rep2 <- best.vars
      for (j in 1:length(best.vars)) {
        if ((length(grep("__2",best.vars[j]))==0)&(length(grep("_t",best.vars[j]))==0)&(length(grep("_c",best.vars[j]))==0)) {
          best.vars.rep2[j] <- best.vars[j]
        }else{
          if (length(grep("_t",best.vars[j]))!=0) {
            new_string <- substr(best.vars[j],1,regexpr("_t",best.vars[j])[1]-1)
          } 
          if (length(grep("__2",best.vars[j]))!=0) {
            new_string <- substr(best.vars[j],1,regexpr("__2",best.vars[j])[1]-1)
          } 
          if (length(grep("_c",best.vars[j]))!=0) {
            new_string <- substr(best.vars[j],1,regexpr("_c",best.vars[j])[1]-1)
          }
          best.vars.rep2[j] <- new_string
        }
      }   
      Res.print.newrow[1,9] <- paste(Resp,"=f(",paste(unique(best.vars.rep2),collapse=","),")",sep="") # Last column is Model: Resp=f(Vars)
      if ((length(grep("__2",best.vars[n]))==0)&(length(grep("_t",best.vars[n]))==0)&(length(grep("_c",best.vars[n]))==0)) {
        Res.print.newrow[1,2]<- best.vars[n] # 2nd column is Variable
        Res.print.newrow[1,3] <- lm.best.text # 3rd column is GModel
        Res.print.newrow[1,5] <- trans[best.vars[n]] # 5th column is IModel
        if (p1.bestModels[best.vars[n],Resp]!='-1') {
          Individual_summary <- summary(p1.allModels[best.vars[n],Resp,p1.bestModels[best.vars[n],Resp]][[1]])
          #Res.print.newrow[1,6] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
          if (p1.bestModels[best.vars[n],Resp]!="nls") {
            #Res.print.newrow[1,7] <- paste(round(Individual_summary$r.squared,digits=3),
            #                               round(Individual_summary$adj.r.squared,digits=3), sep="  ")
            Res.print.newrow[1,6] <- round(Individual_summary$adj.r.squared,digits=3)
          }
        }
        #Res.print.newrow[1,7] <- which(order(summary(lm.best)$coeff[-1,4])==n)
        Res.print <<- rbind(Res.print, Res.print.newrow) # Add new rows to Res.print
        Res.print.local <- rbind(Res.print.local,Res.print.newrow) # Add new rows to Res.print.local (same as Res.print)
      } else {
        if (length(grep("_t",best.vars[n]))!=0) {
          new_string <- substr(best.vars[n],1,regexpr("_t",best.vars[n])[1]-1)
        } 
        if (length(grep("__2",best.vars[n]))!=0) {
          new_string <- substr(best.vars[n],1,regexpr("__2",best.vars[n])[1]-1)
        } 
        if (length(grep("_c",best.vars[n]))!=0) {
          new_string <- substr(best.vars[n],1,regexpr("_c",best.vars[n])[1]-1)
        }
        Res.print.newrow[1,2]<- new_string # 2nd column is Variable
        Res.print.newrow[1,3] <- lm.best.text # 3rd column is GModel
        Res.print.newrow[1,5] <- trans[new_string] # 5th column is IModel
        Individual_summary <- summary(p1.allModels[new_string,Resp,p1.bestModels[new_string,Resp]][[1]])
        #Res.print.newrow[1,6] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
        if (p1.bestModels[new_string,Resp]!="nls") {
          #Res.print.newrow[1,7] <- paste(round(Individual_summary$r.squared,digits=3),
          #                               round(Individual_summary$adj.r.squared,digits=3), sep="  ")
          Res.print.newrow[1,6] <- round(Individual_summary$adj.r.squared,digits=3)
        }
        #Res.print.newrow[1,7] <- which(order(summary(lm.best)$coeff[-1,4])==n)
        insert_ind <- 0
        Res.temp <- Res.print[-1,] # Excludes the first row (NA's) of Res.print
        if (nrow(Res.temp)>0) {
          for ( mm in 1:nrow(Res.temp) ) {
            if ( Res.temp[mm,1]==Resp ) {
              if ( Res.temp[mm,2]==new_string ) {
                insert_ind <- 1
              }
            }
          }
        }
        if (insert_ind==0) {
          Res.print <<- rbind(Res.print, Res.print.newrow) # Add new rows to Res.print
          Res.print.local <- rbind(Res.print.local,Res.print.newrow) # Add new rows to Res.print.local (same as Res.print)
        }
      }	
    }
    all_var <- colnames(x[-1])
    status_matrix[i_status,Resp] <- 2 # Resp was originally 1, changes to 2
    status_matrix0 <- matrix(NA,nrow=1,ncol=dim(status_matrix)[2])
    colnames(status_matrix0) <- colnames(status_matrix) # status_matrix0 is first filled with NA's
    for ( sig_ind in 1:length(all_var) ) { # sig_ind starts from 1
      if ( (all_var[sig_ind]!=colnames(status_matrix)[1])&(max(status_matrix[,all_var[sig_ind]])==0) ) {
        status_temp <- matrix(NA,nrow=1,ncol=dim(status_matrix)[2])
        colnames(status_temp) <- colnames(status_matrix)
        status_temp[1,] <- status_matrix[i_status,] # Fill in the 1st row of status_temp
        status_temp[1,all_var[sig_ind]] <- 1 # Columns with all_var are filled with 1
        new_sig <- all_var[-sig_ind] 
        if (length(new_sig)!=0) {
          for ( i_new in 1:length(new_sig) ) {
            if (max(status_matrix[,new_sig[i_new]])!=0) {
              status_temp[1,new_sig[i_new]] <- 2
            }
          }
        }
        status_temp[1,1] <- 0 # Make sure the 1st column is 0
        status_matrix0 <- rbind(status_matrix0,status_temp) # Add status_temp to status_matrix0. status_matrix0 will later be merged with status_matrix
      }
    }
    status_matrix <<- rbind(status_matrix,status_matrix0[-1,]) # Add rows to status_matrix
  }  
    # The 1st round when i_status is 1 only considers the main response
    # Then adds rows so that other variables (all_var) are also considered
    # As it finishes one row and goes to the next row, Resp changes from 1 to 2
  
  
  ###############################
  ## Main scripts; 'Multiple.relation' functions need to be called
  ############################### 
  
  # Apply principle 1 on the data to find the best models
  if ( !is.null(exogenous) ) {
    exogenous.p1 <- exogenous
  } else {
    exogenous.p1 <- 1
  }
  if ( !is.null(endogenous) ) {
    endogenous.p1 <- endogenous
  } else {
    endogenous.p1 <- 2
  }
  
  # Check for strength type problem
  
  if (str) {
    p1.result <- netSEMp1(x, str = TRUE) # Principle 1 results. Consider strength type problem
  } else {
    p1.result <- netSEMp1(x) # Don't consider strength type problem
  }
  
  p1.bestModels <- p1.result$bestModels
  p1.allModels <- p1.result$allModels
  
  # End apply principle 1
  
  nVar <- ncol(x)    # Total number of variables
  nRVar <- nVar - 1  # Total number of possible endogenous variables
  
  nRes <- 9  # Number of cells in print variable; see below
  Res.print <- matrix(NA, nrow = 1, ncol = nRes) # Output for Principle 2
  colnames(Res.print) <- c("endogenous",  "Variable",  "GModel", 
                           "GAdj-R2",       "IModel",  "IAdj-R2",
                           "Gcoeff",       "GModelB", "Model" )
  
  status_matrix <- matrix(0, nrow=1, ncol=ncol(x))
  colnames(status_matrix) <- colnames(x)[c(2,1,3:ncol(x))] # Response variable comes to the 1st column
  status_matrix[1,2] <- 1
  i_status <- 1
  
  
  
  
  while (i_status<=nrow(status_matrix)) { # Each row of status_matrix provides new status_index to run Multiple.relation (assign new Resp and Var)
    Multiple.relation(x) # Need to keep this as separate script
    i_status <- i_status + 1
  }
  
  Res.print <- Res.print[-1,]
  res <- list(data = x, res.print = Res.print) # Outputs are three tables
  class(res) <- c("netSEMp2","list")
  invisible(res)
}
