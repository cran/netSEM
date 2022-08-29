##' Plot netSEMp2 result
##' plot.netSEMp2 plots a network structural equation model diagram, fitted under principle 2, based on best functional form for each selected pairwise variable.
##' 
##' @title Plotting of netSEM diagram
##' 
##' @import DiagrammeRsvg
##' @import magrittr
##' @import svglite
##' @import rsvg
##' @import png
##' @importFrom DiagrammeR grViz
##' 
##' @param x An object of class "netSEMp2", the returned list from \code{netSEMp2}.
##' @param cutoff A threshold value for adjusted R-squared. The maximum number of cutoff is 3.
##' @param latent The latent variable that corresponds to the mechanic variable. The default is NULL.
##' @param plot.save True/False, it saves the network diagram plot as png file. The default is false.
##' @param filename A character string naming a file to save as a png file.
##' @param style True/False, it plots the first interval in the network diagram with dotted weak line. The default is True.
##' @param ... A S3 generic/method consistency.
##' 
##' @return An html style plot of multiple regression relationship pathway diagram between exogenous variables and an endogenous variable. 
##' Arrows show relationships between each variable with given statistical relations along the connection lines.
##' 
##' @export
##'
##' @seealso \link[netSEM]{netSEMp2}
##'
##' @examples
##' # Load acrylic data set
##' data(acrylic)
##' 
##' # Build a netSEMp2 model
##' ans2 <- netSEMp2(acrylic, criterion = "AIC")
##' ans2_BIC <- netSEMp2(acrylic, criterion = "BIC")
##' 
##' # Plot the network model 
##' plot(ans2, cutoff = c(0.3,0.6,0.8))
##' plot(ans2_BIC , cutoff = c(0.3,0.6,0.8))
##' 
##' \dontrun{
##' # Drop Relationship lower than minimum cutoff
##' plot(ans2, cutoff = c(0.3,0.6,0.8), style = FALSE)
##' 
##' # plot network model with latent argument labels
##' plot(ans2, cutoff = c(0.3, 0.6, 0.8), 
##'      latent = c('IAD1' = 'FundAbsEdge', 
##'                 'IAD2' = 'UVStab', 
##'                 'IAD2p' = 'UVStab', 
##'                 'IAD3' = 'YelMet'))
##' 
##' # plot the network diagram and save file
##' #plot(ans, cutoff = c(0.3,0.6,0.8), plot.save = TRUE, filename = "acrylic-netSEMp2"))
##' }

plot.netSEMp2 <- function(x, cutoff = c(0.3,0.6,0.8),latent = NULL, plot.save = FALSE, 
                          filename = NULL, style = TRUE, ...) {
 
  # Create empty data frames for storing information for each cutoff region --- 
  rtp1 <- x$res.print

  rtp1.a <- rtp1[-c(1:nrow(rtp1)),]
  rtp1.b <- rtp1[-c(1:nrow(rtp1)),]
  rtp1.c <- rtp1[-c(1:nrow(rtp1)),]
  rtp1.d <- rtp1[-c(1:nrow(rtp1)),]
  
  # Set up Cutoff Regions -----------------------------------------------------
  count <- length(cutoff)
  ## Based on the number of cutoff, to separate different regions
  ## Number of cutoff = 1, have two regions
  if (count == 1) {
    c1 <- cutoff[1]  
    rtp1.a <- rtp1[rtp1[, "GAdj-R2"] <  c1 , ]
    rtp1.b <- rtp1[rtp1[, "GAdj-R2"] >= c1 , ]
  }
  ## Number of cutoff = 2, have three regions
  if (count == 2) {
    c1 <- cutoff[1] ; c2 <- cutoff[2] 
    rtp1.a <- rtp1[rtp1[, "GAdj-R2"] < c1, ]
    rtp1.b <- rtp1[rtp1[, "GAdj-R2"] >= c1 & rtp1[, "GAdj-R2"] < c2, ]
    rtp1.c <- rtp1[rtp1[, "GAdj-R2"] >= c2 , ]
  }
  ## Number of cutoff = 3, have four regions
  if (count == 3) {
    c1 <- cutoff[1] ; c2 <- cutoff[2] ; c3 <- cutoff[3] 
    rtp1.a <- rtp1[rtp1[, "GAdj-R2"] < c1, ]
    rtp1.b <- rtp1[rtp1[, "GAdj-R2"] >= c1 & rtp1[, "GAdj-R2"] < c2, ]
    rtp1.c <- rtp1[rtp1[, "GAdj-R2"] >= c2 & rtp1[, "GAdj-R2"] < c3, ]
    rtp1.d <- rtp1[rtp1[, "GAdj-R2"] >= c3 , ]
  }
  
  ## Names of Resp and Vars
  name <- colnames(x$data)
  conp1.m1 <- sapply(3:length(name), function(i) { 
    paste0(name[i], "[fillcolor=khaki]") 
  } 
  )
  
  
  ## Assign the latent variable to each mechanic variable
  if (!is.null(latent)) {
    conp1.La1 <- sapply(1:length(latent), function(i) { 
      paste0(latent[i], "[fillcolor=LightBlue]") 
    } 
    )
    
    conp1.L1 <- sapply(1:length(latent), function(i) {
      paste0(name[i+2], "->", latent[i], "[ fontsize=80, 
             arrowhead=none]")
    }
    )
    ## Style is True represents plot the first interval with dotted weak line
    if (style) {
      ## Set up parameters for the first region
      if (dim(rtp1.a)[1] > 0) {
        conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
          paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@", i, "', 
                 fontsize=80, style=dotted]")
        }
        )
        conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
          paste0("[", i, "]", ": ", "paste('", 
                 #colnames(rtp1.a)[3], ": ", rtp1.a[i,3], "\\n", 
                 #colnames(rtp1.a)[4], ": ", rtp1.a[i,4], "\\n", 
                 #colnames(rtp1.a)[5], ": ", rtp1.a[i,5], "\\n", 
                 #colnames(rtp1.a)[9], ": ", 
                 rtp1.a[i,9],   
                 "')"
          )
        }
        )
      }
      ## Set up parameters for the second region
      count.1 <- nrow(rtp1.a)
      if (dim(rtp1.b)[1] > 0) {
        conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
          paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, 
                 "', fontsize=80]")
        }
        )
        conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
          paste0("[", i + count.1, "]", ": ", "paste('", 
                 #colnames(rtp1.b)[3], ": ", rtp1.b[i,3], "\\n", 
                 #colnames(rtp1.b)[4], ": ", rtp1.b[i,4], "\\n", 
                 #colnames(rtp1.b)[5], ": ", rtp1.b[i,5], "\\n", 
                 #colnames(rtp1.b)[9], ": ", 
                 rtp1.b[i,9],               
                 "')"
          )
        }
        )
      }
      ## Set up parameters for the third region
      count.2 <- count.1 + nrow(rtp1.b)
      if (dim(rtp1.c)[1] > 0) {
        conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
          paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@" ,i + count.2, 
                 "', fontsize=80]")
        }
        )
        conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
          paste0("[", i + count.2, "]", ": ", "paste('", 
                 #colnames(rtp1.c)[3], ": ", rtp1.c[i,3], "\\n", 
                 #colnames(rtp1.c)[4], ": ", rtp1.c[i,4], "\\n", 
                 #colnames(rtp1.c)[5], ": ", rtp1.c[i,5], "\\n", 
                 #colnames(rtp1.c)[9], ": ", 
                 rtp1.c[i,9],                
                 "')"
          )
        }
        )
      }
      ## Set up parameters for the four region
      count.3 <- count.2 + nrow(rtp1.c)
      if (dim(rtp1.d)[1] > 0) {
        conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
          paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, 
                 "', fontsize=80]")
        }
        )
        conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
          paste0("[", i + count.3, "]", ": ", "paste('", 
                 #colnames(rtp1.d)[3], ": ", rtp1.d[i,3], "\\n", 
                 #colnames(rtp1.d)[4], ": ", rtp1.d[i,4], "\\n", 
                 #colnames(rtp1.d)[5], ": ", rtp1.d[i,5], "\\n", 
                 #colnames(rtp1.d)[9], ": ", 
                 rtp1.d[i,9], 
                 "')"
          )
        }
        )
      }
      
      ## Test whether each region has model/parameter or not
      A <- rep(0,4)
      if (exists("conp1.a1") == T) { A[1] <- 1 }
      if (exists("conp1.b1") == T) { A[2] <- 1 }
      if (exists("conp1.c1") == T) { A[3] <- 1 }
      if (exists("conp1.d1") == T) { A[4] <- 1 }
      
      ## All of four regions have models
      ### Plotting starts here. Bug could start here ###
      
      if (sum(A) == 4 ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.a1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                         paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n", 
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## Three regions have models
      if (sum(A) == 3 ) {
        ## Only region 'a' does not have models, plot other three regions
        if(exists("conp1.a1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.b2, collapse = "\n"), "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n", 
                           paste(conp1.d2, collapse = "\n"), "\n"
          )
        }
        ## Only region 'b' does not have models, plot other three regions
        if (exists("conp1.b1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2, collapse = "\n"), "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n", 
                           paste(conp1.d2, collapse = "\n"), "\n" 
          )
        }
        ## Only region 'c' does not have models, plot other three regions
        if (exists("conp1.c1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n",
                           paste(conp1.a2, collapse = "\n"), "\n", 
                           paste(conp1.b2, collapse = "\n"), "\n", 
                           paste(conp1.d2, collapse = "\n"), "\n" 
          )
        }
        ## Only region 'd' does not have models, plot other three regions
        if (exists("conp1.d1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2, collapse = "\n"), "\n", 
                           paste(conp1.b2, collapse = "\n"), "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n" 
          )
        }
      }
      ## Two regions have models
      if (sum(A) == 2 ) {
        ## Only region 'a' and 'b' does not have models, plot other two regions
        if (exists("conp1.a1") == F & 
            exists("conp1.b1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n", 
                           paste(conp1.d2,collapse = "\n"), "\n" 
          )
        }
        ## Only region 'a' and 'c' does not have models, plot other two regions
        if (exists("conp1.a1") == F & 
            exists("conp1.c1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}","\n",
                           paste(conp1.b2, collapse = "\n"), "\n", 
                           paste(conp1.d2, collapse = "\n"), "\n" 
          )
        }
        ## Only region 'a' and 'd' does not have models, plot other two regions
        if (exists("conp1.a1") == F & 
            exists("conp1.d1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.b2, collapse = "\n"), "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n" 
          )
        }
        ## Only region 'b' and 'c' does not have models, plot other two regions
        if (exists("conp1.b1") == F & 
            exists("conp1.c1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2, collapse = "\n"), "\n", 
                           paste(conp1.d2, collapse = "\n"), "\n" 
          )
        }
        ## Only region 'b' and 'd' does not have models, plot other two regions
        if (exists("conp1.b1") == F & 
            exists("conp1.d1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2, collapse = "\n"), "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n" 
          )
        }
        ## Only region 'c' and 'd' does not have models, plot other two regions
        if (exists("conp1.c1") == F & 
            exists("conp1.d1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2, collapse = "\n"), "\n", 
                           paste(conp1.b2, collapse = "\n"), "\n" 
          )
        }
      }
      ## Only one region has the model
      if (sum(A) == 1 ) {
        ## Region 'a', 'b' and 'c' does not have models, plot the other region
        if (exists("conp1.a1") == F & 
            exists("conp1.b1") == F & 
            exists("conp1.c1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.d2, collapse = "\n"), "\n" 
          )
        }
        ## Region 'a', 'b' and 'd' does not have models, plot the other region
        if (exists("conp1.a1") == F & 
            exists("conp1.b1") == F & 
            exists("conp1.d1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4,
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n" 
          )
        }
        ## Region 'a', 'c' and 'd' does not have models, plot the other region
        if (exists("conp1.a1") == F & 
            exists("conp1.c1") == F & 
            exists("conp1.d1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T,
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.b2, collapse = "\n"), "\n" 
          )
        }
        ## Region 'b', 'c' and 'd' does not have models, plot the other region
        if (exists("conp1.b1") == F & 
            exists("conp1.c1") == F & 
            exists("conp1.d1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2, collapse = "\n"), "\n" 
          )
        }
      }
      
  } else {
    
    ## Set up parameters for the first region
    if (dim(rtp1.a)[1] > 0) {
      conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
        paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@", i, "', 
               fontsize=80, style=dotted]")
      }
      )
      # conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
      #  paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",rtp1.a[i,3],"\\n",
      #         colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
      #}
      #)
    }
    ## Set up parameters for the second region
    count.1 <- 0
    if (dim(rtp1.b)[1] > 0) {
      conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
        paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, "', 
               fontsize=80]")
      }
      )
      conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
        paste0("[", i + count.1, "]", ": ", "paste('", 
               #colnames(rtp1.b)[3], ": ", rtp1.b[i,3], "\\n", 
               #colnames(rtp1.b)[4], ": ", rtp1.b[i,4], "\\n", 
               #colnames(rtp1.b)[5], ": ", rtp1.b[i,5], "\\n", 
               #colnames(rtp1.b)[9], ": ", 
               rtp1.b[i,9], 
               "')"
        )
      }
      )
    }
    ## Set up parameters for the third region
    count.2 <- count.1 + nrow(rtp1.b)
    if (dim(rtp1.c)[1] > 0) {
      conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
        paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@", i + count.2, "', 
               fontsize=80]")
      }
      )
      conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
        paste0("[", i + count.2, "]", ": ", "paste('", 
               #colnames(rtp1.c)[3], ": ", rtp1.c[i,3], "\\n", 
               #colnames(rtp1.c)[4], ": ", rtp1.c[i,4], "\\n", 
               #colnames(rtp1.c)[5], ": ", rtp1.c[i,5], "\\n", 
               #colnames(rtp1.c)[9], ": ",
               rtp1.c[i,9], 
               "')"
        )
      }
      )
    }
    ## Set up parameters for the four region
    count.3 <- count.2 + nrow(rtp1.c)
    if (dim(rtp1.d)[1] > 0) {
      conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
        paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, "', 
               fontsize=80]")
      }
      )
      conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
        paste0("[", i + count.3, "]", ": ", "paste('", 
               #colnames(rtp1.d)[3], ": ", rtp1.d[i,3], "\\n", 
               #colnames(rtp1.d)[4], ": ", rtp1.d[i,4], "\\n", 
               #colnames(rtp1.d)[5], ": ", rtp1.d[i,5], "\\n", 
               #colnames(rtp1.d)[9], ": ", 
               rtp1.d[i,9], 
               "')"
        )
      }
      )
    }
    
    ## Test whether each region has model/parameter or not
    A <- rep(0,4)
    if (exists("conp1.a1") == T) { A[1] <- 1 }
    if (exists("conp1.b1") == T) { A[2] <- 1 }
    if (exists("conp1.c1") == T) { A[3] <- 1 }
    if (exists("conp1.d1") == T) { A[4] <- 1 }
    
    ## All of four regions have models
    if (sum(A) == 4 ) {
      c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                       nodesep = 1.8, ranksep =0.9, layout = dot, 
                       rankdir = LR]", "\n", 
                       "node[fontname=Helvetica, shape=box, fixedsize=T, 
                       fontsize=100, color=dodgerblue, style=filled, penwidth=6, 
                       width=4, height=3]", "\n", 
                       paste0(name[1]), "[fillcolor=violet]", "\n", 
                       "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                       fontsize=100, color=dodgerblue, style=filled, 
                       penwidth=6]", "\n", 
                       paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                       "node[fontname=Helvetica, shape=box, fixedsize=F, 
                       color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                       paste0(conp1.m1, collapse = ";"), "\n",
                       "node[fontname=Helvetica, shape=square, fixedsize=F, 
                       color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                       paste0(conp1.La1,collapse = ";"), "\n", 
                       "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                       paste0(conp1.L1, collapse = "\n"), "\n",
                       # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                       # paste0(conp1.a1, collapse = "\n"), "\n", 
                       "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                       paste0(conp1.b1, collapse = "\n"), "\n", 
                       "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                       paste0(conp1.c1, collapse = "\n"), "\n", 
                       "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                       paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                       # paste(conp1.a2,collapse="\n") ,"\n",
                       paste(conp1.b2, collapse = "\n"), "\n", 
                       paste(conp1.c2, collapse = "\n"), "\n", 
                       paste(conp1.d2, collapse = "\n"), "\n" 
      )
      
    }
    ## Three regions have models
    if (sum(A) == 3 ) {
      ## only region 'a' does not have models, plot other three regions
      if (exists("conp1.a1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4,
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n", 
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## Only region 'b' does not have models, plot other three regions
      if (exists("conp1.b1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]","\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, 
                         width=4,height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1,collapse="\n"),"\n",
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                         # paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n", 
                         paste(conp1.d2, collapse = "\n"), "\n"
        )
      }
      ## Only region 'c' does not have models, plot other three regions
      if (exists("conp1.c1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                         # paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n", 
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## Only region 'd' does not have models, plot other three regions
      if (exists("conp1.d1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                         # paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n" 
        )
      }
    }
    ## Two regions have models
    if (sum(A) == 2 ) {
      ## Only region 'a' and 'b' does not have models, plot other two regions
      if (exists("conp1.a1") == F & 
          exists("conp1.b1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n", 
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## Only region 'a' and 'c' does not have models, plot other two regions
      if (exists("conp1.a1") == F & 
          exists("conp1.c1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n", 
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## Only region 'a' and 'd' does not have models, plot other two regions
      if (exists("conp1.a1") == F & 
          exists("conp1.d1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T,
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse=";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n" 
        )
      }
      ## Only region 'b' and 'c' does not have models, plot other two regions
      if (exists("conp1.b1") == F & 
          exists("conp1.c1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                         # paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## Only region 'b' and 'd' does not have models, plot other two regions
      if (exists("conp1.b1") == F & 
          exists("conp1.d1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                         #paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n" 
        )
      }
      ## Only region 'c' and 'd' does not have models, plot other two regions
      if (exists("conp1.c1") == F & 
          exists("conp1.d1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse = "\n"), "\n", 
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                         # paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n" 
        )
      }
    }
    ## Only one region has the model
    if (sum(A) == 1 ) {
      ## Region 'a', 'b' and 'c' does not have models, plot the other region
      if (exists("conp1.a1") == F & 
          exists("conp1.b1") == F & 
          exists("conp1.c1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                         paste0(conp1.d1, collapse = "\n"), "\n", "}","\n",
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## Region 'a', 'b' and 'd' does not have models, plot the other region
      if (exists("conp1.a1") == F & 
          exists("conp1.b1") == F & 
          exists("conp1.d1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n" 
        )
      }
      ## Region 'a', 'c' and 'd' does not have models, plot the other region
      if (exists("conp1.a1") == F & 
          exists("conp1.c1") == F & 
          exists("conp1.d1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n" 
        )
      }
      ## Region 'b', 'c' and 'd' does not have models, plot the other region
      if (exists("conp1.b1") == F & 
          exists("conp1.c1") == F & 
          exists("conp1.d1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=box, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=violet]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n", "}",
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n", "\n" 
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse = "\n"), "\n", 
                         # paste(conp1.a2, collapse = "\n"), "\n" 
        )
    }
    }
    
    }
    
    
    } 
    else{
      
      ## Style is True represents plot the first interval with dotted weak line
      if (style) {
        ## Set up parameters for the first region
        if (dim(rtp1.a)[1] > 0) {
          conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
            paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@", i, "', 
                   fontsize=80, style=dotted]")
          }
          )
          conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
            paste0("[", i, "]", ": ", "paste('", 
                   #colnames(rtp1.a)[3], ": ", rtp1.a[i,3], "\\n", 
                   #colnames(rtp1.a)[4], ": ", rtp1.a[i,4], "\\n", 
                   #colnames(rtp1.a)[5], ": ", rtp1.a[i,5], "\\n", 
                   #colnames(rtp1.a)[9], ": ", 
                   rtp1.a[i,9],   
                   "')"
            )
          }
          )
        }
        ## Set up parameters for the second region
        count.1 <- nrow(rtp1.a)
        if (dim(rtp1.b)[1] > 0) {
          conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
            paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, 
                   "', fontsize=80]")
          }
          )
          conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
            paste0("[", i + count.1, "]", ": ", "paste('", 
                   #colnames(rtp1.b)[3], ": ", rtp1.b[i,3], "\\n", 
                   #colnames(rtp1.b)[4], ": ", rtp1.b[i,4], "\\n", 
                   #colnames(rtp1.b)[5], ": ", rtp1.b[i,5], "\\n", 
                   #colnames(rtp1.b)[9], ": ", 
                   rtp1.b[i,9],               
                   "')"
            )
          }
          )
        }
        ## Set up parameters for the third region
        count.2 <- count.1 + nrow(rtp1.b)
        if (dim(rtp1.c)[1] > 0) {
          conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
            paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@" ,i + count.2, 
                   "', fontsize=80]")
          }
          )
          conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
            paste0("[", i + count.2, "]", ": ", "paste('", 
                   #colnames(rtp1.c)[3], ": ", rtp1.c[i,3], "\\n", 
                   #colnames(rtp1.c)[4], ": ", rtp1.c[i,4], "\\n", 
                   #colnames(rtp1.c)[5], ": ", rtp1.c[i,5], "\\n", 
                   #colnames(rtp1.c)[9], ": ", 
                   rtp1.c[i,9],                
                   "')"
            )
          }
          )
        }
        ## Set up parameters for the four region
        count.3 <- count.2 + nrow(rtp1.c)
        if (dim(rtp1.d)[1] > 0) {
          conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
            paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, 
                   "', fontsize=80]")
          }
          )
          conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
            paste0("[", i + count.3, "]", ": ", "paste('", 
                   #colnames(rtp1.d)[3], ": ", rtp1.d[i,3], "\\n", 
                   #colnames(rtp1.d)[4], ": ", rtp1.d[i,4], "\\n", 
                   #colnames(rtp1.d)[5], ": ", rtp1.d[i,5], "\\n", 
                   #colnames(rtp1.d)[9], ": ", 
                   rtp1.d[i,9], 
                   "')"
            )
          }
          )
        }
        
        ## Test whether each region has model/parameter or not
        A <- rep(0,4)
        if (exists("conp1.a1") == T) { A[1] <- 1 }
        if (exists("conp1.b1") == T) { A[2] <- 1 }
        if (exists("conp1.c1") == T) { A[3] <- 1 }
        if (exists("conp1.d1") == T) { A[4] <- 1 }
        ## All of four regions have models
        if (sum(A) == 4 ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=box, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=violet]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1, collapse = ";"), "\n", 
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.b1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                           paste0(conp1.c1, collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2, collapse = "\n"), "\n", 
                           paste(conp1.b2, collapse = "\n"), "\n", 
                           paste(conp1.c2, collapse = "\n"), "\n", 
                           paste(conp1.d2, collapse = "\n"), "\n" 
          )
        }
        ## Three regions have models
        if (sum(A) == 3 ) {
          ## Only region 'a' does not have models, plot other three regions
          if(exists("conp1.a1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.b2, collapse = "\n"), "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n", 
                             paste(conp1.d2, collapse = "\n"), "\n"
            )
          }
          ## Only region 'b' does not have models, plot other three regions
          if (exists("conp1.b1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.a2, collapse = "\n"), "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n", 
                             paste(conp1.d2, collapse = "\n"), "\n" 
            )
          }
          ## Only region 'c' does not have models, plot other three regions
          if (exists("conp1.c1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n",
                             paste(conp1.a2, collapse = "\n"), "\n", 
                             paste(conp1.b2, collapse = "\n"), "\n", 
                             paste(conp1.d2, collapse = "\n"), "\n" 
            )
          }
          ## Only region 'd' does not have models, plot other three regions
          if (exists("conp1.d1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.a2, collapse = "\n"), "\n", 
                             paste(conp1.b2, collapse = "\n"), "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n" 
            )
          }
        }
        ## Two regions have models
        if (sum(A) == 2 ) {
          ## Only region 'a' and 'b' does not have models, plot other two regions
          if (exists("conp1.a1") == F & 
              exists("conp1.b1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n", 
                             paste(conp1.d2,collapse = "\n"), "\n" 
            )
          }
          ## Only region 'a' and 'c' does not have models, plot other two regions
          if (exists("conp1.a1") == F & 
              exists("conp1.c1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}","\n",
                             paste(conp1.b2, collapse = "\n"), "\n", 
                             paste(conp1.d2, collapse = "\n"), "\n" 
            )
          }
          ## Only region 'a' and 'd' does not have models, plot other two regions
          if (exists("conp1.a1") == F & 
              exists("conp1.d1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.b2, collapse = "\n"), "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n" 
            )
          }
          ## Only region 'b' and 'c' does not have models, plot other two regions
          if (exists("conp1.b1") == F & 
              exists("conp1.c1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.a2, collapse = "\n"), "\n", 
                             paste(conp1.d2, collapse = "\n"), "\n" 
            )
          }
          ## Only region 'b' and 'd' does not have models, plot other two regions
          if (exists("conp1.b1") == F & 
              exists("conp1.d1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.a2, collapse = "\n"), "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n" 
            )
          }
          ## Only region 'c' and 'd' does not have models, plot other two regions
          if (exists("conp1.c1") == F & 
              exists("conp1.d1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.a2, collapse = "\n"), "\n", 
                             paste(conp1.b2, collapse = "\n"), "\n" 
            )
          }
        }
        ## Only one region has the model
        if (sum(A) == 1 ) {
          ## Region 'a', 'b' and 'c' does not have models, plot the other region
          if (exists("conp1.a1") == F & 
              exists("conp1.b1") == F & 
              exists("conp1.c1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.d2, collapse = "\n"), "\n" 
            )
          }
          ## Region 'a', 'b' and 'd' does not have models, plot the other region
          if (exists("conp1.a1") == F & 
              exists("conp1.b1") == F & 
              exists("conp1.d1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4,
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n" 
            )
          }
          ## Region 'a', 'c' and 'd' does not have models, plot the other region
          if (exists("conp1.a1") == F & 
              exists("conp1.c1") == F & 
              exists("conp1.d1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T,
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.b2, collapse = "\n"), "\n" 
            )
          }
          ## Region 'b', 'c' and 'd' does not have models, plot the other region
          if (exists("conp1.b1") == F & 
              exists("conp1.c1") == F & 
              exists("conp1.d1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=box, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1, collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.a2, collapse = "\n"), "\n" 
            )
          }
        }
        
        } else {
          
          ## Set up parameters for the first region
          if (dim(rtp1.a)[1] > 0) {
            conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
              paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@", i, "', 
                     fontsize=80, style=dotted]")
            }
            )
            # conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
            #  paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",rtp1.a[i,3],"\\n",
            #         colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
            #}
            #)
          }
          ## Set up parameters for the second region
          count.1 <- 0
          if (dim(rtp1.b)[1] > 0) {
            conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
              paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, "', 
                     fontsize=80]")
            }
            )
            conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
              paste0("[", i + count.1, "]", ": ", "paste('", 
                     #colnames(rtp1.b)[3], ": ", rtp1.b[i,3], "\\n", 
                     #colnames(rtp1.b)[4], ": ", rtp1.b[i,4], "\\n", 
                     #colnames(rtp1.b)[5], ": ", rtp1.b[i,5], "\\n", 
                     #colnames(rtp1.b)[9], ": ", 
                     rtp1.b[i,9], 
                     "')"
              )
            }
            )
          }
          ## Set up parameters for the third region
          count.2 <- count.1 + nrow(rtp1.b)
          if (dim(rtp1.c)[1] > 0) {
            conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
              paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@", i + count.2, "', 
                     fontsize=80]")
            }
            )
            conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
              paste0("[", i + count.2, "]", ": ", "paste('", 
                     #colnames(rtp1.c)[3], ": ", rtp1.c[i,3], "\\n", 
                     #colnames(rtp1.c)[4], ": ", rtp1.c[i,4], "\\n", 
                     #colnames(rtp1.c)[5], ": ", rtp1.c[i,5], "\\n", 
                     #colnames(rtp1.c)[9], ": ",
                     rtp1.c[i,9], 
                     "')"
              )
            }
            )
          }
          ## Set up parameters for the four region
          count.3 <- count.2 + nrow(rtp1.c)
          if (dim(rtp1.d)[1] > 0) {
            conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
              paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, "', 
                     fontsize=80]")
            }
            )
            conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
              paste0("[", i + count.3, "]", ": ", "paste('", 
                     #colnames(rtp1.d)[3], ": ", rtp1.d[i,3], "\\n", 
                     #colnames(rtp1.d)[4], ": ", rtp1.d[i,4], "\\n", 
                     #colnames(rtp1.d)[5], ": ", rtp1.d[i,5], "\\n", 
                     #colnames(rtp1.d)[9], ": ", 
                     rtp1.d[i,9], 
                     "')"
              )
            }
            )
          }
          
          ## Test whether each region has model/parameter or not
          A <- rep(0,4)
          if (exists("conp1.a1") == T) { A[1] <- 1 }
          if (exists("conp1.b1") == T) { A[2] <- 1 }
          if (exists("conp1.c1") == T) { A[3] <- 1 }
          if (exists("conp1.d1") == T) { A[4] <- 1 }
          
          ## All of four regions have models
          if (sum(A) == 4 ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, penwidth=6, 
                             width=4, height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=violet]", "\n", 
                             "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1, collapse = ";"), "\n", 
                             # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             # paste0(conp1.a1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.b1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                             paste0(conp1.c1, collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                             # paste(conp1.a2,collapse="\n") ,"\n",
                             paste(conp1.b2, collapse = "\n"), "\n", 
                             paste(conp1.c2, collapse = "\n"), "\n", 
                             paste(conp1.d2, collapse = "\n"), "\n" 
            )
            
          }
          ## Three regions have models
          if (sum(A) == 3 ) {
            ## Only region 'a' does not have models, plot other three regions
            if (exists("conp1.a1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4,
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n", 
                               paste(conp1.d2, collapse = "\n"), "\n" 
              )
            }
            ## Only region 'b' does not have models, plot other three regions
            if (exists("conp1.b1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]","\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, 
                               width=4,height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1,collapse="\n"),"\n",
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                               # paste(conp1.a2, collapse = "\n"), "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n", 
                               paste(conp1.d2, collapse = "\n"), "\n"
              )
            }
            ## Only region 'c' does not have models, plot other three regions
            if (exists("conp1.c1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                               # paste(conp1.a2, collapse = "\n"), "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n", 
                               paste(conp1.d2, collapse = "\n"), "\n" 
              )
            }
            ## Only region 'd' does not have models, plot other three regions
            if (exists("conp1.d1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                               # paste(conp1.a2, collapse = "\n"), "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n" 
              )
            }
          }
          ## Two regions have models
          if (sum(A) == 2 ) {
            ## Only region 'a' and 'b' does not have models, plot other two regions
            if (exists("conp1.a1") == F & 
                exists("conp1.b1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n", 
                               paste(conp1.d2, collapse = "\n"), "\n" 
              )
            }
            ## Only region 'a' and 'c' does not have models, plot other two regions
            if (exists("conp1.a1") == F & 
                exists("conp1.c1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n", 
                               paste(conp1.d2, collapse = "\n"), "\n" 
              )
            }
            ## Only region 'a' and 'd' does not have models, plot other two regions
            if (exists("conp1.a1") == F & 
                exists("conp1.d1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T,
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse=";"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n" 
              )
            }
            ## Only region 'b' and 'c' does not have models, plot other two regions
            if (exists("conp1.b1") == F & 
                exists("conp1.c1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                               # paste(conp1.a2, collapse = "\n"), "\n", 
                               paste(conp1.d2, collapse = "\n"), "\n" 
              )
            }
            ## Only region 'b' and 'd' does not have models, plot other two regions
            if (exists("conp1.b1") == F & 
                exists("conp1.d1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                               #paste(conp1.a2, collapse = "\n"), "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n" 
              )
            }
            ## Only region 'c' and 'd' does not have models, plot other two regions
            if (exists("conp1.c1") == F & 
                exists("conp1.d1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                               # paste(conp1.a2, collapse = "\n"), "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n" 
              )
            }
          }
          ## Only one region has the model
          if (sum(A) == 1 ) {
            ## Region 'a', 'b' and 'c' does not have models, plot the other region
            if (exists("conp1.a1") == F & 
                exists("conp1.b1") == F & 
                exists("conp1.c1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}","\n",
                               paste(conp1.d2, collapse = "\n"), "\n" 
              )
            }
            ## Region 'a', 'b' and 'd' does not have models, plot the other region
            if (exists("conp1.a1") == F & 
                exists("conp1.b1") == F & 
                exists("conp1.d1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n" 
              )
            }
            ## Region 'a', 'c' and 'd' does not have models, plot the other region
            if (exists("conp1.a1") == F & 
                exists("conp1.c1") == F & 
                exists("conp1.d1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", "}", "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n" 
              )
            }
            ## Region 'b', 'c' and 'd' does not have models, plot the other region
            if (exists("conp1.b1") == F & 
                exists("conp1.c1") == F & 
                exists("conp1.d1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=box, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=violet]", "\n", 
                               "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", "}", "\n" 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse = "\n"), "\n", 
                               # paste(conp1.a2, collapse = "\n"), "\n" 
              )
          }
            }
          
        }
      }
  
  
  
  
  ## This generates syntax to run "grViz" for plotting using above syntax
  ## "LR" is left to right flow
  semplot <- DiagrammeR::grViz(c.plot)
  
  if (plot.save) {
    ### Export graph as follows
    semplot %>% export_svg %>% charToRaw %>% rsvg %>% 
      png::writePNG(paste0(filename, ".png"))
    
  }
  
  return(semplot)
  }
