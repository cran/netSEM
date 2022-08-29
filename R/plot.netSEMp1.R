##' Plot netSEMp1 result
##' plot.netSEMp1 plots a network structural equation network model diagram, fitted under principle 1, based on best functional form for each selected pairwise variable.
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
##' @param x An object of class "netSEMp1", the returned list from \code{netSEMp1}. Plotting uses the first element of this list (table) in which the first column of it is endogenous variable, second column is variable and other columns are corresponding best functional form, r-squared, adj-r-squared, P-value1, P-value2 and P-value3.
##' @param cutoff A threshold value for adjusted R-squared. The maximum number of cutoff is 3.
##' @param latent The latent variable that corresponds to the mechanic variable. The default is NULL.
##' @param plot.save True/False, it saves the network diagram plot as a png file. The default is false.
##' @param filename A character string naming a file to save as a png file.
##' @param style True/False, it plots the first interval in the network diagram with dotted weak line. The default is True.
##' @param ... A S3 generic/method consistency.
##' 
##' @return An html style plot of pairwise relationship pathway diagram between exogenous variables and an endogenous variable. 
##' Arrows show relationships between each variable with given statistical relations along the connection lines.
##' 
##' @export
##'
##' @seealso \link[netSEM]{netSEMp1}
##'
##' @examples
##' # Load acrylic data set
##' data(acrylic)
##' 
##' # Build a netSEM model
##' ans <- netSEMp1(acrylic)
##' 
##' # Plot the network model 
##' plot(ans,cutoff = c(0.3,0.6,0.8))
##' 
##' # Plot the network diagram with latent argument labels
##' plot(ans, cutoff = c(0.3, 0.6, 0.8), 
##'      latent = c('IAD1' = 'FundAbsEdge', 
##'                 'IAD2' = 'UVStab', 
##'                 'IAD2p' = 'UVStab', 
##'                 'IAD3' = 'YelMet'))
##' 
##' # Drop relationships lower than minimum cutoff value
##' plot(ans, cutoff = c(0.3,0.6,0.8), style = FALSE)
##' 
##' \dontrun{
##' # Save plot 
##' plot(ans, cutoff = c(0.3, 0.6, 0.8), plot.save = TRUE, filename = 'acrylic-netSEMp1')
##' }

plot.netSEMp1 <- function(x, cutoff = c(0.2, 0.5, 0.8),latent = NULL, plot.save = FALSE, 
                          filename = NULL, style = TRUE, ... ) {
  
  rtp1 <- x$table[,1:5]
  rtp1[, -c(1:3)] <- round(rtp1[, -c(1:3)], 3)
  
  rtp1.a <- rtp1[-c(1:nrow(rtp1)), ]
  rtp1.b <- rtp1[-c(1:nrow(rtp1)), ]
  rtp1.c <- rtp1[-c(1:nrow(rtp1)), ]
  rtp1.d <- rtp1[-c(1:nrow(rtp1)), ]
  
  ### set up cutoffs ############
  count <- length(cutoff)
  ## based on the number of cutoff, to separate different regions
  ## number of cutoff = 1, have two regions
  if (count == 1) {
    c1 <- cutoff[1]  
    rtp1.a <- rtp1[rtp1[, "adj-R-Sqr"] <  c1 , ]
    rtp1.b <- rtp1[rtp1[, "adj-R-Sqr"] >= c1 , ]
  }
  ## number of cutoff = 2, have three regions
  if (count == 2) {
    c1 <- cutoff[1] ; c2 <- cutoff[2] 
    rtp1.a <- rtp1[rtp1[, "adj-R-Sqr"] < c1, ]
    rtp1.b <- rtp1[rtp1[, "adj-R-Sqr"] >= c1 & rtp1[, "adj-R-Sqr"] < c2, ]
    rtp1.c <- rtp1[rtp1[, "adj-R-Sqr"] >= c2 , ]
  }
  ## number of cutoff = 3, have four regions
  if (count == 3) {
    c1 <- cutoff[1] ; c2 <- cutoff[2] ; c3 <- cutoff[3] 
    rtp1.a <- rtp1[rtp1[, "adj-R-Sqr"] < c1, ]
    rtp1.b <- rtp1[rtp1[, "adj-R-Sqr"] >= c1 & rtp1[, "adj-R-Sqr"] < c2, ]
    rtp1.c <- rtp1[rtp1[, "adj-R-Sqr"] >= c2 & rtp1[, "adj-R-Sqr"] < c3, ]
    rtp1.d <- rtp1[rtp1[, "adj-R-Sqr"] >= c3 , ]
  }
  
  name <- colnames(x$data)
  conp1.m1 <- sapply(3:length(name), function(i) { 
    paste0(name[i], "[fillcolor=khaki]") } )
  
  ## assign the latent variable to each mechanic variable
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
    
    ## style is True represents plot the first interval with dotted weak line
    if (style) {
      
      ## set up parameters for the first region
      if (dim(rtp1.a)[1] > 0) {
        conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
          paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@",i,"', 
                 fontsize=80, style=dotted]")
        }
        )
        conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
          paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3], ":", 
                 rtp1.a[i,3], "\\n",
                 colnames(rtp1.a)[5], ": ", rtp1.a[i,5], "')")
        }
        )
      }
      ## set up parameters for the second region
      count.1 <- nrow(rtp1.a)
      if (dim(rtp1.b)[1] > 0) {
        conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
          paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, 
                 "', fontsize=80]")
        }
        )
        conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
          paste0("[", i + count.1, "]", ": ", "paste('", colnames(rtp1.b)[3], 
                 ":", rtp1.b[i,3], "\\n",
                 colnames(rtp1.b)[5], ": ", rtp1.b[i,5], "')")
        }
        )
      }
      ## set up parameters for the third region
      count.2 <- count.1 + nrow(rtp1.b)
      if (dim(rtp1.c)[1] > 0) {
        conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
          paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@", i + count.2, 
                 "', fontsize=80]")
        }
        )
        conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
          paste0("[", i + count.2, "]", ": ", "paste('", colnames(rtp1.c)[3], 
                 ":", rtp1.c[i,3], "\\n", colnames(rtp1.c)[5], ": ", 
                 rtp1.c[i,5], "')")
        }
        )
      }
      ## set up parameters for the forth region
      count.3 <- count.2 + nrow(rtp1.c)
      if (dim(rtp1.d)[1] > 0) {
        conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
          paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, 
                 "', fontsize=80]")
        }
        )
        conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
          paste0("[", i + count.3, "]", ": ", "paste('", colnames(rtp1.d)[3], 
                 ":", rtp1.d[i,3], "\\n", colnames(rtp1.d)[5], ": ", 
                 rtp1.d[i,5], "')")
        }
        )
      }
      
      ## test whether each region has model/parameter or not
      A <- rep(0,4)
      if (exists("conp1.a1") == T) { A[1] <- 1 }
      if (exists("conp1.b1") == T) { A[2] <- 1 }
      if (exists("conp1.c1") == T) { A[3] <- 1 }
      if (exists("conp1.d1") == T) { A[4] <- 1 }
      
      ## all of four regions have models
      if (sum(A) == 4 ) {
        c.plot <- paste0("digraph nice {", "\n",
                         "graph [compound = true, nodesep = 1.8, ranksep =0.9, 
                         layout = dot, rankdir = LR]", "\n", 
                         "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6, width=4, height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.m1,collapse = ";"), "\n", 
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
      ## three regions have models
      if (sum(A) == 3 ) {
        ## only region 'a' does not have models, plot other three regions
        if (exists("conp1.a1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                           paste(conp1.c2,collapse = "\n"), "\n",
                           paste(conp1.d2,collapse = "\n"), "\n"
          )
        }
        ## only region 'b' does not have models, plot other three regions
        if (exists("conp1.b1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname = Helvetica, 
                           shape = diamond, fixedsize = T, fontsize = 100, 
                           color = dodgerblue, style = filled, penwidth = 6, 
                           width = 4, height = 3]", "\n", 
                           paste0(name[1]), "[fillcolor = dodgerblue]", "\n", 
                           "node[fontname=Helvetica,shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n",
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1,collapse = ";"), "\n",
                           "node[fontname=Helvetica, shape=square, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.La1,collapse = ";"), "\n", 
                           "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                           paste0(conp1.L1, collapse = "\n"), "\n",
                           "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                           paste0(conp1.a1,collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=10,arrowsize=2]", "\n", 
                           paste0(conp1.c1,collapse = "\n"), "\n", 
                           "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                           paste0(conp1.d1,collapse = "\n"), "\n", "}", "\n", 
                           paste(conp1.a2,collapse = "\n"), "\n",
                           paste(conp1.c2,collapse = "\n"), "\n",
                           paste(conp1.d2,collapse = "\n"), "\n"
          )
        }
        ## only region 'c' does not have models, plot other three regions
        if (exists("conp1.c1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n",
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                           paste0(conp1.m1,collapse = ";"), "\n",
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
        ## only region 'd' does not have models, plot other three regions
        if (exists("conp1.d1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n",
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## two regions have models
      if (sum(A) == 2 ) {
        ## only region 'a' and 'b' does not have models, plot other two regions
        if (exists("conp1.a1") == F &
            exists("conp1.b1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
        ## only region 'a' and 'c' does not have models, plot other two regions
        if (exists("conp1.a1") == F &
            exists("conp1.c1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n",
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]","\n",
                           paste0(name[2]), "[fillcolor=Violet]", "\n",
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
        ## only region 'a' and 'd' does not have models, plot other two regions
        if (exists("conp1.a1") == F &
            exists("conp1.d1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n",
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n",
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
        ## only region 'b' and 'c' does not have models, plot other two regions
        if (exists("conp1.b1") == F &
            exists("conp1.c1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T,
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
        ## only region 'b' and 'd' does not have models, plot other two regions
        if (exists("conp1.b1") == F &
            exists("conp1.d1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]","\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
        ## only region 'c' and 'd' does not have models, plot other two regions
        if (exists("conp1.c1") == F &
            exists("conp1.d1") == F) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n",
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## only one region has the model
      if (sum(A) == 1 ) {
        ## Region 'a', 'b' and 'c' does not have models, plot the other region
        if (exists("conp1.a1") == F &
            exists("conp1.b1") == F &
            exists("conp1.c1") == F ) {
          c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                           nodesep = 1.8, ranksep =0.9, layout = dot, 
                           rankdir = LR]", "\n", "node[fontname=Helvetica, 
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n",
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n",
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n",
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                           shape=diamond, fixedsize=T, fontsize=100, 
                           color=dodgerblue, style=filled, penwidth=6, width=4, 
                           height=3]", "\n",
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n",
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n",
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
    ## set up parameters for the first region
    if (dim(rtp1.a)[1] > 0) {
      conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
        paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@", i, "', 
               fontsize=80, style=dotted]")
      }
      )
      #conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
      #  paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",
      #  rtp1.a[i,3],"\\n",
      #         colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
      #}
      #)
    }
    ## set up parameters for the second region
    count.1 <- 0
    if (dim(rtp1.b)[1] > 0) {
      conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
        paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, "', 
               fontsize=80]")
      }
      )
      conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
        paste0("[", i + count.1, "]", ": ", "paste('", colnames(rtp1.b)[3], 
               ":", rtp1.b[i,3], "\\n", colnames(rtp1.b)[5], ": ", 
               rtp1.b[i,5], "')")
      }
      )
    }
    ## set up parameters for the third region
    count.2 <- count.1 + nrow(rtp1.b)
    if (dim(rtp1.c)[1] > 0) {
      conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
        paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@", i+count.2, 
               "', fontsize=80]")
      }
      )
      conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
        paste0("[", i + count.2, "]", ": ", "paste('", colnames(rtp1.c)[3], 
               ":", rtp1.c[i,3], "\\n", colnames(rtp1.c)[5], ": ", 
               rtp1.c[i,5], "')")
      }
      )
    }
    ## set up parameters for the forth region
    count.3 <- count.2 + nrow(rtp1.c)
    if (dim(rtp1.d)[1] > 0) {
      conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
        paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, 
               "', fontsize=80]")
      }
      )
      conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
        paste0("[", i + count.3, "]", ": ", "paste('", colnames(rtp1.d)[3], 
               ":", rtp1.d[i,3], "\\n", colnames(rtp1.d)[5], ": ", 
               rtp1.d[i,5], "')")
      }
      )
    }
    
    ## test whether each region has model/parameter or not
    A <- rep(0,4)
    if (exists("conp1.a1") == T) { A[1] <- 1 }
    if (exists("conp1.b1") == T) { A[2] <- 1 }
    if (exists("conp1.c1") == T) { A[3] <- 1 }
    if (exists("conp1.d1") == T) { A[4] <- 1 }
    
    ## all of four regions have models
    if (sum(A) == 4 ) {
      c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                       nodesep = 1.8, ranksep =0.9, layout = dot, 
                       rankdir = LR]", "\n", "node[fontname=Helvetica, 
                       shape=diamond, fixedsize=T, fontsize=100, 
                       color=dodgerblue, style=filled, penwidth=6, width=4, 
                       height=3]", "\n",
                       paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                       "node[fontname=Helvetica, shape=box, fixedsize=T, 
                       fontsize=100, color=dodgerblue, style=filled, 
                       penwidth=6]", "\n",
                       paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
    ## three regions have models
    if (sum(A) == 3 ) {
      ## only region 'a' does not have models, plot other three regions
      if (exists("conp1.a1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n",
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## only region 'b' does not have models, plot other three regions
      if (exists("conp1.b1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n",
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## only region 'c' does not have models, plot other three regions
      if (exists("conp1.c1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep = 0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n",
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## only region 'd' does not have models, plot other three regions
      if (exists("conp1.d1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, 
                         width=4, height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                         paste0(conp1.c1, collapse = "\n"), "\n", "}","\n",
                         paste(conp1.b2, collapse = "\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n" 
        )
      }
    }
    ## two regions have models
    if (sum(A) == 2 ) {
      ## only region 'a' and 'b' does not have models, plot other two regions
      if (exists("conp1.a1") == F &
          exists("conp1.b1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## only region 'a' and 'c' does not have models, plot other two regions
      if (exists("conp1.a1") == F &
          exists("conp1.c1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n",
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## only region 'a' and 'd' does not have models, plot other two regions
      if (exists("conp1.a1") == F &
          exists("conp1.d1") == F ) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
      ## only region 'b' and 'c' does not have models, plot other two regions
      if (exists("conp1.b1") == F &
          exists("conp1.c1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                         # paste(conp1.a2, collapse="\n"), "\n",
                         paste(conp1.d2, collapse = "\n"), "\n" 
        )
      }
      ## only region 'b' and 'd' does not have models, plot other two regions
      if (exists("conp1.b1") == F &
          exists("conp1.d1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, 
                         width=4, height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n", 
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse="\n"), "\n", 
                         "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                         paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                         # paste(conp1.a2, collapse="\n"), "\n", 
                         paste(conp1.c2, collapse = "\n"), "\n" 
        )
      }
      ## only region 'c' and 'd' does not have models, plot other two regions
      if (exists("conp1.c1") == F &
          exists("conp1.d1") == F) {
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n", 
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1, collapse="\n"), "\n", 
                         "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         paste0(conp1.b1, collapse = "\n"), "\n", "}","\n",
                         # paste(conp1.a2, collapse = "\n"), "\n", 
                         paste(conp1.b2, collapse = "\n"), "\n" 
        )
      }
    }
    ## only one region has the model
    if (sum(A) == 1 ) {
      ## Region 'a', 'b' and 'c' does not have models, plot the other region
      if (exists("conp1.a1") == F &
          exists("conp1.b1") == F &
          exists("conp1.c1") == F ) { 
        c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                         nodesep = 1.8, ranksep =0.9, layout = dot, 
                         rankdir = LR]", "\n", "node[fontname=Helvetica, 
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                         shape=diamond, fixedsize=T, fontsize=100, 
                         color=dodgerblue, style=filled, penwidth=6, width=4, 
                         height=3]", "\n", 
                         paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=T, 
                         fontsize=100, color=dodgerblue, style=filled, 
                         penwidth=6]", "\n", 
                         paste0(name[2]), "[fillcolor=Violet]", "\n", 
                         "node[fontname=Helvetica, shape=box, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                         paste0(conp1.m1, collapse = ";"), "\n", 
                         "node[fontname=Helvetica, shape=square, fixedsize=F, 
                         color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                         paste0(conp1.La1,collapse = ";"), "\n", 
                         "edge[color=black, penwidth=8, arrowsize=2]", "\n", 
                         paste0(conp1.L1, collapse = "\n"), "\n",
                         # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                         # paste0(conp1.a1,collapse="\n"),"\n", "}", "\n" 
                         # paste(conp1.a2, collapse="\n"), "\n"
        )
      }
    }
    
  }
    
    } else {
      # latent <- sapply(3:length(name), function(i) {
      #     paste0("M", i-2)
      #   }
      #   )
      # conp1.La1 <- sapply(1:length(latent), function(i) { 
      #   paste0(latent[i], "[fillcolor=LightBlue]") 
      # } 
      # )
      # 
      # conp1.L1 <- sapply(1:length(latent), function(i) {
      #   paste0(name[i+2], "->", latent[i], "[ fontsize=80, 
      #            arrowhead=none]")
      # }
      # )
      
      
      ## style is True represents plot the first interval with dotted weak line
      if (style) {
        ## set up parameters for the first region
        if (dim(rtp1.a)[1] > 0) {
          conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
            paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@",i,"', 
                   fontsize=80, style=dotted]")
          }
          )
          conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
            paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3], ":", 
                   rtp1.a[i,3], "\\n",
                   colnames(rtp1.a)[5], ": ", rtp1.a[i,5], "')")
          }
          )
        }
        ## set up parameters for the second region
        count.1 <- nrow(rtp1.a)
        if (dim(rtp1.b)[1] > 0) {
          conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
            paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, 
                   "', fontsize=80]")
          }
          )
          conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
            paste0("[", i + count.1, "]", ": ", "paste('", colnames(rtp1.b)[3], 
                   ":", rtp1.b[i,3], "\\n",
                   colnames(rtp1.b)[5], ": ", rtp1.b[i,5], "')")
          }
          )
        }
        ## set up parameters for the third region
        count.2 <- count.1 + nrow(rtp1.b)
        if (dim(rtp1.c)[1] > 0) {
          conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
            paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@", i + count.2, 
                   "', fontsize=80]")
          }
          )
          conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
            paste0("[", i + count.2, "]", ": ", "paste('", colnames(rtp1.c)[3], 
                   ":", rtp1.c[i,3], "\\n", colnames(rtp1.c)[5], ": ", 
                   rtp1.c[i,5], "')")
          }
          )
        }
        ## set up parameters for the forth region
        count.3 <- count.2 + nrow(rtp1.c)
        if (dim(rtp1.d)[1] > 0) {
          conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
            paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, 
                   "', fontsize=80]")
          }
          )
          conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
            paste0("[", i + count.3, "]", ": ", "paste('", colnames(rtp1.d)[3], 
                   ":", rtp1.d[i,3], "\\n", colnames(rtp1.d)[5], ": ", 
                   rtp1.d[i,5], "')")
          }
          )
        }
        
        ## test whether each region has model/parameter or not
        A <- rep(0,4)
        if (exists("conp1.a1") == T) { A[1] <- 1 }
        if (exists("conp1.b1") == T) { A[2] <- 1 }
        if (exists("conp1.c1") == T) { A[3] <- 1 }
        if (exists("conp1.d1") == T) { A[4] <- 1 }
        
        ## all of four regions have models
        if (sum(A) == 4 ) {
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9, 
                           layout = dot, rankdir = LR]", "\n", 
                           "node[fontname=Helvetica, shape=diamond, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6, width=4, height=3]", "\n", 
                           paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=T, 
                           fontsize=100, color=dodgerblue, style=filled, 
                           penwidth=6]", "\n", 
                           paste0(name[2]), "[fillcolor=Violet]", "\n", 
                           "node[fontname=Helvetica, shape=box, fixedsize=F, 
                           color=dodgerblue, fontsize=80, penwidth=6]", "\n",
                           paste0(conp1.m1,collapse = ";"), "\n", 
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
        ## three regions have models
        if (sum(A) == 3 ) {
          ## only region 'a' does not have models, plot other three regions
          if (exists("conp1.a1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                             paste(conp1.c2,collapse = "\n"), "\n",
                             paste(conp1.d2,collapse = "\n"), "\n"
            )
          }
          ## only region 'b' does not have models, plot other three regions
          if (exists("conp1.b1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname = Helvetica, 
                             shape = diamond, fixedsize = T, fontsize = 100, 
                             color = dodgerblue, style = filled, penwidth = 6, 
                             width = 4, height = 3]", "\n", 
                             paste0(name[1]), "[fillcolor = dodgerblue]", "\n", 
                             "node[fontname=Helvetica,shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1,collapse = ";"), "\n",
                             "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                             paste0(conp1.a1,collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=10,arrowsize=2]", "\n", 
                             paste0(conp1.c1,collapse = "\n"), "\n", 
                             "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                             paste0(conp1.d1,collapse = "\n"), "\n", "}", "\n", 
                             paste(conp1.a2,collapse = "\n"), "\n",
                             paste(conp1.c2,collapse = "\n"), "\n",
                             paste(conp1.d2,collapse = "\n"), "\n"
            )
          }
          ## only region 'c' does not have models, plot other three regions
          if (exists("conp1.c1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n",
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=F, 
                             color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                             paste0(conp1.m1,collapse = ";"), "\n",
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
          ## only region 'd' does not have models, plot other three regions
          if (exists("conp1.d1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n",
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
        ## two regions have models
        if (sum(A) == 2 ) {
          ## only region 'a' and 'b' does not have models, plot other two regions
          if (exists("conp1.a1") == F &
              exists("conp1.b1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
          ## only region 'a' and 'c' does not have models, plot other two regions
          if (exists("conp1.a1") == F &
              exists("conp1.c1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n",
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]","\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n",
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
          ## only region 'a' and 'd' does not have models, plot other two regions
          if (exists("conp1.a1") == F &
              exists("conp1.d1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n",
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
          ## only region 'b' and 'c' does not have models, plot other two regions
          if (exists("conp1.b1") == F &
              exists("conp1.c1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T,
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
          ## only region 'b' and 'd' does not have models, plot other two regions
          if (exists("conp1.b1") == F &
              exists("conp1.d1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]","\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
          ## only region 'c' and 'd' does not have models, plot other two regions
          if (exists("conp1.c1") == F &
              exists("conp1.d1") == F) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
        ## only one region has the model
        if (sum(A) == 1 ) {
          ## Region 'a', 'b' and 'c' does not have models, plot the other region
          if (exists("conp1.a1") == F &
              exists("conp1.b1") == F &
              exists("conp1.c1") == F ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n", 
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n", 
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n",
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n",
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n",
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
          ## set up parameters for the first region
          if (dim(rtp1.a)[1] > 0) {
            conp1.a1 <- sapply(1:nrow(rtp1.a), function(i) {
              paste0(rtp1.a[i,2], "->", rtp1.a[i,1], "[label='@@", i, "', 
                     fontsize=80, style=dotted]")
            }
            )
            #conp1.a2 <- sapply(1:nrow(rtp1.a), function(i) {
            #  paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",
            #  rtp1.a[i,3],"\\n",
            #         colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
            #}
            #)
          }
          ## set up parameters for the second region
          count.1 <- 0
          if (dim(rtp1.b)[1] > 0) {
            conp1.b1 <- sapply(1:nrow(rtp1.b), function(i) {
              paste0(rtp1.b[i,2], "->", rtp1.b[i,1], "[label='@@", i + count.1, "', 
                     fontsize=80]")
            }
            )
            conp1.b2 <- sapply(1:nrow(rtp1.b), function(i) {
              paste0("[", i + count.1, "]", ": ", "paste('", colnames(rtp1.b)[3], 
                     ":", rtp1.b[i,3], "\\n", colnames(rtp1.b)[5], ": ", 
                     rtp1.b[i,5], "')")
            }
            )
          }
          ## set up parameters for the third region
          count.2 <- count.1 + nrow(rtp1.b)
          if (dim(rtp1.c)[1] > 0) {
            conp1.c1 <- sapply(1:nrow(rtp1.c), function(i) {
              paste0(rtp1.c[i,2], "->", rtp1.c[i,1], "[label='@@", i+count.2, 
                     "', fontsize=80]")
            }
            )
            conp1.c2 <- sapply(1:nrow(rtp1.c), function(i) {
              paste0("[", i + count.2, "]", ": ", "paste('", colnames(rtp1.c)[3], 
                     ":", rtp1.c[i,3], "\\n", colnames(rtp1.c)[5], ": ", 
                     rtp1.c[i,5], "')")
            }
            )
          }
          ## set up parameters for the forth region
          count.3 <- count.2 + nrow(rtp1.c)
          if (dim(rtp1.d)[1] > 0) {
            conp1.d1 <- sapply(1:nrow(rtp1.d), function(i) {
              paste0(rtp1.d[i,2], "->", rtp1.d[i,1], "[label='@@", i + count.3, 
                     "', fontsize=80]")
            }
            )
            conp1.d2 <- sapply(1:nrow(rtp1.d), function(i) {
              paste0("[", i + count.3, "]", ": ", "paste('", colnames(rtp1.d)[3], 
                     ":", rtp1.d[i,3], "\\n", colnames(rtp1.d)[5], ": ", 
                     rtp1.d[i,5], "')")
            }
            )
          }
          
          ## test whether each region has model/parameter or not
          A <- rep(0,4)
          if (exists("conp1.a1") == T) { A[1] <- 1 }
          if (exists("conp1.b1") == T) { A[2] <- 1 }
          if (exists("conp1.c1") == T) { A[3] <- 1 }
          if (exists("conp1.d1") == T) { A[4] <- 1 }
          
          ## all of four regions have models
          if (sum(A) == 4 ) {
            c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                             nodesep = 1.8, ranksep =0.9, layout = dot, 
                             rankdir = LR]", "\n", "node[fontname=Helvetica, 
                             shape=diamond, fixedsize=T, fontsize=100, 
                             color=dodgerblue, style=filled, penwidth=6, width=4, 
                             height=3]", "\n",
                             paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                             "node[fontname=Helvetica, shape=box, fixedsize=T, 
                             fontsize=100, color=dodgerblue, style=filled, 
                             penwidth=6]", "\n",
                             paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
          ## three regions have models
          if (sum(A) == 3 ) {
            ## only region 'a' does not have models, plot other three regions
            if (exists("conp1.a1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n",
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
            ## only region 'b' does not have models, plot other three regions
            if (exists("conp1.b1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n",
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
            ## only region 'c' does not have models, plot other three regions
            if (exists("conp1.c1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep = 0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n",
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
            ## only region 'd' does not have models, plot other three regions
            if (exists("conp1.d1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, 
                               width=4, height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", "}","\n",
                               paste(conp1.b2, collapse = "\n"), "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n" 
              )
            }
          }
          ## two regions have models
          if (sum(A) == 2 ) {
            ## only region 'a' and 'b' does not have models, plot other two regions
            if (exists("conp1.a1") == F &
                exists("conp1.b1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
            ## only region 'a' and 'c' does not have models, plot other two regions
            if (exists("conp1.a1") == F &
                exists("conp1.c1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n",
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
            ## only region 'a' and 'd' does not have models, plot other two regions
            if (exists("conp1.a1") == F &
                exists("conp1.d1") == F ) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
            ## only region 'b' and 'c' does not have models, plot other two regions
            if (exists("conp1.b1") == F &
                exists("conp1.c1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               "edge[color=black, penwidth=14, arrowsize=2]", "\n", 
                               paste0(conp1.d1, collapse = "\n"), "\n", "}", "\n", 
                               # paste(conp1.a2, collapse="\n"), "\n",
                               paste(conp1.d2, collapse = "\n"), "\n" 
              )
            }
            ## only region 'b' and 'd' does not have models, plot other two regions
            if (exists("conp1.b1") == F &
                exists("conp1.d1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, 
                               width=4, height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse="\n"), "\n", 
                               "edge[color=black, penwidth=10, arrowsize=2]", "\n", 
                               paste0(conp1.c1, collapse = "\n"), "\n", "}", "\n", 
                               # paste(conp1.a2, collapse="\n"), "\n", 
                               paste(conp1.c2, collapse = "\n"), "\n" 
              )
            }
            ## only region 'c' and 'd' does not have models, plot other two regions
            if (exists("conp1.c1") == F &
                exists("conp1.d1") == F) {
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1, collapse="\n"), "\n", 
                               "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               paste0(conp1.b1, collapse = "\n"), "\n", "}","\n",
                               # paste(conp1.a2, collapse = "\n"), "\n", 
                               paste(conp1.b2, collapse = "\n"), "\n" 
              )
            }
          }
          ## only one region has the model
          if (sum(A) == 1 ) {
            ## Region 'a', 'b' and 'c' does not have models, plot the other region
            if (exists("conp1.a1") == F &
                exists("conp1.b1") == F &
                exists("conp1.c1") == F ) { 
              c.plot <- paste0("digraph nice {", "\n", "graph [compound = true, 
                               nodesep = 1.8, ranksep =0.9, layout = dot, 
                               rankdir = LR]", "\n", "node[fontname=Helvetica, 
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
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
                               shape=diamond, fixedsize=T, fontsize=100, 
                               color=dodgerblue, style=filled, penwidth=6, width=4, 
                               height=3]", "\n", 
                               paste0(name[1]), "[fillcolor=dodgerblue]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=T, 
                               fontsize=100, color=dodgerblue, style=filled, 
                               penwidth=6]", "\n", 
                               paste0(name[2]), "[fillcolor=Violet]", "\n", 
                               "node[fontname=Helvetica, shape=box, fixedsize=F, 
                               color=dodgerblue, fontsize=80, penwidth=6]", "\n", 
                               paste0(conp1.m1, collapse = ";"), "\n", 
                               # "edge[color=black, penwidth=5, arrowsize=2]", "\n", 
                               # paste0(conp1.a1,collapse="\n"),"\n", "}", "\n" 
                               # paste(conp1.a2, collapse="\n"), "\n"
              )
            }
          }
          
          }
      }
  
  
  
  
  ## This generates syntax to run "grViz" for plotting using above syntax
  ## "LR" is left to right flow
  semplot <- DiagrammeR::grViz(c.plot)
  if (plot.save) {
    ### export graph as follows
    semplot %>% export_svg %>% charToRaw %>% rsvg %>% 
      png::writePNG(paste0(filename,".png"))
  }
  return(semplot)
  }
