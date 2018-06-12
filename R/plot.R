##' Plot netSEM result
##' plot.netSEM plots a structural equation network model diagram based on best functional form for each selected pairwise variable.
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
##' @param x An object of class "netSEM", the returned list from \code{netSEMm}. Plotting uses the first element of this list (table) in which the first column of it is endogenous variable, second column is variable and other columns are corresponding best functional form, r-squared, adj-r-squared, P-value1, P-value2 and P-value3.
##' @param res An object of class "subsetData", the returned list from \code{subsetData}. Stronger Solid lines represent relationship with higher adjusted R-sqr and weak dotted lines with less than the first cutoff.
##' @param plot.save True/False, it saves the network diagram plot as png file. The default is false.
##' @param filename A character string naming a file to save as a png file.
##' @param style True/False, it plots the first interval in the network diagram with dotted weak line. The default is True.
##' @param label True/False, it use label to express model and Adj-R2 in the path between variables. The default is True.
##' @param ... A S3 generic/method consistency.
##' 
##' @return An html style plot of pairwise relationship pathway diagram between exogenous variables and an endogenous variable. 
##' Arrows show relationships between each variable with given statistical relations along the connection lines.
##' 
##' @export
##'
##' @seealso \link[netSEM]{netSEMm}
##'
##' @examples
##' # Load acrylic data set
##' data(acrylic)
##' # Build a semi-gSEM model
##' ans <- netSEMm(acrylic)
##' # Subset dataset with three cutoff
##' res <- subsetData(ans,cutoff=c(0.3,0.6,0.8))
##' # Plot the network model 
##' plot(ans,res)
##' # plot the network diagram and save as 'semplot.png' file
##' #plot(ans,res,plot.save=TRUE,filename=c("semplot"))
##' 
##' 

plot.netSEM <- function(x,res,plot.save=FALSE,filename=NULL,style=TRUE,label=TRUE,... ){
  
  rtp1.a <- res[[1]]
  rtp1.b <- res[[2]]
  rtp1.c <- res[[3]]
  rtp1.d <- res[[4]]
  
  name <- colnames(x$data)
  conp1.m1 <- sapply(3:length(name), function(i){ paste0(name[i],"[fillcolor=khaki]") } )
  
  if(label){
    ## style is True represents plot the first interval with dotted weak line
    if (style) {
      
      #conp1.aR2 <- sapply(1:(nrow(rtp1.a) + nrow(rtp1.b) + nrow(rtp1.c) + nrow(rtp1.d)), function(i){ paste0("'@@",i,"'"," [fillcolor=white]") })
      
      if (dim(rtp1.a)[1] > 0) {
        conp1.a1 <- sapply(1:nrow(rtp1.a), function(i){
          paste0(rtp1.a[i,2],"->",rtp1.a[i,1],"[label='@@",i,"',fontsize=80,style=dotted]")
        }
        )
        conp1.a2 <- sapply(1:nrow(rtp1.a), function(i){
          paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",rtp1.a[i,3],"\\n",
                 colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
        }
        )
      }
      
      count.1 <- nrow(rtp1.a)
      if (dim(rtp1.b)[1] > 0) {
        conp1.b1 <- sapply(1:nrow(rtp1.b), function(i){
          paste0(rtp1.b[i,2],"->",rtp1.b[i,1],"[label='@@",i + count.1,"',fontsize=80]")
        }
        )
        conp1.b2 <- sapply(1:nrow(rtp1.b), function(i){
          paste0("[",i + count.1,"]",": ","paste('", colnames(rtp1.b)[3],":",rtp1.b[i,3],"\\n",
                 colnames(rtp1.b)[5],": ",rtp1.b[i,5],"')")
        }
        )
      }
      
      count.2 <- count.1 + nrow(rtp1.b)
      if (dim(rtp1.c)[1] > 0) {
        conp1.c1 <- sapply(1:nrow(rtp1.c), function(i){
          paste0(rtp1.c[i,2],"->",rtp1.c[i,1],"[label='@@",i + count.2,"',fontsize=80]")
        }
        )
        conp1.c2 <- sapply(1:nrow(rtp1.c), function(i){
          paste0("[",i + count.2,"]",": ","paste('", colnames(rtp1.c)[3],":",rtp1.c[i,3],"\\n",
                 colnames(rtp1.c)[5],": ",rtp1.c[i,5],"')")
        }
        )
      }
      
      count.3 <- count.2 + nrow(rtp1.c)
      if (dim(rtp1.d)[1] > 0) {
        conp1.d1 <- sapply(1:nrow(rtp1.d), function(i){
          paste0(rtp1.d[i,2],"->",rtp1.d[i,1],"[label='@@",i + count.3,"',fontsize=80]")
        }
        )
        conp1.d2 <- sapply(1:nrow(rtp1.d), function(i){
          paste0("[",i + count.3,"]",": ","paste('", colnames(rtp1.d)[3],":",rtp1.d[i,3],"\\n",
                 colnames(rtp1.d)[5],": ",rtp1.d[i,5],"')")
        }
        )
      }

      A <- rep(0,4)
      if (exists("conp1.a1") == T) { A[1] <- 1 }
      if (exists("conp1.b1") == T) { A[2] <- 1 }
      if (exists("conp1.c1") == T) { A[3] <- 1 }
      if (exists("conp1.d1") == T) { A[4] <- 1 }
      
      
      if (sum(A) == 4 ) {
        
        c.plot <- paste0("digraph nice {", "\n",
                         "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                         "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                         paste0(name[1]),"[fillcolor=lightcoral]","\n",
                         "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                         paste0(name[2]),"[fillcolor=Palegreen]","\n",
                         "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                         paste0(conp1.m1,collapse = ";"),"\n",
                         #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                         #paste0(conp1.aR2,collapse=";"),"\n",
                         "edge[color=black,penwidth=5,arrowsize=2]","\n",
                         paste0(conp1.a1,collapse = "\n"),"\n",
                         "edge[color=black,penwidth=5,arrowsize=2]","\n",
                         paste0(conp1.b1,collapse = "\n"),"\n",
                         "edge[color=black,penwidth=10,arrowsize=2]","\n",
                         paste0(conp1.c1,collapse = "\n"),"\n",
                         "edge[color=black,penwidth=14,arrowsize=2]","\n",
                         paste0(conp1.d1,collapse = "\n"),"\n",
                         "}","\n",
                         paste(conp1.a2,collapse = "\n") ,"\n",
                         paste(conp1.b2,collapse = "\n") ,"\n",
                         paste(conp1.c2,collapse = "\n") ,"\n",
                         paste(conp1.d2,collapse = "\n") ,"\n"
        )
        
      }
      
      if (sum(A) == 3 ) {
        if (exists("conp1.a1") == F) {
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse = ";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 2 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 1 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n"
          )
        }
      }
      
      
    }else{
      
      #conp1.aR2 <- sapply(1:(nrow(rtp1.b) + nrow(rtp1.c) + nrow(rtp1.d)), function(i){ paste0("'@@",i,"'"," [fillcolor=white]") })
      
      if (dim(rtp1.a)[1] > 0) {
        conp1.a1 <- sapply(1:nrow(rtp1.a), function(i){
          paste0(rtp1.a[i,2],"->",rtp1.a[i,1],"[label='@@",i,"',fontsize=80,style=dotted]")
        }
        )
        #conp1.a2 <- sapply(1:nrow(rtp1.a), function(i){
        #  paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",rtp1.a[i,3],"\\n",
        #         colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
        #}
        #)
      }
      
      count.1 <- 0
      if (dim(rtp1.b)[1] > 0) {
        conp1.b1 <- sapply(1:nrow(rtp1.b), function(i){
          paste0(rtp1.b[i,2],"->",rtp1.b[i,1],"[label='@@",i+count.1,"',fontsize=80]")
        }
        )
        conp1.b2 <- sapply(1:nrow(rtp1.b), function(i){
          paste0("[",i + count.1,"]",": ","paste('", colnames(rtp1.b)[3],":",rtp1.b[i,3],"\\n",
                 colnames(rtp1.b)[5],": ",rtp1.b[i,5],"')")
        }
        )
      }
      
      count.2 <- count.1 + nrow(rtp1.b)
      if (dim(rtp1.c)[1] > 0) {
        conp1.c1 <- sapply(1:nrow(rtp1.c), function(i){
          paste0(rtp1.c[i,2],"->",rtp1.c[i,1],"[label='@@",i+count.2,"',fontsize=80]")
        }
        )
        conp1.c2 <- sapply(1:nrow(rtp1.c), function(i){
          paste0("[",i + count.2,"]",": ","paste('", colnames(rtp1.c)[3],":",rtp1.c[i,3],"\\n",
                 colnames(rtp1.c)[5],": ",rtp1.c[i,5],"')")
        }
        )
      }
      
      count.3 <- count.2 + nrow(rtp1.c)
      if (dim(rtp1.d)[1] > 0) {
        conp1.d1 <- sapply(1:nrow(rtp1.d), function(i){
          paste0(rtp1.d[i,2],"->",rtp1.d[i,1],"[label='@@",i+count.3,"',fontsize=80]")
        }
        )
        conp1.d2 <- sapply(1:nrow(rtp1.d), function(i){
          paste0("[",i + count.3,"]",": ","paste('", colnames(rtp1.d)[3],":",rtp1.d[i,3],"\\n",
                 colnames(rtp1.d)[5],": ",rtp1.d[i,5],"')")
        }
        )
      }
      
      A <- rep(0,4)
      if (exists("conp1.a1") == T) { A[1] <- 1 }
      if (exists("conp1.b1") == T) { A[2] <- 1 }
      if (exists("conp1.c1") == T) { A[3] <- 1 }
      if (exists("conp1.d1") == T) { A[4] <- 1 }
      
      
      if (sum(A) == 4 ) {
        
        c.plot <- paste0("digraph nice {", "\n",
                         "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                         "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                         paste0(name[1]),"[fillcolor=lightcoral]","\n",
                         "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                         paste0(name[2]),"[fillcolor=Palegreen]","\n",
                         "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                         paste0(conp1.m1,collapse=";"),"\n",
                         #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                         #paste0(conp1.aR2,collapse=";"),"\n",
                         #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                         #paste0(conp1.a1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=5,arrowsize=2]","\n",
                         paste0(conp1.b1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=10,arrowsize=2]","\n",
                         paste0(conp1.c1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=14,arrowsize=2]","\n",
                         paste0(conp1.d1,collapse="\n"),"\n",
                         "}","\n",
                         #paste(conp1.a2,collapse="\n") ,"\n",
                         paste(conp1.b2,collapse="\n") ,"\n",
                         paste(conp1.c2,collapse="\n") ,"\n",
                         paste(conp1.d2,collapse="\n") ,"\n"
        )
        
      }
      
      if (sum(A) == 3 ) {
        if(exists("conp1.a1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 2 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 1 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 1.8, ranksep =0.9,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           #"node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           #paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "}","\n"
                           #paste(conp1.a2,collapse="\n") ,"\n"
          )
        }
      }
      
    }
  }else{
    ## style is True represents plot the first interval with dotted weak line
    if (style) {
      
      conp1.aR2 <- sapply(1:(nrow(rtp1.a) + nrow(rtp1.b) + nrow(rtp1.c) + nrow(rtp1.d)), function(i){ paste0("'@@",i,"'"," [fillcolor=white]") })
      
      if (dim(rtp1.a)[1] > 0) {
        conp1.a1 <- sapply(1:nrow(rtp1.a), function(i){
          paste0(rtp1.a[i,2],"->'@@",i,"'","->",rtp1.a[i,1]," [style=dotted]")
        }
        )
        conp1.a2 <- sapply(1:nrow(rtp1.a), function(i){
          paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",rtp1.a[i,3],"\\n",
                 colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
        }
        )
      }
      
      count.1 <- nrow(rtp1.a)
      if (dim(rtp1.b)[1] > 0) {
        conp1.b1 <- sapply(1:nrow(rtp1.b), function(i){
          paste0(rtp1.b[i,2],"->'@@",i + count.1,"'","->",rtp1.b[i,1])
        }
        )
        conp1.b2 <- sapply(1:nrow(rtp1.b), function(i){
          paste0("[",i + count.1,"]",": ","paste('", colnames(rtp1.b)[3],":",rtp1.b[i,3],"\\n",
                 colnames(rtp1.b)[5],": ",rtp1.b[i,5],"')")
        }
        )
      }
      
      count.2 <- count.1 + nrow(rtp1.b)
      if (dim(rtp1.c)[1] > 0) {
        conp1.c1 <- sapply(1:nrow(rtp1.c), function(i){
          paste0(rtp1.c[i,2],"->'@@",i + count.2,"'","->",rtp1.c[i,1])
        }
        )
        conp1.c2 <- sapply(1:nrow(rtp1.c), function(i){
          paste0("[",i + count.2,"]",": ","paste('", colnames(rtp1.c)[3],":",rtp1.c[i,3],"\\n",
                 colnames(rtp1.c)[5],": ",rtp1.c[i,5],"')")
        }
        )
      }
      
      count.3 <- count.2 + nrow(rtp1.c)
      if (dim(rtp1.d)[1] > 0) {
        conp1.d1 <- sapply(1:nrow(rtp1.d), function(i){
          paste0(rtp1.d[i,2],"->'@@",i + count.3,"'","->",rtp1.d[i,1])
        }
        )
        conp1.d2 <- sapply(1:nrow(rtp1.d), function(i){
          paste0("[",i + count.3,"]",": ","paste('", colnames(rtp1.d)[3],":",rtp1.d[i,3],"\\n",
                 colnames(rtp1.d)[5],": ",rtp1.d[i,5],"')")
        }
        )
      }
      
      
      ## This generates syntax to run "grViz" for plotting using above syntax
      ## "LR" is left to right flow
      
      A <- rep(0,4)
      if (exists("conp1.a1") == T) { A[1] <- 1 }
      if (exists("conp1.b1") == T) { A[2] <- 1 }
      if (exists("conp1.c1") == T) { A[3] <- 1 }
      if (exists("conp1.d1") == T) { A[4] <- 1 }
      
      
      if (sum(A) == 4 ) {
        
        c.plot <- paste0("digraph nice {", "\n",
                         "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                         "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                         paste0(name[1]),"[fillcolor=lightcoral]","\n",
                         "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                         paste0(name[2]),"[fillcolor=Palegreen]","\n",
                         "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                         paste0(conp1.m1,collapse=";"),"\n",
                         "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                         paste0(conp1.aR2,collapse=";"),"\n",
                         "edge[color=black,penwidth=5,arrowsize=2]","\n",
                         paste0(conp1.a1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=5,arrowsize=2]","\n",
                         paste0(conp1.b1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=10,arrowsize=2]","\n",
                         paste0(conp1.c1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=14,arrowsize=2]","\n",
                         paste0(conp1.d1,collapse="\n"),"\n",
                         "}","\n",
                         paste(conp1.a2,collapse="\n") ,"\n",
                         paste(conp1.b2,collapse="\n") ,"\n",
                         paste(conp1.c2,collapse="\n") ,"\n",
                         paste(conp1.d2,collapse="\n") ,"\n"
        )
        
      }
      
      if (sum(A) == 3 ) {
        if(exists("conp1.a1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 2 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 1 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.a1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.a2,collapse="\n") ,"\n"
          )
        }
      }
      
      
    }else{
      
      conp1.aR2 <- sapply(1:(nrow(rtp1.b) + nrow(rtp1.c) + nrow(rtp1.d)), function(i){ paste0("'@@",i,"'"," [fillcolor=white]") })
      
      if (dim(rtp1.a)[1] > 0) {
        conp1.a1 <- sapply(1:nrow(rtp1.a), function(i){
          paste0(rtp1.a[i,2],"->'@@",i,"'","->",rtp1.a[i,1]," [style=dotted]")
        }
        )
        #conp1.a2 <- sapply(1:nrow(rtp1.a), function(i){
        #  paste0("[",i,"]",": ","paste('", colnames(rtp1.a)[3],":",rtp1.a[i,3],"\\n",
        #         colnames(rtp1.a)[5],": ",rtp1.a[i,5],"')")
        #}
        #)
      }
      
      count.1 <- 0
      if (dim(rtp1.b)[1] > 0) {
        conp1.b1 <- sapply(1:nrow(rtp1.b), function(i){
          paste0(rtp1.b[i,2],"->'@@",i + count.1,"'","->",rtp1.b[i,1])
        }
        )
        conp1.b2 <- sapply(1:nrow(rtp1.b), function(i){
          paste0("[",i + count.1,"]",": ","paste('", colnames(rtp1.b)[3],":",rtp1.b[i,3],"\\n",
                 colnames(rtp1.b)[5],": ",rtp1.b[i,5],"')")
        }
        )
      }
      
      count.2 <- count.1 + nrow(rtp1.b)
      if (dim(rtp1.c)[1] > 0) {
        conp1.c1 <- sapply(1:nrow(rtp1.c), function(i){
          paste0(rtp1.c[i,2],"->'@@",i + count.2,"'","->",rtp1.c[i,1])
        }
        )
        conp1.c2 <- sapply(1:nrow(rtp1.c), function(i){
          paste0("[",i + count.2,"]",": ","paste('", colnames(rtp1.c)[3],":",rtp1.c[i,3],"\\n",
                 colnames(rtp1.c)[5],": ",rtp1.c[i,5],"')")
        }
        )
      }
      
      count.3 <- count.2 + nrow(rtp1.c)
      if (dim(rtp1.d)[1] > 0) {
        conp1.d1 <- sapply(1:nrow(rtp1.d), function(i){
          paste0(rtp1.d[i,2],"->'@@",i + count.3,"'","->",rtp1.d[i,1])
        }
        )
        conp1.d2 <- sapply(1:nrow(rtp1.d), function(i){
          paste0("[",i + count.3,"]",": ","paste('", colnames(rtp1.d)[3],":",rtp1.d[i,3],"\\n",
                 colnames(rtp1.d)[5],": ",rtp1.d[i,5],"')")
        }
        )
      }
      
      
      ## This generates syntax to run "grViz" for plotting using above syntax
      ## "LR" is left to right flow
      
      A <- rep(0,4)
      if (exists("conp1.a1") == T) { A[1] <- 1 }
      if (exists("conp1.b1") == T) { A[2] <- 1 }
      if (exists("conp1.c1") == T) { A[3] <- 1 }
      if (exists("conp1.d1") == T) { A[4] <- 1 }
      
      
      if (sum(A) == 4 ) {
        
        c.plot <- paste0("digraph nice {", "\n",
                         "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                         "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                         paste0(name[1]),"[fillcolor=lightcoral]","\n",
                         "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                         paste0(name[2]),"[fillcolor=Palegreen]","\n",
                         "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                         paste0(conp1.m1,collapse=";"),"\n",
                         "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                         paste0(conp1.aR2,collapse=";"),"\n",
                         #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                         #paste0(conp1.a1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=5,arrowsize=2]","\n",
                         paste0(conp1.b1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=10,arrowsize=2]","\n",
                         paste0(conp1.c1,collapse="\n"),"\n",
                         "edge[color=black,penwidth=14,arrowsize=2]","\n",
                         paste0(conp1.d1,collapse="\n"),"\n",
                         "}","\n",
                         #paste(conp1.a2,collapse="\n") ,"\n",
                         paste(conp1.b2,collapse="\n") ,"\n",
                         paste(conp1.c2,collapse="\n") ,"\n",
                         paste(conp1.d2,collapse="\n") ,"\n"
        )
        
      }
      
      if (sum(A) == 3 ) {
        if(exists("conp1.a1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 2 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           #paste(conp1.a2,collapse="\n") ,"\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
      }
      
      if (sum(A) == 1 ) {
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.c1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=14,arrowsize=2]","\n",
                           paste0(conp1.d1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.d2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.b1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=10,arrowsize=2]","\n",
                           paste0(conp1.c1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.c2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.a1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F ){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           "edge[color=black,penwidth=5,arrowsize=2]","\n",
                           paste0(conp1.b1,collapse="\n"),"\n",
                           "}","\n",
                           paste(conp1.b2,collapse="\n") ,"\n"
          )
        }
        if(exists("conp1.b1") == F & exists("conp1.c1") == F & exists("conp1.d1") == F){
          c.plot <- paste0("digraph nice {", "\n",
                           "graph [compound = true, nodesep = 2, ranksep =0.25,layout = dot,rankdir = LR]","\n",
                           "node[fontname=Helvetica,shape=diamond,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6,width=4,height=3]","\n",
                           paste0(name[1]),"[fillcolor=lightcoral]","\n",
                           "node[fontname=Helvetica,shape=oval,fixedsize=T,fontsize=100,color=dodgerblue,style=filled,penwidth=6]","\n",
                           paste0(name[2]),"[fillcolor=Palegreen]","\n",
                           "node[fontname=Helvetica,shape=box,fixedsize=F,color=dodgerblue,fontsize=80,penwidth=6]","\n",
                           paste0(conp1.m1,collapse=";"),"\n",
                           "node[fontname=Helvetica,shape=egg,fixedsize=T,color=dodgerblue,fontsize=40,penwidth=5]","\n",
                           paste0(conp1.aR2,collapse=";"),"\n",
                           #"edge[color=black,penwidth=5,arrowsize=2]","\n",
                           #paste0(conp1.a1,collapse="\n"),"\n",
                           "}","\n"
                           #paste(conp1.a2,collapse="\n") ,"\n"
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
    semplot %>% export_svg %>% charToRaw %>% rsvg %>% png::writePNG(paste0(filename,".png"))
    
  }
  
  return(semplot)
}
