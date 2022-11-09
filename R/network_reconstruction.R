#' @title convert ODE results(ODE_solving2) to basic network plot table
#' @param result list result from qsODE_parallel
#' @return a list with basic information to plot network
#' @export
network_conversion <- function(result){
  n = ncol(result$fit)

  effect.mean = apply(result$fit,2,mean)[4:n]
  effect.predict.mean = apply(result$predict,2,mean)[4:n]

  temp = matrix(NA,nrow = n-3, ncol=3)
  colnames(temp) = c("From", "To", "Effect")
  temp[,1] = colnames(result$fit)[4:n]
  temp[,2] = colnames(result$fit)[4]
  temp[,3] = effect.predict.mean
  if (nrow(temp)==2) {
    temp = t(data.frame(temp[-1,]))
  } else{
    temp = data.frame(temp[-1,])
  }
  output <- list(ind.name = colnames(result$fit)[4],
                 dep.name = colnames(result$fit)[5:n],
                 ODE.par = result$ODE.value,
                 ind.par = result$LOP_par[,3],
                 dep.par = result$LOP_par[,4:(n-1)],
                 effect.mean = effect.predict.mean,
                 effect.all = result$fit,
                 edge = temp,
                 ind.effect = effect.predict.mean[1])
  return(output)
}


#' @title generate network plot
#' @import igraph
#' @importFrom stats aggregate
#' @param result list result from network_conversion
#' @param title text for plot title
#' @return network plot
#' @export
network_plot <- function(result, title = NULL){
  #ind effect control node size
  extra <- sapply(result,"[[", "ind.effect")

  #edge
  after <- data.frame(do.call(rbind, lapply(result, "[[", "edge")))
  after$Effect = as.numeric(after$Effect)
  rownames(after) = NULL


  after$edge.colour = NA
  for (i in 1:nrow(after)) {
    if(after$Effect[i]>=0){
      after$edge.colour[i] = "#FE433C"
    } else{
      after$edge.colour[i] = "#0095EF"
    }
  }

  #nodes
  nodes <- data.frame(unique(after[,2]),unique(after[,2]),extra)
  colnames(nodes) <- c("id","name","ind_effect")
  nodes$influence <- aggregate(Effect ~ To, data = after, sum)[,2]
  nodes$node.colour = NA
  for (i in 1:nrow(nodes)) {
    if(nodes$influence[i]>=0){
      nodes$node.colour[i] = "#FFC4C4"
    } else{
      nodes$node.colour[i] = "#89CFFD"
    }
  }

  #normalization
  normalization <- function(x){(x-min(x))/(max(x)-min(x))+0.2}

  #final plot
  after[,3] <- normalization(abs(after[,3]))
  nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  #after[,3] <- abs(after[,3])
  #nodes[,3:4] <- abs(nodes[,3:4])

  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )

  #layout
  l <- layout_randomly(net)

  plot.igraph(net,
              vertex.label=V(net)$name,
              vertex.label.color="black",
              vertex.shape="circle",
              vertex.label.cex=V(net)$ind_effect*1.5,
              vertex.size=V(net)$ind_effect*20+5,
              edge.curved=0.05,
              edge.color=E(net)$edge.colour,
              edge.frame.color=E(net)$edge.colour,
              edge.width=E(net)$Effect*5,
              vertex.color=V(net)$node.colour,
              layout=l,
              main=title,
              margin=c(-.05,-.05,-.05,-.05)
  )
}

