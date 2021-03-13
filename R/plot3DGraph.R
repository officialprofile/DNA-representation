#' @title Plot 3D scatter plot
#' @description Creates 3D scatter plot
#' @usage plot3DGraph(seq, genbank = FALSE)

plot3DGraph <- function(seqs, genbank = FALSE, xlab = 'X', ylab = 'Y', zlab = 'Z', main = '3D Graph',
                        colorset = c('#da186f', '#681f1c', '#ffa600', '#bc5090', '#003f5c'),
                        radius = 5){
  graphs <- matrix(,0,5)
  colnames(graphs) <- c('X', 'Y', 'Z', 'freq', 'nr')
  for (seq in seqs){
    number <- which(seq == seqs)
    graph <- as.data.frame(dGraph(seq = seq, dim = 3, genbank = genbank)$graph)
    graphs <- rbind(graphs, cbind(graph, list('nr' = rep(number, nrow(graph)))))
  }
  graphs_shuffled <- graphs[sample(nrow(graphs)), ] # for the overlaps

  plot3d(graphs_shuffled$X, graphs_shuffled$Y, graphs_shuffled$Z,
         col = colorset[graphs_shuffled$nr], type = "s",
         radius = radius, xlab = xlab, ylab = ylab, zlab = zlab)
}



