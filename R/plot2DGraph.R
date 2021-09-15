#' @title Plot 2D scatter plot
#' @description Creates scatter plot
#' @usage plot2DGraph(seq, dim = 2, genbank = FALSE)
#' @return

plot2DGraph <- function(seqs, genbank = FALSE, xlab = 'X', ylab = 'Y', main = '',
                        colorset = c('#da186f', '#681f1c', '#ffa600', '#bc5090', '#003f5c'),
                        show.legend = TRUE, legend.pos = 'topleft'){

  graphs <- matrix(,0,4)
  colnames(graphs) <- c('X', 'Y', 'freq', 'nr')
  for (seq in seqs){
    number <- which(seq == seqs)
    graph <- as.data.frame(dGraph(seq = seq, dim = 2, genbank = genbank)$graph)
    graphs <- rbind(graphs, cbind(graph, list('nr' = rep(number, nrow(graph)))))
  }
  graphs_shuffled <- graphs[sample(nrow(graphs)), ] # for the overlaps
  palette <- adjustcolor(colorset[graphs_shuffled$nr], alpha.f = 0.2)
  plot(graphs_shuffled$X, graphs_shuffled$Y, col = palette, pch = 20,
       cex = sqrt(graphs_shuffled$freq), xlab = xlab, ylab = ylab, main = main)
  if (show.legend){
    legend(legend.pos, legend = seqs, col = colorset[1:length(seqs)], pch=16, pt.cex = 2, cex=1, bty = 'n')
  }
}
