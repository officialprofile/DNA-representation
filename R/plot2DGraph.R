#' @title Plot 2D-Dynamic Graph
#' @description Create plot
#' @usage plot2DGraph(seq, type = 's', dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA type
#' @param type 's' - symbolic, 'wb' - grayscale, 'c' - colors
#' @return

plot2DGraph <- function(seq, dim = 2, type = 's', palette = FALSE, xlab = 'X', ylab = 'Y', main = ''){
  graph <- dGraph(seq, dim)$graph
  plot(graph[,'X'], graph[,'Y'], pch = graph[,'freq'], xlab = xlab, ylab = ylab, main = main)
}
