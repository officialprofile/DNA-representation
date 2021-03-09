#' @title Plot 2D-Dynamic Graph
#' @description Creates scatter plot
#' @usage plot2DGraph(seq, dim = 2, genbank = FALSE)
#' @return

plot2DGraph <- function(seq, dim = 2, genbank = FALSE, xlab = 'X', ylab = 'Y', main = ''){
  graph <- as.data.frame(dGraph(seq = seq, dim = dim, genbank = genbank)$graph)

  options(scipen = 999)  # Turning-off scientific notation like 1e+48
  ggplot2::theme_set(theme_classic())

  gg <- ggplot(graph, aes(x = X, y = Y)) + geom_point(aes(col = freq)) +
    labs(subtitle = ifelse(genbank, seq, ''), y = ylab, x = xlab, title = 'Dynamic graph',
         col = 'Mass', caption = '')

  plot(gg)
}
