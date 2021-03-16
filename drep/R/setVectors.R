#' @description Set vectors
#' @usage setVectors(vec)
#' @param vec list containing coordinates of the unit vectors
#' @return A matrix of size 3x5 (regardless of dimensionality)

setVectors <- function(vec){
  # Error handling
  lengths <- c()
  for (v in names(vec)){
    lengths <- c(lengths, length(vec[v][[1]]))
  }
  if (min(lengths) != max(lengths)){
    stop('Number of dimensions must be the same for every nucleobase.')
  }

  # Add third dimension for two dimensional case
  if (lengths[1] == 2){
    for (v in names(vec)){
      vec[v][[1]] <- c(vec[v][[1]], 1)
    }
  }

  elements <- matrix(,3,5)
  colnames(elements) <- c('A', 'C', 'G', 'T', 'U')
  rownames(elements) <- c('X', 'Y', 'Z')
  elements[, 'A'] <- vec['A'][[1]]
  elements[, 'C'] <- vec['C'][[1]]
  elements[, 'G'] <- vec['G'][[1]]
  elements[, 'T'] <- vec['T'][[1]]
  elements[, 'U'] <- vec['T'][[1]]

  if (lengths[1] == 2){
    return (elements[1:2, ])
  }
  else {
    return (elements)
  }
}
