#' @title Dynamic Graph
#' @description Create representation
#' @usage dGraph(seq, dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA type
#' @param type 's' - raw sequence, 'f' - fasta
#' @return A matrix of size Nx(dim+1) where N is the number of unique positions on the graph
dGraph <- function(seq, dim = 2){
  # Create vectors
  elements <- setVectors(c(-1,  0,  1), c( 0,  1,  1), c( 1,  0,  1), c( 0, -1,  1))

  # Prepare elements for the walk
  N <- length(seq[[1]])
  vec <- rep(0, dim)
  coordinates <- matrix(, 0, dim)
  colnames(coordinates) <- LETTERS[24:(24+dim-1)]

  # Do the walk
  for (i in 1:N){
    vec <- vec + elements[1:dim, seq[[1]][i]]
    coordinates <- rbind(coordinates, vec)
  }
  rownames(coordinates) <- seq(1:nrow(coordinates))

  x <- coordinates[,LETTERS[24:(24 + dim - 1)]]
  mass <- rep(0, nrow(x))
  for (i in 1:nrow(x)){
    for (j in 1:nrow(x)){
      if (sum(x[i,] == x[j,]) == length(x[1,])){
        mass[i] = mass[i] + 1
      }
    }
  }
  x <- cbind(x[,1:dim], mass)

  return(unique(x))
}
