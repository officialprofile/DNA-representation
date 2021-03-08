#' @title Dynamic Graph
#' @description Create representation
#' @usage dGraph(seq, dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA type
#' @param type 's' - raw sequence, 'f' - fasta
#' @param dim - number of dimensions: 2 (default) or 3
#' @param vec - list of unit vectors associated with nucleobases
#' @return A list containing two matrices: 'coordinates' and 'graph'

dGraph <- function(seq, dim = 2,
                   vec = list('A' = c(-1, 0, 1), 'C' = c(0, 1, 1), 'G' = c(1, 0, 1), 'T' = c(0, -1, 1))){
  # Create vectors
  elements <- setVectors(vec = vec)

  # Prepare the sequence
  seq <- toupper(seq)
  seq <- stringr::str_split(seq, '')
  if (length(unique(seq[[1]])) > 4){
    warning(paste('Number of different characters in the sequence is greater than 4.','\nYour sequence may be incomplete'))
  }

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
  coordinates <- as.data.frame(coordinates)

  dataWithMass <- plyr::count(coordinates)

  return(list('coordinates' = coordinates, 'graph' = dataWithMass))
}
