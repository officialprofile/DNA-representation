#' @title 2D-DYnamic Representation of the DNA/RNA Sequences
#' @param seq sequence given as a single string, e.g. "AGTTGAGGGAG"
#' @param type DNA or RNA

drep <- function(seq, type = 'sequence', dim = 2){

  # Create vectors
  elements <- matrix(,3,5)
  colnames(elements) <- c("A", "C", "G", "T", "U")
  elements[,"A"] <- c(-1,  0,  1)
  elements[,"C"] <- c( 0,  1,  1)
  elements[,"G"] <- c( 1,  0,  1)
  elements[,"T"] <- c( 0, -1,  1)
  elements[,"U"] <- c( 0, -1,  1)

  # Prepare the sequence
  seq <- toupper(seq)
  seq <- str_split(seq, '')
  if (length(unique(seq[[1]])) > 4){
    warning('The number of different characters in the sequence is greater than 4.')
  }

  # Prepare elements for the walk
  N <- length(seq[[1]])
  vec <- rep(0,dim)
  coordinates <- matrix(vec,1,dim)
  colnames(coordinates) <- LETTERS[24:(24+dim-1)]

  # Do the walk
  for (i in 1:N){
    vec <- vec + elements[1:dim, seq[[1]][i]]
    coordinates <- rbind(coordinates, vec)
  }
  rownames(coordinates) <- seq(1:nrow(coordinates))

  return (coordinates[,LETTERS[24:(24+dim-1)]])
}
