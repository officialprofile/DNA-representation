#' @title Dynamic Representation of the DNA/RNA Sequences
#' @description Create representation
#' @usage dwalk(seq, type = 's', dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA type
#' @param type 's' - raw sequence, 'f' - fasta
#' @return A matrix of size 2xN or 3xN, where N is the length of the sequence, giving position of every step of the walk.

dWalk <- function(seq, type = 's', dim = 2){
  # Basic error handling
  if (!(dim == 2 || dim == 3)){
    stop(' dim must be equal to 2 or 3')
  }
  if (!(type == 's' || type == 'f')){
    stop(' type must be equal to "s" or "f"')
  }
  if (type == 's' && (typeof(seq) != 'character' || length(seq) > 1)){
    stop(' incorrect sequence')
  }

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

  return(coordinates[,LETTERS[24:(24 + dim - 1)]])
}
