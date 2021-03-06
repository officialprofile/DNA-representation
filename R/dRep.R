#' @title Dynamic Representation of the DNA/RNA Sequences
#' @description Create representation
#' @usage drep(seq, type = 's', dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA type
#' @param type 's' - raw sequence, 'f' - fasta
#' @return A matrix of size 2xN or 3xN, where N is the length of the sequence, giving position of every step of the walk.

dRep <- function(seq, type = 's', dims = 2){
  # Basic error handling
  if (!(dims == 2 || dims == 3)){
    stop(' dim must be equal to 2 or 3')
  }
  if (!(type == 's' || type == 'f')){
    stop(' type must be equal to "s" or "f"')
  }
  if (type == 's' && (typeof(seq) != 'character' || length(seq) > 1)){
    stop(' incorrect sequence')
  }

  # Create vectors
  elements <- setVectors(c(-1,  0,  1), c( 0,  1,  1), c( 1,  0,  1), c( 0, -1,  1))

  data <- dGraph(seq, dim = dims)$graph
  return(print(data))
}
