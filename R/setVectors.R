#' @description Create representation
#' @usage dwalk(seq, type = 's', dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA type
#' @param type 's' - raw sequence, 'f' - fasta
#' @return A matrix of size 2xN or 3xN, where N is the length of the sequence, giving position of every step of the walk.

setVectors <- function(A, C, G, TU){
  lengths <- c()
  for (vec in c(A, C, G, TU)){
    lengths <- c(lengths, length(vec))
  }
  if (lengths[1] == 2){
    for (vec in c(A, C, G, TU)){
      A <- c(A, 1)
    }
  }

  elements <- matrix(,3,5)
  colnames(elements) <- c("A", "C", "G", "T", "U")
  elements[,"A"] <- A
  elements[,"C"] <- C
  elements[,"G"] <- G
  elements[,"T"] <- TU
  elements[,"U"] <- TU

  if (lengths[1] == 2){
    return (elements[1:2, ])
  }
  else {
    return (elements)
  }

}
