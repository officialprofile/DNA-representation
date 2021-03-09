#' @title Dynamic Graph
#' @description Create representation
#' @usage dGraph(seq, dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA type
#' @param genbank FALSE (default) - raw sequence, TRUE - accession number from GenBank
#' @param dim - number of dimensions: 2 (default) or 3
#' @param vec - list of unit vectors associated with nucleobases
#' @return A list containing two matrices: 'coordinates' and 'graph'

dGraph <- function(seq, dim = 2, genbank = FALSE,
                   vec = list('A' = c(-1, 0, 1), 'C' = c(0, 1, 1), 'G' = c(1, 0, 1), 'T' = c(0, -1, 1))){
  # Create vectors
  elements <- setVectors(vec = vec)

  # Prepare the sequence
  if (genbank == FALSE) seq <- stringr::str_split(seq, '')
  if (genbank == TRUE) seq <- ape::read.GenBank(seq, species.names = TRUE, as.character = TRUE)
  seq <- toupper(seq[[1]])

  # Prepare elements for the walk
  N <- length(seq)
  pos <- rep(0, dim) # Position
  coordinates <- matrix(, 0, dim)
  colnames(coordinates) <- LETTERS[24:(24+dim-1)]

  # Do the walk
  for (i in 1:N){
    if (length(intersect(seq[i], c('A', 'C', 'G', 'T', 'U')))){
      pos <- pos + elements[1:dim, seq[i]]
      coordinates <- rbind(coordinates, pos)
    }
  }
  rownames(coordinates) <- seq(1:nrow(coordinates))
  coordinates <- as.data.frame(coordinates)

  dataWithMass <- plyr::count(coordinates)

  return(list('coordinates' = coordinates, 'graph' = dataWithMass))
}
