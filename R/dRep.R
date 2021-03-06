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

  graph <- dGraph(seq, dim = dims)$coordinates

  # Characterize Dynamic Graph
  if (dims == 2){
    data <- matrix(,1,14)
    colnames(data) <- c("len","mi_x","mi_y","sqrt","I_xx","I_yy","I_xy","I_yx","I_11",
                     "I_22","Dx1","Dx2","Dy1","Dy2")

    data[1, "len"] <- length(graph)
    data[1, "mi_x"] <- sum(graph[,1])/length(graph)
    data[1, "mi_y"] <- sum(graph[,2])/length(graph)
    data[1, "sqrt"] <- sqrt(data[1, "mi_x"]^2 + data[1, "mi_y"]^2)

    graph_centered <- graph
    graph_centered[, 'X'] <- graph_centered[, 'X'] - data[1, "mi_x"]
    graph_centered[, 'Y'] <- graph_centered[, 'Y'] - data[1, "mi_y"]

    data[1, "I_xx"] <- sum(graph_centered[, 'Y']^2) # Moments of interia
    data[1, "I_yy"] <- sum(graph_centered[, 'X']^2)
    data[1, "I_xy"] <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])
    data[1, "I_yx"] <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])

    char_mi <- matrix(,2,2) # We need matrix to calculate eigenvalues
    char_mi[1,1] <- data[1, "I_xx"]
    char_mi[1,2] <- data[1, "I_xy"]
    char_mi[2,1] <- data[1, "I_yx"]
    char_mi[2,2] <- data[1, "I_yy"]
    data[1, "I_11"] <- eigen(char_mi)$values[1] # Principal moments of interia
    data[1, "I_22"] <- eigen(char_mi)$values[2]

    data[1, "Dx1"] <- data[1, "mi_x"] / eigen(char_mi)$values[1] # Descriptors
    data[1, "Dy1"] <- data[1, "mi_y"] / eigen(char_mi)$values[1]
    data[1, "Dx2"] <- data[1, "mi_x"] / eigen(char_mi)$values[2]
    data[1, "Dy2"] <- data[1, "mi_y"] / eigen(char_mi)$values[2]
  }
  return(as.data.frame(data))
}
