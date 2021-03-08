#' @title Dynamic Representation of the DNA/RNA Sequences
#' @description dRep reates characterization of the dynamic graph that is built by dGraph function. The characterization can be used in an alignment-free analysis of the DNA/RNA sequences.
#' @usage drep(seq, type = 's', dim = 2)
#' @param seq sequence given as a single string, e.g. "AGTGGG" or as a FASTA file
#' @param type 's' - raw sequence (default), 'f' - fasta
#' @param dim number of dimensions: 2 (default) or 3
#' @param vec list of unit vectors associated with nucleobases
#' @return Dataframe with one row and 14 columns (dim = 2) or 26 columns (dim = 3) with descriptors
#' @example dRep('ACGCATGCGGCGAGTG', dim = 2)

dRep <- function(seq, type = 's', dim = 2,
                 vec = list('A' = c(-1, 0, 1), 'C' = c(0, 1, 1), 'G' = c(1, 0, 1), 'T' = c(0, -1, 1))){
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

  # Create Dynamic Graph
  graph <- dGraph(seq, dim = dim, vec = vec)$coordinates

  # Characterize Dynamic Graph
  if (dim == 2){
    data <- matrix(,1,14)
    colnames(data) <- c('len', 'mi_x', 'mi_y', 'sqrt',
                        'I_xx', 'I_yy', 'I_xy', 'I_yx', 'I_11', 'I_22',
                        'Dx1', 'Dx2', 'Dy1', 'Dy2')

    data[1, 'len'] <- length(graph)
    data[1, 'mi_x'] <- sum(graph[, 'X']) / length(graph[, 'X'])
    data[1, 'mi_y'] <- sum(graph[, 'Y']) / length(graph[, 'Y'])
    data[1, 'sqrt'] <- sqrt(data[1, 'mi_x']^2 + data[1, 'mi_y']^2)

    graph_centered <- graph
    graph_centered[, 'X'] <- graph_centered[, 'X'] - data[1, 'mi_x']
    graph_centered[, 'Y'] <- graph_centered[, 'Y'] - data[1, 'mi_y']

    data[1, 'I_xx'] <- sum(graph_centered[, 'Y']^2) # Moments of inertia
    data[1, 'I_yy'] <- sum(graph_centered[, 'X']^2)
    data[1, 'I_xy'] <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])
    data[1, 'I_yx'] <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])

    char_mi <- matrix(,2,2) # We need matrix to calculate eigenvalues
    char_mi[1,1] <- data[1, 'I_xx']
    char_mi[1,2] <- data[1, 'I_xy']
    char_mi[2,1] <- data[1, 'I_yx']
    char_mi[2,2] <- data[1, 'I_yy']
    data[1, 'I_11'] <- eigen(char_mi)$values[1] # Principal moments of inertia
    data[1, 'I_22'] <- eigen(char_mi)$values[2]

    data[1, 'Dx1'] <- data[1, 'mi_x'] / eigen(char_mi)$values[1] # Descriptors
    data[1, 'Dy1'] <- data[1, 'mi_y'] / eigen(char_mi)$values[1]
    data[1, 'Dx2'] <- data[1, 'mi_x'] / eigen(char_mi)$values[2]
    data[1, 'Dy2'] <- data[1, 'mi_y'] / eigen(char_mi)$values[2]
  }

  if (dim == 3){
    data <- matrix(,1,26)
    colnames(data) <- c('len', 'mi_x', 'mi_y', 'mi_z', 'sqrt',
                        'I_xx', 'I_xy', 'I_xz', 'I_yx', 'I_yy', 'I_yz', 'I_zx', 'I_zy', 'I_zz',  'I_11', 'I_22', 'I_33',
                        'Dx1', 'Dx2', 'Dx3', 'Dy1', 'Dy2', 'Dy3', 'Dz1', 'Dz2', 'Dz3')

    data[1, 'len'] <- length(graph)
    data[1, 'mi_x'] <- sum(graph[, 'X']) / length(graph[, 'X'])
    data[1, 'mi_y'] <- sum(graph[, 'Y']) / length(graph[, 'Y'])
    data[1, 'mi_z'] <- sum(graph[, 'Z']) / length(graph[, 'Z'])
    data[1, 'sqrt'] <- sqrt(data[1, 'mi_x']^2 + data[1, 'mi_y']^2 + data[1, 'mi_z']^2)

    graph_centered <- graph
    graph_centered[, 'X'] <- graph_centered[, 'X'] - data[1, 'mi_x']
    graph_centered[, 'Y'] <- graph_centered[, 'Y'] - data[1, 'mi_y']
    graph_centered[, 'Z'] <- graph_centered[, 'Z'] - data[1, 'mi_z']

    data[1, 'I_xx'] <- sum(graph_centered[, 'Y']^2 + graph_centered[, 'Z']^2) # Moments of inertia
    data[1, 'I_yy'] <- sum(graph_centered[, 'X']^2 + graph_centered[, 'Z']^2)
    data[1, 'I_zz'] <- sum(graph_centered[, 'X']^2 + graph_centered[, 'Y']^2)
    data[1, 'I_xy'] <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])
    data[1, 'I_yx'] <- -sum(graph_centered[, 'Y'] * graph_centered[, 'X'])
    data[1, 'I_xz'] <- -sum(graph_centered[, 'X'] * graph_centered[, 'Z'])
    data[1, 'I_zx'] <- -sum(graph_centered[, 'Z'] * graph_centered[, 'X'])
    data[1, 'I_zy'] <- -sum(graph_centered[, 'Z'] * graph_centered[, 'Y'])
    data[1, 'I_yz'] <- -sum(graph_centered[, 'Y'] * graph_centered[, 'Z'])

    char_mi <- matrix(,3,3) # We need matrix to calculate eigenvalues
    char_mi[1, 1] <- data[1, 'I_xx']
    char_mi[1, 2] <- data[1, 'I_xy']
    char_mi[1, 3] <- data[1, 'I_xz']
    char_mi[2, 1] <- data[1, 'I_yx']
    char_mi[2, 2] <- data[1, 'I_yy']
    char_mi[2, 3] <- data[1, 'I_yz']
    char_mi[3, 1] <- data[1, 'I_zx']
    char_mi[3, 2] <- data[1, 'I_zy']
    char_mi[3, 3] <- data[1, 'I_zz']
    data[1, 'I_11'] <- eigen(char_mi)$values[1] # Principal moments of inertia
    data[1, 'I_22'] <- eigen(char_mi)$values[2]
    data[1, 'I_33'] <- eigen(char_mi)$values[3]

    data[1, 'Dx1'] <- data[1, 'mi_x'] / eigen(char_mi)$values[1] # Descriptors
    data[1, 'Dy1'] <- data[1, 'mi_y'] / eigen(char_mi)$values[1]
    data[1, 'Dz1'] <- data[1, 'mi_z'] / eigen(char_mi)$values[1]
    data[1, 'Dx2'] <- data[1, 'mi_x'] / eigen(char_mi)$values[2]
    data[1, 'Dy2'] <- data[1, 'mi_y'] / eigen(char_mi)$values[2]
    data[1, 'Dz2'] <- data[1, 'mi_z'] / eigen(char_mi)$values[2]
    data[1, 'Dx3'] <- data[1, 'mi_x'] / eigen(char_mi)$values[3]
    data[1, 'Dy3'] <- data[1, 'mi_y'] / eigen(char_mi)$values[3]
    data[1, 'Dz3'] <- data[1, 'mi_z'] / eigen(char_mi)$values[3]
  }

  return(as.data.frame(data))
}
