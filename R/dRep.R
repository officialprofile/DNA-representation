#' @title Numerical Representation of DNA Sequences
#' @usage dRep(seqs, dim = 2, genbank = FALSE, vec = list('A' = c(-1, 0, 1),
#' 'C' = c(0, 1, 1), 'G' = c(1, 0, 1), 'T' = c(0, -1, 1)))
#' @param seqs vector of sequences given as strings or GenBank accession numbers
#' @param genbank FALSE (default) - sequence, TRUE - accession number from GenBank
#' @param dim number of dimensions: 2 (default) or 3
#' @param vec list of unit vectors associated with nucleobases
#' @return Dataframe with N row and 14 columns (dim = 2) or 26 columns (dim = 3), where each row represents one of N given sequences
#' @example dRep('ACGCATGCGGCGAGCCAATG', dim = 3)
#' @example dRep(c('KX369547', 'HQ234498'), dim = 2, genbank = TRUE)

dRep <- function(seqs, dim = 2, genbank = FALSE,
                 vec = list('A' = c(-1, 0, 1), 'C' = c(0, 1, 1), 'G' = c(1, 0, 1), 'T' = c(0, -1, 1))){
  # Prepare dataset
  if (dim == 2){
    dataset <- matrix(,0,14)
    colnames(dataset) <- c('len', 'mi_x', 'mi_y', 'sqrt',
                        'I_xx', 'I_yy', 'I_xy', 'I_yx', 'I_11', 'I_22',
                        'Dx1', 'Dx2', 'Dy1', 'Dy2')
  }
  if (dim == 3){
    dataset <- matrix(,0,26)
    colnames(dataset) <- c('len', 'mi_x', 'mi_y', 'mi_z', 'sqrt',
                        'I_xx', 'I_xy', 'I_xz', 'I_yx', 'I_yy', 'I_yz', 'I_zx', 'I_zy', 'I_zz',  'I_11', 'I_22', 'I_33',
                        'Dx1', 'Dx2', 'Dx3', 'Dy1', 'Dy2', 'Dy3', 'Dz1', 'Dz2', 'Dz3')
  }

  for (seq in seqs){
    graph <- dGraph(seq, genbank = genbank, dim = dim, vec = vec)$coordinates
    data <- list()

    if (dim == 2){
      data$len <- nrow(graph)
      data$mi_x <- sum(graph[, 'X']) / length(graph[, 'X'])
      data$mi_y <- sum(graph[, 'Y']) / length(graph[, 'Y'])
      data$sqrt <- sqrt(data$mi_x^2 + data$mi_y^2)

      graph_centered <- graph
      graph_centered[, 'X'] <- graph_centered[, 'X'] - data$mi_x
      graph_centered[, 'Y'] <- graph_centered[, 'Y'] - data$mi_y

      data$I_xx <- sum(graph_centered[, 'Y']^2) # Moments of inertia
      data$I_yy <- sum(graph_centered[, 'X']^2)
      data$I_xy <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])
      data$I_yx <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])

      char_mi <- matrix(,2,2) # We need matrix to calculate eigenvalues
      char_mi[1,1] <- data$I_xx
      char_mi[1,2] <- data$I_xy
      char_mi[2,1] <- data$I_yx
      char_mi[2,2] <- data$I_yy
      data$I_11 <- eigen(char_mi)$values[1] # Principal moments of inertia
      data$I_22 <- eigen(char_mi)$values[2]

      data$Dx1 <- data$mi_x / data$I_11
      data$Dy1 <- data$mi_y / data$I_11
      data$Dx2 <- data$mi_x / data$I_22
      data$Dy2 <- data$mi_y / data$I_22
    }

    if (dim == 3){
      data$len <- nrow(graph)
      data$mi_x <- sum(graph[, 'X']) / length(graph[, 'X'])
      data$mi_y <- sum(graph[, 'Y']) / length(graph[, 'Y'])
      data$mi_z <- sum(graph[, 'Z']) / length(graph[, 'Z'])
      data$sqrt <- sqrt(data$mi_x^2 + data$mi_y^2 + data$mi_z^2)

      graph_centered <- graph
      graph_centered[, 'X'] <- graph_centered[, 'X'] - data$mi_x
      graph_centered[, 'Y'] <- graph_centered[, 'Y'] - data$mi_y
      graph_centered[, 'Z'] <- graph_centered[, 'Z'] - data$mi_z

      data$I_xx <- sum(graph_centered[, 'Y']^2 + graph_centered[, 'Z']^2) # Moments of inertia
      data$I_yy <- sum(graph_centered[, 'X']^2 + graph_centered[, 'Z']^2)
      data$I_zz <- sum(graph_centered[, 'X']^2 + graph_centered[, 'Y']^2)
      data$I_xy <- -sum(graph_centered[, 'X'] * graph_centered[, 'Y'])
      data$I_yx <- -sum(graph_centered[, 'Y'] * graph_centered[, 'X'])
      data$I_xz <- -sum(graph_centered[, 'X'] * graph_centered[, 'Z'])
      data$I_zx <- -sum(graph_centered[, 'Z'] * graph_centered[, 'X'])
      data$I_zy <- -sum(graph_centered[, 'Z'] * graph_centered[, 'Y'])
      data$I_yz <- -sum(graph_centered[, 'Y'] * graph_centered[, 'Z'])

      char_mi <- matrix(,3,3) # We need matrix to calculate eigenvalues
      char_mi[1, 1] <- data$I_xx
      char_mi[1, 2] <- data$I_xy
      char_mi[1, 3] <- data$I_xz
      char_mi[2, 1] <- data$I_yx
      char_mi[2, 2] <- data$I_yy
      char_mi[2, 3] <- data$I_yz
      char_mi[3, 1] <- data$I_zx
      char_mi[3, 2] <- data$I_zy
      char_mi[3, 3] <- data$I_zz

      data$I_11 <- eigen(char_mi)$values[1] # Principal moments of inertia
      data$I_22 <- eigen(char_mi)$values[2]
      data$I_33 <- eigen(char_mi)$values[3]

      data$Dx1 <- data$mi_x / data$I_11
      data$Dy1 <- data$mi_y / data$I_11
      data$Dz1 <- data$mi_z / data$I_11
      data$Dx2 <- data$mi_x / data$I_22
      data$Dy2 <- data$mi_y / data$I_22
      data$Dz2 <- data$mi_z / data$I_22
      data$Dx3 <- data$mi_x / data$I_33
      data$Dy3 <- data$mi_y / data$I_33
      data$Dz3 <- data$mi_z / data$I_33
    }

    dataset <- rbind(dataset, data)
  }

  if (genbank == FALSE) rownames(dataset) <- 1:length(seqs)
  if (genbank == TRUE) rownames(dataset) <- seqs

  return(as.data.frame(dataset))
}
