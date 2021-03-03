#' @title 2D-DYnamic Representation of the DNA/RNA Sequences
#' @param seq sequence given as a single string, e.g. "AGTTGAGGGAG"
#' @param type DNA or RNA

drep <- function(seq, type = 'DNA'){
  g=c(1, 0); c=c(0, 1); a=c(-1, 0); t=c(0, -1); u=c(0, -1)
  seq <- tolower(seq)
  seq <- str_split(seq, '')
  if (length(table(seq)) > 4){
    print('The number of different characters in the sequence is greater than 4. Your sequence is probably incomplete and some errors may occur.')
  }
  N=length(seq)
  xy = matrix(,N,2)
  v=c(0,0)

  if(type == 'DNA'){
    for (i in 1:N){ # Do the walk
      if (seq[i] == "a") v=v+a
      if (seq[i] == "c") v=v+c
      if (seq[i] == "g") v=v+g
      if (seq[i] == "t") v=v+t
      xy[i,]=v
    }
  }
  if (type== 'RNA'){
    for (i in 1:N){ # Do the walk
      if (seq[i] == "a") v=v+a
      if (seq[i] == "c") v=v+c
      if (seq[i] == "g") v=v+g
      if (seq[i] == "u") v=v+u
      xy[i,]=v
    }
  }
  #else {
  #  stop('Incorrect type of string.')
  #}
  return (xy)
}
