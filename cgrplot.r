 cgrplot <- function(seq_number)
{
  mat <- matrix(0, nchar(fasta_filtered[[seq_number]]), 2)  ## empty matrix 
  
  
  x_base <- 1/2  ## start x value
  y_base <- 1/2  ## start y value
  
  for(i in 1:nchar(fasta_filtered[[seq_number]]))
  { 
    base <- substr(fasta_filtered[[seq_number]], start = i, stop = i)
    
    if(base=='a'|| base=='A')
    {
      x_base=x_base/2
      y_base=y_base/2
    }else if(base=='t'||base=='T')
    {
      x_base=(x_base+1)/2
      y_base=y_base/2
    }else if(base=='g'||base=='G')
    {
      x_base=(x_base+1)/2
      y_base=(y_base+1)/2
    }  else if(base=='c'||base=='C')
    {
      x_base=x_base/2
      y_base=(y_base+1)/2
    }
    
    mat[i,1] <- x_base
    mat[i,2] <- y_base 
  }
  return(mat)
}


