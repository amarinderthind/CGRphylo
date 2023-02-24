
cgat <- function(k_mer, seq_number, len_trim)
{

matDim<- sqrt(4^k_mer)    ##freq matrix dimension
mat <- matrix(0, matDim, matDim)  ## empty matrix 
 
x_base <- matDim/2  ## start x value
y_base <- matDim/2  ## start y value

for(i in 1:len_trim)
{ 
  base <- substr(fasta_filtered[[seq_number]], start = i, stop = i)
  
  if(base=='a'|| base=='A')
  {
    x_base=x_base/2
    y_base=y_base/2
  }else if(base=='t'||base=='T')
  {
    x_base=(x_base+matDim)/2
    y_base=y_base/2
  }else if(base=='g'||base=='G')
  {
    x_base=(x_base+matDim)/2
    y_base=(y_base+matDim)/2
  }  else if(base=='c'||base=='C')
  {
    x_base=x_base/2
    y_base=(y_base+matDim)/2
  }
  for(x in 1:matDim){
  for(y in 1:matDim){

        if(x_base > (x-1) & x_base < x){
          if(y_base> (y-1) & y_base< y){
          
                            mat[y,x] <- mat[y,x]+1 
        }}}}} 
  return(mat)
}

 
