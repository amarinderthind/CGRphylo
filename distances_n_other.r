## function for creating distances
matrixDistance <- function( matrixA,   matrixB, distance_type='Euclidean')
{
  #print(paste("Applying",distance_type,"distance", sep=" "))
  distance <- 0
  for(i in 1:dim(matrixA)[1])
  {
    for(j in 1:dim(matrixA)[1])
    {
      if(distance_type =='Euclidean' || distance_type =='S_Euclidean' )
      {
        dist=(matrixA[i,j]-matrixB[i,j])^2
        
      }else   ##Manhattan
      {
        dist=abs(matrixA[i,j]-matrixB[i,j])
        
      } 
      
      distance=distance+dist    
    }
  }
  
  if(distance_type=='Euclidean')
  {
    distance <- sqrt(distance)
    return(distance)
    exit
    
  }else{ 
    #print("no square distance")
    return(distance)
    exit
  }
}
## function for creating metafile
create_meta <- function(fastafile,  N_filter =50)
{
  library(stringr) ##for str_count
  sequence <- names(fastafile) ##substr(names(fastafile),1,12)
  n_new <- 1
  meta <- data.frame(matrix(nrow = length(fastafile), ncol=7))
  fastafile_new <- list() #vector() 
  sequence_new <- vector()
  
  for(n in 1:length(fastafile)) {   
    
    if(str_count(fastafile[[n]], "n")<=N_filter){
      
      fastafile_new[[n_new]] <- fastafile[[n]] 
      meta[n_new,1] <- nchar(fastafile[[n]])
      meta[n_new,2] <- str_count(fastafile[[n]], "n")
      meta[n_new,3] <- str_count(fastafile[[n]], "a")
      meta[n_new,4] <- str_count(fastafile[[n]], "g")
      meta[n_new,5] <- str_count(fastafile[[n]], "t")
      meta[n_new,6] <- str_count(fastafile[[n]], "c")
      meta[n_new,7] <-   (meta[n_new,4]+meta[n_new,6])/meta[n_new,1]*100
      sequence_new[n_new] <- sequence[n]  
      # print(paste("processing sequence : ",sequence[n],"Total length of the sequence : ",nchar(fastafile[[n]]), sep=" "))
      n_new <- n_new+1
      
    }else{ print(paste("Sequnce",sequence[n]," is filetred"), sep = '\t')}
  }
  names(fastafile_new) <- sequence_new
  meta <- meta[!(is.na(meta$X1)),] ##as.data.frame(t(as.data.frame(length_n)))
  colnames(meta) <-c("length","n_base","a_base","g_base","t_base","c_base","GC_content")
  meta$name <- sequence_new
  row.names(meta) <- sequence_new
  return(meta)
}
## function for creating  filtered fasta file 

fastafile_new <- function(fastafile,  N_filter =50)
{
  library(stringr) ##for str_count
  sequence <- names(fastafile) ##substr(names(fastafile),1,12)
  n_new <- 1
  fastafile_new <- list() #vector() 
  sequence_new <- vector()
  for(n in 1:length(fastafile)) {   
    
    if(str_count(fastafile[[n]], "n")<=N_filter){
      
      fastafile_new[[n_new]] <- fastafile[[n]] 
      sequence_new[n_new] <- sequence[n]  
      n_new <- n_new+1
      
    }
  }
  names(fastafile_new) <- sequence_new
  return(fastafile_new)
}

### function for mega distance file 

saveMegaDistance <- function(Mega_file_name, distance)
{
sequence_names <- colnames(distance)
dist_mega <- distance
dist_mega[upper.tri(dist_mega,diag=TRUE)] <- NA  #lower.tri

dist_mega <- sapply(as.data.frame(dist_mega), as.character)
dist_mega[is.na(dist_mega)] <- " "

lin1 <- paste0("#mega","\n","!Title: Concatenated Files;", "\n","!Format DataType=Distance DataFormat=LowerLeft;","\n", sep = '')
write.table(lin1,Mega_file_name, col.names = FALSE, row.names = FALSE,quote = FALSE)
#write.table(paste('#',names1$V1,names$V2,names$V3, sep = ''),file_name, sep="\t",append = TRUE, col.names = FALSE, row.names = FALSE,quote = FALSE)
write.table(paste('#',sequence_names, sep = ''),Mega_file_name, sep="\t",append = TRUE, col.names = FALSE, row.names = FALSE,quote = FALSE)
write.table("\n",Mega_file_name, sep="\t",append = TRUE, col.names = FALSE, row.names = FALSE,quote = FALSE)
write.table(dist_mega,Mega_file_name, sep="\t", append = TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
}

savePhylipDistance <- function(PhylipFile_name, distance , mode='relaxed')
{
  sequence_names <- colnames(distance)
  dist_phylip <- distance
  dist_phylip[upper.tri(dist_phylip,diag=TRUE)] <- NA  #lower.tri
  
  dist_phylip <- sapply(as.data.frame(dist_phylip), as.character)
  dist_phylip[is.na(dist_phylip)] <- " "
  if(mode=='relaxed'){
  rownames(dist_phylip) <-   substr(sequence_names,1,10)}else{ rownames(dist_phylip) <-   substr(sequence_names,1,250)}
  
  write.table(length(sequence_names),PhylipFile_name, col.names = FALSE, row.names = FALSE,quote = FALSE)
  write.table(dist_phylip,PhylipFile_name, sep=" ", append = TRUE,col.names = FALSE, row.names = TRUE, quote = FALSE)
}
