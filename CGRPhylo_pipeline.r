############################
##### CGRPhylo Pipeline #####
#############################

setwd("./")
file <- "example.fasta"

library("seqinr")
fastafile <- seqinr::read.fasta(file = file, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)

source('cgat_function.r') 
source('distances.r')


##########################################################################
######## Filtering of the sequencing with specified number of 'N'  bases
##########################################################################

library(stringr) ##for str_count

sequence <- names(fastafile) ##substr(names(fastafile),1,12)
N_filter <- 50
n_new <- 1
fastafile_new <- list() 
sequence_new <- vector()
length_n <- list()

for(n in 1:length(fastafile)) {   
  
  if(str_count(fastafile[[n]], "n")<=N_filter){
    
    fastafile_new[[n_new]] <- fastafile[[n]] 
    sequence_new[n_new] <- sequence[n]  
    print(paste("processing sequence : ",sequence[n],"Total length of the sequence : ",nchar(fastafile[[n]]), sep=" "))
    print(paste(" A: ",str_count(fastafile[[n]], "a")," N: ",str_count(fastafile[[n]], "n")," G: ",str_count(fastafile[[n]], "g"), " C: ",str_count(fastafile[[n]], "c"),sep=" "))
    
    length_n[n_new] <- nchar(fastafile[[n]])
    n_new <- n_new+1
    
  }else{ print(paste("Filtering ",sequence[n]), sep = '\t')}
}
names(fastafile_new) <- sequence_new

#write fasta file from filtered sequences
#seqinr::write.fasta(sequences=fastafile_new,names =sequence_new,file.out =paste("Filtered_N",N_filter,"_new_seq.fasta",sep = ''))

 
#######################################################################################
############ Create frequency matrix object for each sequence; one by one #############
######################################################################################

k_mer <- 6  ## define the value of K
Freq_mat_obj <- list()
sequence_new <- names(fastafile_new) 

source('cgat_function.r')

for(n in 1:length(fastafile_new)) {   
  
  ##skip passing whole data to the function ##increase speed  ## "fastafile_new" name is fixed 
  Freq_mat_obj[[n]] <- cgat(k_mer,n, nchar(fastafile[[n]])) # executing one seq at a time #k-mer,seq_length,trimmed_length
  print(paste("processing sequence : ",n , sequence_new[n], sep=" "))
  
}
names(Freq_mat_obj) <- sequence_new

#################################################################################
############### calculates distances ############################################
################################################################################

j <- length(sequence_new)

distance <- matrix(0,j,j)  ## defining the empty matrix from distance calculations
row.names(distance) <- sequence_new[1:j]
colnames(distance) <- sequence_new[1:j]

for(n in 1:j) { 
  for(cr in 1:j) {
distance[n,cr] <- matrixDistance(Freq_used[[n]],Freq_used[[cr]],distance_type ='Euclidean' )  ##'Euclidean','S_Euclidean' and Manhattan
}}

distance <- round(distance, 5)

#plot(hclust(as.dist(distance)))

