############################
##### CGRPhylo Pipeline #####
#############################

setwd("./")
file <- "Input_recom_SARS_cov2.fasta"

library("seqinr")
fastafile <- seqinr::read.fasta(file = file, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)

source('cgat_function.r')
source('distances_n_other.r')
source('cgrplot.r')

##########################################################################
######## Filtering of the sequencing with a specified number of 'N'  bases
##########################################################################

library(stringr)

N_filter <- 50  ## filter sequence with n bases > this value
fasta_filtered <- fastafile_new(fastafile, N_filter) ## create filtered sequence file using fastafile_new function

#write fasta file from filtered sequences
seqinr::write.fasta(sequences=fasta_filtered,names =names(fasta_filtered),file.out=paste("recombinant_XBB.1_Filter",N_filter,".fasta",sep = ''))

#write fasta file from filtered sequences
#seqinr::write.fasta(sequences=fastafile_new,names =sequence_new,file.out =paste("Filtered_N",N_filter,"_new_seq.fasta",sep = ''))

#######################################################################################
############ Create frequency matrix object for each sequence; one by one #############  
######################################################################################

library(parallel)

k_mer <- 6  ## define the value of K
sequence_new <- names(fasta_filtered)
num_cores <- detectCores() - 1  # Use one less than max to avoid freezing system

# Define the wrapper function
process_sequence <- function(n) {
  seq_name <- sequence_new[n]
  message(paste("Processing sequence:", n, seq_name))
  result <- cgat(k_mer, n, len_trim)
  return(result)
}

# Run in parallel using mclapply (works on Linux/macOS)
Freq_mat_obj <- mclapply(1:length(fasta_filtered), process_sequence, mc.cores = num_cores)

# Assign names to the result list
names(Freq_mat_obj) <- sequence_new

#######################################################################################
############ Create frequency matrix object for each sequence; one by one ############# WITHOUT PARALLEL PROCESSING
######################################################################################

k_mer <- 6  ## define the value of K
Freq_mat_obj <- list()
sequence_new <- names(fastafile_new) 

source('cgat_function.r')

for(n in 1:length(fastafile_new)) {   
  
  ##skip passing whole data to the function ##increase speed  ## "fastafile_new" name is Locked 
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

