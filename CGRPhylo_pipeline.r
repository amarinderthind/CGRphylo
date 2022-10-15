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

##########################################################################################
################### Charactersticks of the sequences #####################################
##########################################################################################
 
lk <- as.data.frame(t(as.data.frame(length_n)))
lk$name <- sequence_new
lk <-as.data.frame(lk)

sd(lk$V1)
boxplot(lk$V1)
#boxplot(log2(lk$V1+1))
median(lk$V1)
range(lk$V1)

len_trim <- min(lk$V1)

#write.csv(lk,paste("length_and_names",N_filter,".csv"))

#######################################################################################
############ Create frequency matrix object for each sequence; one by one #############
######################################################################################

k_mer <- 6  ## define the value of K
Freq_mat_obj <- list()
sequence_new <- names(fastafile_new) 

source('C:/Users/athind/OneDrive - University of Wollongong/CGAT_covid19/cgat_function.r')

for(n in 1:length(fastafile_new)) {   
  
  ##skip passing whole data to the function ##increase speed  ## "fastafile_new" name is fixed 
  Freq_mat_obj[[n]] <- cgat(k_mer,n, len_trim) # executing one seq at a time #k-mer,seq_length,trimmed_length
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

#####################################################
### converting distances to mega format #############
#####################################################

Mega_file_name <- paste('Test_new__TRIM_After_N_filtering_',N_filter,'_mega_distance_file_for_k_',k_mer,'.meg')

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

########################################################
########################################################

library(ape)
#?ape::nj

my_dist_mat2 <- as.dist(distance)
##Minimum evolution tree
plot(fastme.ols(my_dist_mat2))
#my_dist_mat2

my_nj <- ape::nj(my_dist_mat2)
plot(my_nj)
plot(my_nj, "unrooted")
plot(my_nj,type = "cladogram")

tree <-write.tree(my_nj, file = "", append = FALSE, digits = 10, tree.names = FALSE)
ape::write.tree(my_nj, file='filename.txt') #for Newick format
##ape::write.nexus(my_nj, file='filename.nex') ##for Nexus format


tree <- ape::read.tree("filename.txt")

#tree <- ape::read.nexus("filename.nex")
plot(my_nj)
ape::nodelabels(bg = c(3, 3, 2))


##library
##BiocManager::install("treeio")
library('treeio')
tree <- treeio::read.newick("filename.txt")
##tree <-treeio::read.nexus("filename.nex")

library("ggtree")
ggplot(my_nj, aes(x, y)) + geom_tree() + theme_tree()


##Rerooting tree
# library(tidyverse)
# #BiocManager::install("ggtree")
# library(ggtree)
# # build a ggplot with a geom_tree
# ggplot(tree) + geom_tree() + theme_tree()
# 
# # This is convenient shorthand
# ggtree(as.data.frame(distance))


