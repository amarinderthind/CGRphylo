# CGRPhylo Pipeline: Chaos Game Representation for phylogeny

A CGRPhylo pipeline combines the R core module with various packages to compare multiple whole genome sequences using Chaos Game Representation (CGR). CRG core function creates the frequencies object for each sequence which can be used to calculate distances among sequences. Later, CGR-based distance matrices can be converted to a phylogeny tree using neighbour-joining (NJ) or other methods. A major advantage of the CGRphylo R pipeline is its ability to handle large DNA sequences (per a user machine) and its effectiveness at classifying very similar sequences.

### How to start
Line-to-line Rscript is available in cgrPhlyo.r (and cgrPhylo.rmd) script. You can find out what to aspect in the o/p by following cgrPhlyo.pdf (or cgrPhylo.html).

### Background
Chaos Game Representation (CGR) is an iterative mapping technique to construct a two dimensional representation of genomic sequences (Jeffrey, 1990). CGRs have been conventionally used to visualize the large nucleotide sequences. However, apart from visualization, CGRs can be used to compare DNA sequences, construct cladograms and address various biological problems. It effectively classify very similar sequences (e.g., sub-subtypes of HIV-1 sequences) and deliver better resolved classification trees when compared to standard Maximum Likelihood methods. 

It efficiently classifies sequences based on both inter-species and intra-species variation in a computationally less intense manner. It analyses whole genome variations using an alignment free and scale invariant method resulting in trees that can be used to interpret similarity between multiple whole genome sequences, even when they are closely related.

###
![figure1](https://user-images.githubusercontent.com/45668229/195962013-fef235d1-6987-4b98-bab9-7d6083f01e5e.png)

 Word Frequency is a frequencies of all different k-letter words corresponding to CGR map. Following figure shows various K-letter words (above) and their calculated frequencies (below) at k=3 
 
 <p align="center">
 <img src="https://user-images.githubusercontent.com/45668229/186616875-97dcc3aa-0d9d-4f1b-a0f9-4db0f96f8390.png"  width="45%" height="400">&nbsp; &nbsp; &nbsp; &nbsp;
<img src="https://user-images.githubusercontent.com/45668229/186616980-55dcef85-6164-496f-86c1-add97badadcf.png"  width="45%" height="400">
</p>

##### Steps of the pipeline

Input required is a set of two or more genome sequences in FASTA format. The other input required by user is “word length (K value)” for which the frequencies of all the words in the sequences are calculated. User can also specify the Out-group for construction of Neighbor Joining Tree. 

These Sequences are then used to calculate frequencies of words at user specified word length. , each sequence gets it's own frequency matrices. There is also an option to visualize the CGR plots. Based on these word frequencies the pair-wise distance (Euclidean-default,or  Square euclidean or Manhattan) is calculated between all pairs of sequences in input data. Package integrate other packages where distance matrices can be converted to NJ tree or other and visualized. We provide the option to convert distance matrix into MEGA and phylip distance formats, which is widely used program for visualization. 

##### Input file requirements

The DNA sequences to be analyzed can be uploaded in the form fasta file. All the sequences must be in FASTA format. Check the example file for details.  

```
setwd(".")
source('cgat_function.r')
source('distances_n_other.r')

file <- "recom_3.fasta"

library("seqinr")
fastafile <- seqinr::read.fasta(file = file, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
```

##### Filtering and trimming, if required 
```
library(stringr)

N_filter <- 50  ## filter sequence with n bases > this value
fasta_filtered <- fastafile_new(fastafile, N_filter) ## create filtered sequence file

#write fasta file from filtered sequences
seqinr::write.fasta(sequences=fasta_filtered,names =names(fasta_filtered),file.out=paste("recombinant_XBB.1_Filter",N_filter,".fasta",sep = ''))
```

##### Sequence length and GC content /Meta info

`create_meta` function extract various types of information from the sequences and store into dataframe. 

```
meta <- create_meta(fastafile, N_filter) ## create seq features information
print(paste("std dev for seq length is",sd(meta$length),sep=" "))
print(paste("Median of the seq length is",median(meta$length),sep=" "))
print("Range of the seq length")
range(meta$length)

#boxplot(meta$length, ylab="Sequence length") ## overall

dotchart(meta$length, labels = meta$name, xlab = "Sequence length", pch = 21, bg = "green", pt.cex = 1, cex = 0.7)

dotchart(meta$GC_content, labels = meta$name, xlab = "GC content", pch = 21, bg = "green", pt.cex = 1, cex = 0.7)

```
 <p align="center">
<img src="https://user-images.githubusercontent.com/45668229/196595195-970288f9-d2b1-48fa-90de-aab72a43981f.jpg" width=45% height="400">&nbsp; &nbsp; &nbsp; &nbsp;
<img src="https://user-images.githubusercontent.com/45668229/196830695-c8702ce3-a2b8-4be0-b112-856e364b930a.png" width=45% height="400">
 
</p>

```
## box plot for each strains (In this example first part of the name is strain name)

meta$strains <- as.character(lapply(meta$name, function(x) strsplit(x, '_')[[1]][1])) ## split strains names

library(ggplot2)
ggplot(meta, aes(x = strains, y =length, color = strains )) + geom_boxplot()+ ylab("Sequence length (group level)")+coord_flip() 

#boxplot(log2(meta$length+1))
len_trim <- min(meta$length)

## Save meta
#write.csv(meta,paste("length_and_names",N_filter,".csv"))
```
 
<img src="https://user-images.githubusercontent.com/45668229/196326195-f5a3172d-78f1-4bc8-a225-51d80a993834.png" width="1000" height="400">


##### Visualization of CGR plot
CGRs for each sequence can be visualized by selecting the sequence. `cgrplot` function creates the 'x' and 'y' cordinates for each base pair (to plot on CRG plot).

```
source('cgrplot.r')
cgr1 <- cgrplot(1) ## enter the number of sequence from  "fasta_filtered"


plot(cgr1[,1],cgr1[,2], main=paste("CGR plot of", names(fasta_filtered)[1],sep=''),
     xlab = "", ylab = "", cex=0.2, pch = 4, frame = TRUE) 

#compare 2 cgr plots
cgr2 <- cgrplot(4)

library('RColorBrewer')
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))

plot(cgr1[,1],cgr1[,2], main=paste("CGR plot of ", names(fasta_filtered)[1],sep=''),
     xlab = "", ylab = "", cex.main=0.5, cex=0.4, pch = 4, frame = TRUE) 
plot(cgr2[,1],cgr2[,2], main=paste("CGR plot of ", names(fasta_filtered)[2],sep=''),
     xlab = "", ylab = "",cex.main=0.5,cex=0.4,pch = 4, frame = TRUE) 
```

![CGR_2plots](https://user-images.githubusercontent.com/45668229/196325788-e054df7d-2689-4e77-89c7-53c9f6797a6c.png)

##### Create frequency object for sequences for specific "Word Length"
The clustering of the sequences is based on the distances calculated from the frequencies of DNA words. The word length to be used for the calculation can be specified. This default word length used is 6.  `cgat` function do this job.

```
k_mer <- 6  ## define the value of K
Freq_mat_obj <- list()
sequence_new <- names(fasta_filtered) 

for(n in 1:length(fasta_filtered)) {   
  
  ##skip passing whole data to the function ##increase speed  ## "fasta_filtered" name is fixed
  
  Freq_mat_obj[[n]] <- cgat(k_mer,n, len_trim) # executing one seq at a time    #k-mer,seq_length,trimmed_length
  print(paste("processing sequence : ",n , sequence_new[n], sep=" "))
  
}

names(Freq_mat_obj) <- sequence_new
```
## Calculate distances 
User can use any of Euclidean (default),or  Square euclidean or Manhattan distance for the distance matrices. `matrixDistance` function takes inputs in following way for distance calculations.

```
j <- length(sequence_new)

distance <- matrix(0,j,j)  ## defining the empty matrix from distance calculations
row.names(distance) <- sequence_new[1:j]
colnames(distance) <- sequence_new[1:j]

for(n in 1:j) { 
  for(cr in 1:j) {
distance[n,cr] <- matrixDistance(Freq_mat_obj[[n]],Freq_mat_obj[[cr]],distance_type ='Euclidean' )  ##'Euclidean','S_Euclidean' and Manhattan
}}

distance <- round(distance, 5)
distance[1:5,1:3]

#plot(hclust(as.dist(distance)))
```
###  Saving Results in different output formats (mega, phylip, newick, Nexus)

For tree visulization with third party tools, results can be saved into mega or phylip format using `saveMegaDistance` and `savePhylipDistance` respectively.

```
Mega_file_name <- paste('Recombi_N_filtering_',N_filter,'_mega_distance_file_for_k_',k_mer,'.meg')

### save 
saveMegaDistance(Mega_file_name, distance) ## inputs are file name (with path) and distance 

PhylipFile_name <- paste('Phylip_N_filtering_',N_filter,'distances_k_',k_mer,'.txt')

 ## inputs are file name (with path) and distance 
savePhylipDistance(PhylipFile_name, distance, mode= 'original')  ## relaxed or original/other

##typical phylip allows 10 charcter for texa name ##new relaxed formart allows 250
## http://www.phylo.org/index.php/help/relaxed_phylip
```

Distance Matrix (as shown below) contains pairwise distance matrix generated using input sequences. By default CGRPhylo use Euclidean distance method to calculate pairwise distances between multiple whole genome sequences.


![tt](https://user-images.githubusercontent.com/45668229/195969269-dd01ab1c-d94b-4e52-9abc-bb8c14474524.png)

Creating and Saving tree into newick, Nexus format (pipeline use `ape` and `treeio` Bioconductor package)

Outtree,  users can save outtree files created by NJ method of Bioconductor package 'ape'. These files contains information about trees generated in standard NEWICK format or other. These files are compatible with standard bioinformatics tools like TreeView , MEGA etc. to create, view , edit and customize trees. 

```
tree <-write.tree(my_nj, file = "", append = FALSE, digits = 10, tree.names = FALSE)

ape::write.tree(my_nj, file='Newick_NJ_tree.txt') #for Newick format #write out trees in a text file
ape::write.nexus(my_nj, file='Nexus_NJ_tree.nex') ##for Nexus format


## reading tree from file 
 newick.tree <- ape::read.tree("Newick_NJ_tree.txt")  
 nexus.tree <- ape::read.nexus("Nexus_NJ_tree.nex")
```

Following tree is created this pipeline and later modified/label using mega software. 

 <p align="center">
<img src="https://user-images.githubusercontent.com/45668229/195969475-b84d614e-e43d-4e67-ba25-416a5a102f44.png" width=60% height=400>
</p>
