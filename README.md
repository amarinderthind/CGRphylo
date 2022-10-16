# CGRPhylo Pipeline: Chaos Game Representation for phylogeny

### Overview
Chaos Game Representation (CGR) is an iterative mapping technique to construct a two dimensional representation of genomic sequences (Jeffrey, 1990). 

CGRs have been conventionally used to visualize the large nucleotide sequences. However, apart from visualization, CGRs can be used to compare DNA sequences, construct cladograms and address various biological problems. It effectively classify very similar sequences (e.g., sub-subtypes of HIV-1 sequences) and deliver better resolved classification trees when compared to standard Maximum Likelihood methods. 

In this CGR package (Written in R), core module is integrated with various packages to provide a complete pipeline for multiple whole genome sequence comparison using Chaos Game Representation (CGR) based approach. CGR based distance matrices can be converted to phylogency tree. The most important idea for the development of this pipeline is its ability to handle large DNA sequences (As user machine).

It efficiently classifies sequences based on both inter-species and intra-species variation in a computationally less intense manner. It analyses whole genome variations using an alignment free and scale invariant method resulting in trees that can be used to interpret similarity between multiple whole genome sequences, even when they are closely related.

###
![figure1](https://user-images.githubusercontent.com/45668229/195962013-fef235d1-6987-4b98-bab9-7d6083f01e5e.png)

#### Steps of the pipeline

Input required is a set of two or more genome sequences in FASTA format. The other input required by user is “word length (K value)” for which the frequencies of all the words in the sequences are calculated. User can also specify the Out-group for construction of Neighbor Joining Tree. 

These Sequences are then used to calculate frequencies of words at user specified word length. , each sequence gets it's own frequency matrices. There is also an option to visualize the CGR plots. Based on these word frequencies the pair-wise distance (Euclidean-default,or  Square euclidean or Manhattan) is calculated between all pairs of sequences in input data. Package integrate other packages where distance matrices can be converted to NJ tree or other and visualized. We provide the option to convert distance matrix into MEGA and phylip distance formats, which is widely used program for visualization. 

##### Input file requirements

The DNA sequences to be analyzed can be uploaded in the form fasta file. All the sequences must be in FASTA format.  

#### Filtering and trimming, if required 


##### Selection of "Word Length"
The clustering of the sequences is based on the distances calculated from the frequencies of DNA words. The word length to be used for the calculation can be specified. This default word length used is 6.  

##### Visualization of CGR plot
CGRs for each sequence can be visualized by selecting the sequence.

###  Saving Results in different output formats

The results can be saved at different steps and different format for tree visulization with third party tools.

1) Distance Matrix file- It is (.txt) text file that contains pairwise distance matrix generated using input sequences. We used Euclidean distance method to calculate pairwise distances between multiple whole genome sequences. It contains output in standard PHYLIP distance matrix format. Following is the screenshot of the file for example data.


![tt](https://user-images.githubusercontent.com/45668229/195969269-dd01ab1c-d94b-4e52-9abc-bb8c14474524.png)

 
 2) Word Frequency file- It is a (.txt) text file that contains frequencies of all different k-letter words corresponding to CGR map.

For example – At k=3 in example data set

![image](https://user-images.githubusercontent.com/45668229/186616875-97dcc3aa-0d9d-4f1b-a0f9-4db0f96f8390.png)


![image](https://user-images.githubusercontent.com/45668229/186616980-55dcef85-6164-496f-86c1-add97badadcf.png)

3) Outtree,  users can save outtree files. These files contains information about trees generated in standard NEWICK format or other. These files are compatible with standard bioinformatics tools like TreeView , MEGA etc. to create, view , edit and customize trees.

Following example shows the SARS-CoV-2 different strains cladogram created using CGRPhylo.. 




 ![figure5 copy](https://user-images.githubusercontent.com/45668229/195969475-b84d614e-e43d-4e67-ba25-416a5a102f44.png)

