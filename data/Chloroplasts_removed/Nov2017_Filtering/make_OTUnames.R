# Purpose: Create 2 files with lists of OTUs for filtering the fasta files for phylogenetic analysis
# Date: Monday, November 13th, 2017
# Author: Marian Schmidt 

setwd("~/data/Chloroplasts_removed/Nov2017_Filtering")

########### 1seq in 3 samples
dat_1seq_in_3samps <- read.table("1seq_in_3samples/nochloro_taxonomy_otu_1seqin3samps.tsv")
OTUnames_1seq_in_3samps <- as.vector(row.names(dat_1seq_in_3samps))
length(OTUnames_1seq_in_3samps) # 2951

write(OTUnames_1seq_in_3samps, 
      file = "1seq_in_3samples/OTUnames_1seq_in_3samps.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")



########### 5seqs in 10percent samples
dat_5seqs_in_10perc <- read.table("5seqs_in_10percent_samples/nochloro_taxonomy_otu_5in10percent.tsv")
OTUnames_5seqs_in_10perc <- as.vector(row.names(dat_5seqs_in_10perc))
length(OTUnames_5seqs_in_10perc) # 486

write(OTUnames_5seqs_in_10perc, 
      file = "5seqs_in_10percent_samples/OTUnames_5seqs_in_10perc.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")


# All OTUs in 5seqs_in_10percent_samples are present within 1seq_in_3samples!!!
length(intersect(OTUnames_1seq_in_3samps, OTUnames_5seqs_in_10perc))
