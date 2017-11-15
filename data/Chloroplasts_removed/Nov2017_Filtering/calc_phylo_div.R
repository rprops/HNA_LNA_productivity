# Author: Marian Schmidt 
# Date: November 13th, 2017

###################### Set the working directory 
setwd("/scratch/lsa_fluxm/marschmi/HNA_LNA/fasttree_Nov13_2017/")

###################### Load Libraries
library(phyloseq)
library(picante)

###################### 5seqs_in_10percent_samples ##############################
# Read in the data 
rel_5seqs_in10perc <- read.table("5seqs_in_10percent_samples/nochloro_relative_5in10percent.tsv")
tax_5seqs_in10perc <- read.table("5seqs_in_10percent_samples/nochloro_taxonomy_otu_5in10percent.tsv")
tree_5seqs_in10perc <- read.tree(file = "5seqs_in_10percent_samples/newick_tree_5seqs_in_10perc_rmN.tre")

rel_physeq_5in10 <- merge_phyloseq(otu_table(rel_5seqs_in10perc, taxa_are_rows = FALSE), 
                               tax_table(as.matrix(tax_5seqs_in10perc)),
                               phy_tree(tree_5seqs_in10perc))


# Create the OTU table for picante 
rel_5in10_otu <- matrix(otu_table(rel_physeq_5in10), nrow = nrow(otu_table(rel_physeq_5in10)))
rownames(rel_5in10_otu) <- sample_names(rel_physeq_5in10)
colnames(rel_5in10_otu) <- taxa_names(rel_physeq_5in10)


## Calculate input for SES_MPD  
# Convert the abundance data to standardized abundanced vegan function `decostand' , NOTE: method = "total"
rel_5in10_otu_decostand <- decostand(rel_5in10_otu, method = "total")
# check total abundance in each sample
apply(rel_5in10_otu_decostand, 1, sum)

# check for mismatches/missing species between community data and phylo tree
rel_5in10_otu_decostand_matches <- match.phylo.comm(tree_5seqs_in10perc, rel_5in10_otu_decostand)
# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
phy_5in10 <- rel_5in10_otu_decostand_matches$phy
comm_5in10 <- rel_5in10_otu_decostand_matches$comm

# Calculate the phylogenetic distances
phy_dist_5in10 <- cophenetic(phy_5in10)


## Calculate SES_MPD
###################################### INDEPENDENT SWAP ############################################
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap_5in10 <- ses.mpd(comm_5in10, phy_dist_5in10, null.model = "independentswap", 
                                                    abundance.weighted = FALSE, runs = 999)

WEIGHTED_sesMPD_indepswap_5in10 <- ses.mpd(comm_5in10, phy_dist_5in10, null.model = "independentswap", 
                                                  abundance.weighted = TRUE, runs = 999)

##### Write the data to a file 
write.table(unweighted_sesMPD_indepswap_5in10,
            file = "5seqs_in_10percent_samples/unweighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)
write.table(WEIGHTED_sesMPD_indepswap_5in10,
            file = "5seqs_in_10percent_samples/weighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)


######################################################################################## 
############################### 1seq_in_3samples ######################################
######################################################################################## 

# Read in the data 
rel_1seq_in3samps <- read.table("1seq_in_3samples/nochloro_relative_1seqin3samps.tsv")
tax_1seq_in3samps <- read.table("1seq_in_3samples/nochloro_taxonomy_otu_1seqin3samps.tsv")
tree_1seq_in3samps <- read.tree(file = "1seq_in_3samples/newick_tree_1seqs_in_3samps_rmN.tre")

rel_physeq_1in3 <- merge_phyloseq(otu_table(rel_1seq_in3samps, taxa_are_rows = FALSE), 
                                   tax_table(as.matrix(tax_1seq_in3samps)),
                                   phy_tree(tree_1seq_in3samps))


# Create the OTU table for picante 
rel_1in3_otu <- matrix(otu_table(rel_physeq_1in3), nrow = nrow(otu_table(rel_physeq_1in3)))
rownames(rel_1in3_otu) <- sample_names(rel_physeq_1in3)
colnames(rel_1in3_otu) <- taxa_names(rel_physeq_1in3)

## Calculate input for SES_MPD  
# Convert the abundance data to standardized abundanced vegan function `decostand' , NOTE: method = "total"
rel_1in3_otu_decostand <- decostand(rel_1in3_otu, method = "total")
# check total abundance in each sample
apply(rel_1in3_otu_decostand, 1, sum)

# check for mismatches/missing species between community data and phylo tree
rel_1in3_otu_decostand_matches <- match.phylo.comm(tree_1seq_in3samps, rel_1in3_otu_decostand)
# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
phy_1in3 <- rel_1in3_otu_decostand_matches$phy
comm_1in3 <- rel_1in3_otu_decostand_matches$comm

# Calculate the phylogenetic distances
phy_dist_1in3 <- cophenetic(phy_1in3)


## Calculate SES_MPD
###################################### INDEPENDENT SWAP ############################################
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap_1in3 <- ses.mpd(comm_1in3, phy_dist_1in3, null.model = "independentswap", 
                                                      abundance.weighted = FALSE, runs = 999)

WEIGHTED_sesMPD_indepswap_1in3 <- ses.mpd(comm_1in3, phy_dist_1in3, null.model = "independentswap", 
                                                    abundance.weighted = TRUE, runs = 999)

##### Write the data to a file 
write.table(unweighted_sesMPD_indepswap_1in3,
            file = "1seq_in_3samples/unweighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)
write.table(WEIGHTED_sesMPD_indepswap_1in3,
            file = "1seq_in_3samples/weighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)
