# Author: Marian Schmidt 
# Date: November 26th, 2017
# Purpose: To calculate the phylogenetic diversity of HNA and LNA

###################### Set the working directory 
setwd("/scratch/lsa_fluxm/marschmi/HNA_LNA/fasttree_Nov13_2017/")

###################### Load Libraries
library(phyloseq)
library(picante)
library(tidyr)
library(dplyr)


###################### Load Data
# Read in the data with 2 data objects 
    # 1. ALL_rel_physeq_1in3
    # 2. ALL_rel_physeq_5in10
load("ALL-physeq-for-phylo.RData")


# Create the OTU table for picante 
# 1in3
rel_1in3_otu <- matrix(otu_table(ALL_rel_physeq_1in3), nrow = nrow(otu_table(ALL_rel_physeq_1in3)))
rownames(rel_1in3_otu) <- sample_names(ALL_rel_physeq_1in3)
colnames(rel_1in3_otu) <- taxa_names(ALL_rel_physeq_1in3)

# 5in10
rel_5in10_otu <- matrix(otu_table(ALL_rel_physeq_5in10), nrow = nrow(otu_table(ALL_rel_physeq_5in10)))
rownames(rel_5in10_otu) <- sample_names(ALL_rel_physeq_5in10)
colnames(rel_5in10_otu) <- taxa_names(ALL_rel_physeq_5in10)


## Calculate input for SES_MPD  
# Convert the abundance data to standardized abundanced vegan function `decostand' , NOTE: method = "total"
rel_1in3_otu_decostand <- decostand(rel_1in3_otu, method = "total")
rel_5in10_otu_decostand <- decostand(rel_5in10_otu, method = "total")
# check total abundance in each sample
stopifnot(sum(apply(rel_1in3_otu_decostand, 1, sum)) == nsamples(rel_1in3_otu))
stopifnot(sum(apply(rel_5in10_otu_decostand, 1, sum)) == nsamples(rel_5in10_otu))

# check for mismatches/missing species between community data and phylo tree
rel_1in3_otu_decostand_matches <- match.phylo.comm(phy_tree(ALL_rel_physeq_1in3), rel_1in3_otu_decostand)
rel_5in10_otu_decostand_matches <- match.phylo.comm(phy_tree(ALL_rel_physeq_5in10), rel_5in10_otu_decostand)


# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
phy_1in3 <- rel_1in3_otu_decostand_matches$phy
comm_1in3 <- rel_1in3_otu_decostand_matches$comm

phy_5in10 <- rel_5in10_otu_decostand_matches$phy
comm_5in10 <- rel_5in10_otu_decostand_matches$comm

# Calculate the phylogenetic distances
phy_dist_1in3 <- cophenetic(phy_1in3)
phy_dist_5in10 <- cophenetic(phy_5in10)

## Calculate SES_MPD
###################################### INDEPENDENT SWAP ############################################
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap_1in3 <- ses.mpd(comm_1in3, phy_dist_1in3, null.model = "independentswap", 
                                                    abundance.weighted = FALSE, runs = 999)

unweighted_sesMPD_indepswap_5in10 <- ses.mpd(comm_5in10, phy_dist_5in10, null.model = "independentswap", 
                                             abundance.weighted = FALSE, runs = 999)

##### Write the data to a file 
write.table(unweighted_sesMPD_indepswap_1in3,
            file = "unweighted_MPD_1seqs_in_3samps.tsv", row.names=TRUE)

write.table(unweighted_sesMPD_indepswap_5in10,
            file = "unweighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)