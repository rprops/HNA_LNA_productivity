# Author: Marian Schmidt 
# Date: November 13th, 2017

###################### Set the working directory 
#setwd("/scratch/lsa_fluxm/marschmi/HNA_LNA/fasttree_Nov13_2017/")


# Create the OTU table for picante 
HNA_5in10_otu <- matrix(otu_table(HNA_physeq_5in10), nrow = nrow(otu_table(HNA_physeq_5in10)))
rownames(HNA_5in10_otu) <- sample_names(HNA_physeq_5in10)
colnames(HNA_5in10_otu) <- taxa_names(HNA_physeq_5in10)


## Calculate input for SES_MPD  
# Convert the abundance data to standardized abundanced vegan function `decostand' , NOTE: method = "total"
HNA_5in10_otu_decostand <- decostand(HNA_5in10_otu, method = "total")
# check total abundance in each sample
apply(HNA_5in10_otu_decostand, 1, sum)

# check for mismatches/missing species between community data and phylo tree
HNA_5in10_otu_decostand_matches <- match.phylo.comm(phy_tree(HNA_physeq_5in10), HNA_5in10_otu_decostand)
# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
HNA_phy_5in10 <- HNA_5in10_otu_decostand_matches$phy
HNA_comm_5in10 <- HNA_5in10_otu_decostand_matches$comm

# Calculate the phylogenetic distances
HNA_phy_dist_5in10 <- cophenetic(HNA_phy_5in10)

## Calculate SES_MPD
###################################### INDEPENDENT SWAP ############################################
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap_5in10_HNA <- ses.mpd(HNA_comm_5in10, HNA_phy_dist_5in10, null.model = "independentswap", 
                                                    abundance.weighted = FALSE, runs = 999)

WEIGHTED_sesMPD_indepswap_5in10_HNA <- ses.mpd(HNA_comm_5in10, HNA_phy_dist_5in10, null.model = "independentswap", 
                                                  abundance.weighted = TRUE, runs = 999)



################################################  LNA ANALYSIS 5seqs_in10perc
LNA_5in10 <- read.csv("../../../Scores/lnascores_otus_tuned_thr_0.21_5seq10_rel.csv") %>%
  rename(OTU = X)

LNA_otus <- LNA_5in10$OTU

LNA_physeq_5in10 <- subset_taxa(rel_physeq_5in10,  OTU %in% LNA_otus)

# Create the OTU table for picante 
LNA_5in10_otu <- matrix(otu_table(LNA_physeq_5in10), nrow = nrow(otu_table(LNA_physeq_5in10)))
rownames(LNA_5in10_otu) <- sample_names(LNA_physeq_5in10)
colnames(LNA_5in10_otu) <- taxa_names(LNA_physeq_5in10)


## Calculate input for SES_MPD  
# Convert the abundance data to standardized abundanced vegan function `decostand' , NOTE: method = "total"
LNA_5in10_otu_decostand <- decostand(LNA_5in10_otu, method = "total")
# check total abundance in each sample
apply(LNA_5in10_otu_decostand, 1, sum)

# check for mismatches/missing species between community data and phylo tree
LNA_5in10_otu_decostand_matches <- match.phylo.comm(phy_tree(LNA_physeq_5in10), LNA_5in10_otu_decostand)
# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
LNA_phy_5in10 <- LNA_5in10_otu_decostand_matches$phy
LNA_comm_5in10 <- LNA_5in10_otu_decostand_matches$comm

# Calculate the phylogenetic distances
LNA_phy_dist_5in10 <- cophenetic(LNA_phy_5in10)

## Calculate SES_MPD
###################################### INDEPENDENT SWAP ############################################
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap_5in10_LNA <- ses.mpd(LNA_comm_5in10, LNA_phy_dist_5in10, null.model = "independentswap", 
                                                 abundance.weighted = FALSE, runs = 999)

WEIGHTED_sesMPD_indepswap_5in10_LNA <- ses.mpd(LNA_comm_5in10, LNA_phy_dist_5in10, null.model = "independentswap", 
                                               abundance.weighted = TRUE, runs = 999)

##########  PLOT
unweighted_sesMPD_indepswap_5in10_LNA <- 
  unweighted_sesMPD_indepswap_5in10_LNA %>%
  mutate(type = "LNA")

ARG <- unweighted_sesMPD_indepswap_5in10_HNA %>%
  mutate(type = "HNA") %>% 
  bind_rows(unweighted_sesMPD_indepswap_5in10_LNA)


ggplot(ARG, aes(y = mpd.obs.z, x = type, color = type)) + 
  geom_boxplot(alpha = 0.3, color = "black", aes(fill = type)) +
  theme_classic() +
  geom_jitter() + ggtitle("5 Seqs in 10% of Samples") +
  theme(axis.title.x = element_blank())




read.csv("../../../Scores/lnascores_otus_abun_1seq3_rel.csv") %>%
  rename(OTU = X) %>%
  filter(score > 0.5) %>%
  nrow()










##### Write the data to a file 
write.table(unweighted_sesMPD_indepswap_5in10,
            file = "5seqs_in_10percent_samples/unweighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)
write.table(WEIGHTED_sesMPD_indepswap_5in10,
            file = "5seqs_in_10percent_samples/weighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)

















######################################################################################## 


# Fix the taxonomy names
colnames(tax_table(rel_physeq_1in3)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

################################################
########## ADD THE PROTEOBACTERIA TO THE PHYLA
phy <- data.frame(tax_table(rel_physeq_1in3))
Phylum <- as.character(phy$Phylum)
Class <- as.character(phy$Class)

for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    Phylum[i] <- Class[i]
  } 
}

phy$Phylum <- Phylum
phy$OTU <- row.names(phy)

t <- tax_table(as.matrix(phy))

tax_table(rel_physeq_1in3) <- t
################################################



################################################  HNA ANALYSIS 5seqs_in10perc
HNA_5in10 <- read.csv("../../../Scores/hnascores_otus_tuned_thr_0.36_5seq10_rel.csv") %>%
  rename(OTU = X)

HNA_otus <- HNA_5in10$OTU

HNA_physeq_5in10 <- subset_taxa(rel_physeq_5in10,  OTU %in% HNA_otus)

# Create the OTU table for picante 
HNA_5in10_otu <- matrix(otu_table(HNA_physeq_5in10), nrow = nrow(otu_table(HNA_physeq_5in10)))
rownames(HNA_5in10_otu) <- sample_names(HNA_physeq_5in10)
colnames(HNA_5in10_otu) <- taxa_names(HNA_physeq_5in10)


## Calculate input for SES_MPD  
# Convert the abundance data to standardized abundanced vegan function `decostand' , NOTE: method = "total"
HNA_5in10_otu_decostand <- decostand(HNA_5in10_otu, method = "total")
# check total abundance in each sample
apply(HNA_5in10_otu_decostand, 1, sum)

# check for mismatches/missing species between community data and phylo tree
HNA_5in10_otu_decostand_matches <- match.phylo.comm(phy_tree(HNA_physeq_5in10), HNA_5in10_otu_decostand)
# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
HNA_phy_5in10 <- HNA_5in10_otu_decostand_matches$phy
HNA_comm_5in10 <- HNA_5in10_otu_decostand_matches$comm

# Calculate the phylogenetic distances
HNA_phy_dist_5in10 <- cophenetic(HNA_phy_5in10)

## Calculate SES_MPD
###################################### INDEPENDENT SWAP ############################################
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap_5in10_HNA <- ses.mpd(HNA_comm_5in10, HNA_phy_dist_5in10, null.model = "independentswap", 
                                                 abundance.weighted = FALSE, runs = 999)

WEIGHTED_sesMPD_indepswap_5in10_HNA <- ses.mpd(HNA_comm_5in10, HNA_phy_dist_5in10, null.model = "independentswap", 
                                               abundance.weighted = TRUE, runs = 999)



################################################  LNA ANALYSIS 5seqs_in10perc
LNA_5in10 <- read.csv("../../../Scores/lnascores_otus_tuned_thr_0.21_5seq10_rel.csv") %>%
  rename(OTU = X)

LNA_otus <- LNA_5in10$OTU

LNA_physeq_5in10 <- subset_taxa(rel_physeq_5in10,  OTU %in% LNA_otus)

# Create the OTU table for picante 
LNA_5in10_otu <- matrix(otu_table(LNA_physeq_5in10), nrow = nrow(otu_table(LNA_physeq_5in10)))
rownames(LNA_5in10_otu) <- sample_names(LNA_physeq_5in10)
colnames(LNA_5in10_otu) <- taxa_names(LNA_physeq_5in10)


## Calculate input for SES_MPD  
# Convert the abundance data to standardized abundanced vegan function `decostand' , NOTE: method = "total"
LNA_5in10_otu_decostand <- decostand(LNA_5in10_otu, method = "total")
# check total abundance in each sample
apply(LNA_5in10_otu_decostand, 1, sum)

# check for mismatches/missing species between community data and phylo tree
LNA_5in10_otu_decostand_matches <- match.phylo.comm(phy_tree(LNA_physeq_5in10), LNA_5in10_otu_decostand)
# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
LNA_phy_5in10 <- LNA_5in10_otu_decostand_matches$phy
LNA_comm_5in10 <- LNA_5in10_otu_decostand_matches$comm

# Calculate the phylogenetic distances
LNA_phy_dist_5in10 <- cophenetic(LNA_phy_5in10)

## Calculate SES_MPD
###################################### INDEPENDENT SWAP ############################################
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap_5in10_LNA <- ses.mpd(LNA_comm_5in10, LNA_phy_dist_5in10, null.model = "independentswap", 
                                                 abundance.weighted = FALSE, runs = 999)

WEIGHTED_sesMPD_indepswap_5in10_LNA <- ses.mpd(LNA_comm_5in10, LNA_phy_dist_5in10, null.model = "independentswap", 
                                               abundance.weighted = TRUE, runs = 999)

##########  PLOT
unweighted_sesMPD_indepswap_5in10_LNA <- 
  unweighted_sesMPD_indepswap_5in10_LNA %>%
  mutate(type = "LNA")

ARG <- unweighted_sesMPD_indepswap_5in10_HNA %>%
  mutate(type = "HNA") %>% 
  bind_rows(unweighted_sesMPD_indepswap_5in10_LNA)


ggplot(ARG, aes(y = mpd.obs.z, x = type, color = type)) + 
  geom_boxplot(alpha = 0.3, color = "black", aes(fill = type)) +
  theme_classic() +
  geom_jitter() + ggtitle("5 Seqs in 10% of Samples") +
  theme(axis.title.x = element_blank())

##### Write the data to a file 
write.table(unweighted_sesMPD_indepswap_1in3,
            file = "1seq_in_3samples/unweighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)
write.table(WEIGHTED_sesMPD_indepswap_1in3,
            file = "1seq_in_3samples/weighted_MPD_5seqs_in_10perc.tsv", row.names=TRUE)
