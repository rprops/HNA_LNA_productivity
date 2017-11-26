
###################### Load Libraries
library(phyloseq)
library(picante)
library(tidyverse)

########################################################################################  
################################ 1seq_in_3samples ######################################
########################################################################################  

# Read in the data 
# Read in the relative abundance data 
rel_1seq_in3samps <- read.table("1seq_in_3samples/nochloro_relative_1seqin3samps.tsv")
# Read in the taxonomy data 
tax_1seq_in3samps <- read.table("1seq_in_3samples/nochloro_taxonomy_otu_1seqin3samps.tsv")
# Read in the phylogenetic tree 
tree_1seq_in3samps <- read.tree(file = "1seq_in_3samples/newick_tree_1seqs_in_3samps_rmN.tre")

# Combine all the data together into one object
rel_physeq_1in3 <- merge_phyloseq(otu_table(rel_1seq_in3samps, taxa_are_rows = FALSE), 
                                  tax_table(as.matrix(tax_1seq_in3samps)),
                                  phy_tree(tree_1seq_in3samps))

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


################################################  HNA ANALYSIS 1seq_in_3samples
HNA_1in3 <- read.csv("../../../Scores/hnascores_otus_tuned_thr_0.15_1seq3_rel.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 1in3 OTU names
HNA_otus_1in3 <- HNA_1in3$OTU
# How many OTUs?
length(HNA_otus_1in3) # 187

HNA_physeq_1in3 <- subset_taxa(rel_physeq_1in3,  OTU %in% HNA_otus_1in3)



################################################  LNA ANALYSIS 1seq_in_3samples
LNA_1in3 <- read.csv("../../../Scores/lnascores_otus_tuned_thr_0.18_1seq3_rel.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 1in3 OTU names
LNA_otus_1in3 <- LNA_1in3$OTU
# How many OTUs?
length(LNA_otus_1in3) # 153

LNA_physeq_1in3 <- subset_taxa(rel_physeq_1in3,  OTU %in% LNA_otus_1in3)


################################################  
################################################  WRITE OUT THE PHYLOSEQ OBJECT 
save(list=c("HNA_physeq_1in3","LNA_physeq_1in3"), file=paste0("HNA-LNA-physeq_1in3.RData"))
  # This above line of code will create a single .RData object and with the 2 following phyloseq objects: 
    # 1. HNA_physeq_1in3
    # 2. LNA_physeq_1in3




################################################  
################################################  ONE BIG PHYLOSEQ WITH HNA AND LNA SAMPLES = Phylogenetic Analysis
### Make a phyloseq object with HNA and LNA "Samples" keeping the original taxonomy 
# This object will be for phylogenetic analysis!
rel_physeq_1in3 # Will keep tax_table and phy_tree but REPLACE otu_table

############### LNA
###  MAKE NEW OTU TABLE FOR LNA SAMPLES
otu_LNA_1in3 <- otu_table(LNA_physeq_1in3) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names") %>%
  # Make a new names column with _LNA at the end of each sample name
  mutate(names = paste(names, "_LNA", sep="")) 

############### HNA
###  MAKE NEW OTU TABLE FOR HNA SAMPLES
otu_HNA_1in3 <- otu_table(HNA_physeq_1in3) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names") %>%
  # Make a new names column with _LNA at the end of each sample name
  mutate(names = paste(names, "_HNA", sep="")) 


# Full join with original OTU table 
otu_all_1in3 <- otu_table(rel_physeq_1in3) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names")

## Combine the regular, HNA, and LNA otus all into the same otu table
otu_table_ALL_1in3 <- full_join(otu_all_1in3, otu_LNA_1in3) %>%
  full_join(otu_HNA_1in3) %>%
  # Replace all of the NAs with zeros
  mutate_all(funs(ifelse(is.na(.), 0, .))) %>%
  # Make the names to be rownames
  tibble::column_to_rownames("names")

# compatible with phyloseq
#### MUST BE A NUMBER AND NOT AN NA!
# Combine all the data together into one object
ALL_rel_physeq_1in3 <- merge_phyloseq(otu_table(otu_table_ALL_1in3, taxa_are_rows = FALSE), 
                                       tax_table(as.matrix(tax_1seq_in3samps)),
                                       phy_tree(tree_1seq_in3samps))
### See bottom of file for saving this phyloseq object!






########################################################################################  
############################### 5seq_in_10percent ######################################
########################################################################################  

# Read in the data 
# Read in the relative abundance data 
rel_5seq_in_10perc <- read.table("5seqs_in_10percent_samples/nochloro_relative_5in10percent.tsv")
# Read in the taxonomy data 
tax_5seq_in_10perc <- read.table("5seqs_in_10percent_samples/nochloro_taxonomy_otu_5in10percent.tsv")
# Read in the phylogenetic tree 
tree_5seq_in_10perc <- read.tree(file = "5seqs_in_10percent_samples/newick_tree_5seqs_in_10perc_rmN.tre")

# Combine all the data together into one object
rel_physeq_5in10 <- merge_phyloseq(otu_table(rel_5seq_in_10perc, taxa_are_rows = FALSE), 
                                  tax_table(as.matrix(tax_5seq_in_10perc)),
                                  phy_tree(tree_5seq_in_10perc))

# Fix the taxonomy names
colnames(tax_table(rel_physeq_5in10)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

################################################
########## ADD THE PROTEOBACTERIA TO THE PHYLA
phy <- data.frame(tax_table(rel_physeq_5in10))
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

tax_table(rel_physeq_5in10) <- t
################################################


################################################  HNA ANALYSIS 5seq_in_10percent
HNA_5in10 <- read.csv("../../../Scores/hnascores_otus_tuned_thr_0.36_5seq10_rel.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 1in3 OTU names
HNA_otus_5in10 <- HNA_5in10$OTU
# How many OTUs?
length(HNA_otus_5in10) # 41

HNA_physeq_5in10 <- subset_taxa(rel_physeq_5in10,  OTU %in% HNA_otus_5in10)



################################################  LNA ANALYSIS 5seq_in_10percent
LNA_5in10 <- read.csv("../../../Scores/lnascores_otus_tuned_thr_0.21_5seq10_rel.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 1in3 OTU names
LNA_otus_5in10 <- LNA_5in10$OTU
# How many OTUs?
length(LNA_otus_5in10) # 124

LNA_physeq_5in10 <- subset_taxa(rel_physeq_5in10,  OTU %in% LNA_otus_5in10)


################################################  
################################################  WRITE OUT THE PHYLOSEQ OBJECT 
save(list=c("HNA_physeq_5in10","LNA_physeq_5in10"), file=paste0("HNA-LNA-physeq_5in10.RData"))
  # This above line of code will create a single .RData object and with the 2 following phyloseq objects: 
        # 1. HNA_physeq_5in10
        # 2. LNA_physeq_5in10


################################################  
################################################  ONE BIG PHYLOSEQ WITH HNA AND LNA SAMPLES = Phylogenetic Analysis
### Make a phyloseq object with HNA and LNA "Samples" keeping the original taxonomy 
# This object will be for phylogenetic analysis!
rel_physeq_5in10 # Will keep tax_table and phy_tree but REPLACE otu_table

############### LNA
###  MAKE NEW OTU TABLE FOR LNA SAMPLES
otu_LNA_5in10 <- otu_table(LNA_physeq_5in10) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names") %>%
  # Make a new names column with _LNA at the end of each sample name
  mutate(names = paste(names, "_LNA", sep="")) 

############### HNA
###  MAKE NEW OTU TABLE FOR HNA SAMPLES
otu_HNA_5in10 <- otu_table(HNA_physeq_5in10) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names") %>%
  # Make a new names column with _LNA at the end of each sample name
  mutate(names = paste(names, "_HNA", sep="")) 


# Full join with original OTU table 
otu_all_5in10 <- otu_table(rel_physeq_5in10) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names")

## Combine the regular, HNA, and LNA otus all into the same otu table
otu_table_ALL_5in10 <- full_join(otu_all_5in10, otu_LNA_5in10) %>%
  full_join(otu_HNA_5in10) %>%
  # Replace all of the NAs with zeros
  mutate_all(funs(ifelse(is.na(.), 0, .))) %>%
  # Make the names to be rownames
  tibble::column_to_rownames("names")

# compatible with phyloseq
#### MUST BE A NUMBER AND NOT AN NA!
# Combine all the data together into one object
ALL_rel_physeq_5in10 <- merge_phyloseq(otu_table(otu_table_ALL_5in10, taxa_are_rows = FALSE), 
                                   tax_table(as.matrix(tax_5seq_in_10perc)),
                                   phy_tree(tree_5seq_in_10perc))



################################################  
################################################  WRITE OUT THE PHYLOSEQ OBJECT 
save(list=c("ALL_rel_physeq_5in10","ALL_rel_physeq_1in3"), file=paste0("ALL-physeq-for-phylo.RData"))
# This above line of code will create a single .RData object and with the 2 following phyloseq objects: 
    # 1. ALL_rel_physeq_5in10
    # 2. ALL_rel_physeq_1in3




########################################################################################  
###################################### FIN #############################################
########################################################################################  



