
###################### Load Libraries
library(phyloseq)
library(picante)
library(tidyverse)

#setwd("~/git_repos/HNA_LNA_productivity/data/Chloroplasts_removed/ByLake_Filtering")

########################################################################################  
#################################### 1in3 ##############################################
########################################################################################  

# Read in the data 
# Read in the relative abundance data 
load("1in3/muskegon/muskegon_1in3_physeqs.RData")
rare_muskegon_physeq_1in3_rel
musk_otu_1in3 <- otu_table(rare_muskegon_physeq_1in3_rel)
musk_OTUnames_1in3 <- as.vector(colnames(otu_table(musk_otu_1in3)))
length(musk_OTUnames_1in3) # 1723

#write(musk_OTUnames_1in3, 
#      file = "1in3/muskegon/OTUnames_muskegon_1in3.txt",
#     ncolumns = 1,
#      append = FALSE, sep = "\n")


# Fix the taxonomy names
colnames(tax_table(rare_muskegon_physeq_1in3_rel)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

phy <- data.frame(tax_table(rare_muskegon_physeq_1in3_rel))
phy$OTU <- row.names(phy)
t <- tax_table(as.matrix(phy))
tax_table(rare_muskegon_physeq_1in3_rel) <- t

################################################  HNA ANALYSIS 1seq_in_3samples
musk_HNA_1in3 <- read.csv("../../../Scores/hnascores_otus_tuned_thr_0.23_1seq3_rel_Muskegon.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 1in3 OTU names
musk_HNA_otus_1in3 <- musk_HNA_1in3$OTU
# How many OTUs?
length(musk_HNA_otus_1in3) # 30

musk_HNA_physeq_1in3 <- subset_taxa(rare_muskegon_physeq_1in3_rel,  OTU %in% musk_HNA_otus_1in3)



################################################  LNA ANALYSIS 1seq_in_3samples
musk_LNA_1in3 <- read.csv("../../../Scores/lnascores_otus_tuned_thr_0.25_1seq3_rel_Muskegon.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 1in3 OTU names
musk_LNA_otus_1in3 <- musk_LNA_1in3$OTU
# How many OTUs?
length(musk_LNA_otus_1in3) # 28

musk_LNA_physeq_1in3 <- subset_taxa(rare_muskegon_physeq_1in3_rel,  OTU %in% musk_LNA_otus_1in3)


################################################  
################################################  WRITE OUT THE PHYLOSEQ OBJECT 
#save(list=c("musk_HNA_physeq_1in3","musk_LNA_physeq_1in3"), file=paste0("/1in3/muskegon/HNA-LNA-physeq_1in3.RData"))
# This above line of code will create a single .RData object and with the 2 following phyloseq objects: 
# 1. musk_HNA_physeq_1in3
# 2. musk_LNA_physeq_1in3




########################################################################################  
###################################### 5in10 ###########################################
######################################################################################## 
# 5 in 10
load("5in10/muskegon/muskegon_5in10_physeqs.RData")
rare_muskegon_physeq_5in10_rel
musk_otu_5in10 <- otu_table(rare_muskegon_physeq_5in10_rel)
musk_OTUnames_5in10 <- as.vector(colnames(otu_table(musk_otu_5in10)))
length(musk_OTUnames_5in10) # 482

#write(musk_OTUnames_5in10, 
#      file = "5in10/muskegon/OTUnames_muskegon_5in10.txt",
#      ncolumns = 1,
#      append = FALSE, sep = "\n")

# Fix the taxonomy names
colnames(tax_table(rare_muskegon_physeq_5in10_rel)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

phy <- data.frame(tax_table(rare_muskegon_physeq_5in10_rel))
phy$OTU <- row.names(phy)
t <- tax_table(as.matrix(phy))
tax_table(rare_muskegon_physeq_5in10_rel) <- t

################################################  HNA ANALYSIS 1seq_in_3samples
musk_HNA_5in10 <- read.csv("../../../Scores/hnascores_otus_tuned_thr_0.26_5seq10_rel_Muskegon.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 5in10 OTU names
musk_HNA_otus_5in10 <- musk_HNA_5in10$OTU
# How many OTUs?
length(musk_HNA_otus_5in10) # 31

musk_HNA_physeq_5in10 <- subset_taxa(rare_muskegon_physeq_5in10_rel,  OTU %in% musk_HNA_otus_5in10)



################################################  LNA ANALYSIS 1seq_in_3samples
musk_LNA_5in10 <- read.csv("../../../Scores/lnascores_otus_tuned_thr_0.24_5seq10_rel_Muskegon.csv") %>%
  rename(OTU = X) 

# Create a vector of the HNA 5in10 OTU names
musk_LNA_otus_5in10 <- musk_LNA_5in10$OTU
# How many OTUs?
length(musk_LNA_otus_5in10) # 44

musk_LNA_physeq_5in10 <- subset_taxa(rare_muskegon_physeq_5in10_rel,  OTU %in% musk_LNA_otus_5in10)


################################################  
################################################  WRITE OUT THE PHYLOSEQ OBJECT 
#save(list=c("musk_HNA_physeq_5in10","musk_LNA_physeq_5in10"), file=paste0("/5in10/muskegon/HNA-LNA-physeq_5in10.RData"))
# This above line of code will create a single .RData object and with the 2 following phyloseq objects: 
# 1. musk_HNA_physeq_5in10
# 2. musk_LNA_physeq_5in10



################################################  
################################################  ONE BIG PHYLOSEQ WITH HNA AND LNA SAMPLES = Phylogenetic Analysis
### Make a phyloseq object with HNA and LNA "Samples" keeping the original taxonomy 
# This object will be for phylogenetic analysis!
rel_physeq_5in10 # Will keep tax_table and phy_tree but REPLACE otu_table

############### LNA
###  MAKE NEW OTU TABLE FOR LNA SAMPLES
musk_otu_LNA_5in10 <- otu_table(musk_LNA_physeq_5in10) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names") %>%
  # Make a new names column with _LNA at the end of each sample name
  mutate(names = paste(names, "_LNA", sep="")) 

############### HNA
###  MAKE NEW OTU TABLE FOR HNA SAMPLES
musk_otu_HNA_5in10 <- otu_table(musk_HNA_physeq_5in10) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names") %>%
  # Make a new names column with _LNA at the end of each sample name
  mutate(names = paste(names, "_HNA", sep="")) 


# Full join with original OTU table 
musk_otu_all_5in10 <- otu_table(rare_muskegon_physeq_5in10_rel) %>%
  data.frame() %>%
  # Make a new column from the rownames
  tibble::rownames_to_column("names")

## Combine the regular, HNA, and LNA otus all into the same otu table
musk_otu_table_ALL_5in10 <- full_join(musk_otu_all_5in10, musk_otu_LNA_5in10) %>%
  full_join(musk_otu_HNA_5in10) %>%
  # Replace all of the NAs with zeros
  mutate_all(funs(ifelse(is.na(.), 0, .))) %>%
  # Make the names to be rownames
  tibble::column_to_rownames("names")

# compatible with phyloseq
#### MUST BE A NUMBER AND NOT AN NA!
# Combine all the data together into one object
ALL_rel_physeq_5in10 <- merge_phyloseq(otu_table(musk_otu_table_ALL_5in10, taxa_are_rows = FALSE), 
                                       tax_table(as.matrix(tax_5seq_in_10perc)),
                                       phy_tree(tree_5seq_in_10perc))



################################################  
################################################  WRITE OUT THE PHYLOSEQ OBJECT 
save(list=c("ALL_rel_physeq_5in10","ALL_rel_physeq_1in3"), file=paste0("ALL-physeq-for-phylo.RData"))
# This above line of code will create a single .RData object and with the 2 following phyloseq objects: 
# 1. ALL_rel_physeq_5in10
# 2. ALL_rel_physeq_1in3



