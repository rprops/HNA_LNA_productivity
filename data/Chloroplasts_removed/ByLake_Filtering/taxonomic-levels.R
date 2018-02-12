# Purpose: See how many levels within each taxonomic level.
# Date: February 9th, 2018

# set working directory
setwd("data/Chloroplasts_removed/ByLake_Filtering")

# Load packages
library(phyloseq)
library(tidyverse)

# Load Data
# Muskegon 5in10
load("5in10/muskegon/muskegon_5in10_physeqs.RData")
# Muskegon 1in3
load("1in3/muskegon/muskegon_1in3_physeqs.RData")
# Inland 5in10
load("5in10/inland/inland_5in10_physeqs.RData")
# Inland 1in3
load("1in3/inland/inland_1in3_physeqs.RData")
# michigan 5in10
load("5in10/michigan/michigan_5in10_physeqs.RData")
# michigan 1in3
load("1in3/michigan/michigan_1in3_physeqs.RData")

# Rank1 = Domain
# Rank2 = Phylum
# Rank3 = Class
# Rank4 = Order
# Rank5 = Family 
# Rank6 = Genus
# Rank7 = Species

# Ranks 
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")


################ NUMBER OF TAXA ####################
## MUSKEGON
musk_5in10 <- c(ntaxa(tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank1")),
  ntaxa(tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank2")),
  ntaxa(tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank3")),
  ntaxa(tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank4")),
  ntaxa(tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank5")),
  ntaxa(tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank6")),
  ntaxa(tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank7")),
  ntaxa(rare_muskegon_physeq_5in10_rel))

## MUSKEGON
musk_1in3 <- c(ntaxa(tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank1")),
                ntaxa(tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank2")),
                ntaxa(tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank3")),
                ntaxa(tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank4")),
                ntaxa(tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank5")),
                ntaxa(tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank6")),
                ntaxa(tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank7")),
                ntaxa(rare_muskegon_physeq_1in3_rel))

## INLAND
inland_5in10 <- c(ntaxa(tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank1")),
                  ntaxa(tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank2")),
                  ntaxa(tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank3")),
                  ntaxa(tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank4")),
                  ntaxa(tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank5")),
                  ntaxa(tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank6")),
                  ntaxa(tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank7")),
                  ntaxa(rare_inland_physeq_5in10_rel))

inland_1in3 <- c(ntaxa(tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank1")),
                  ntaxa(tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank2")),
                  ntaxa(tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank3")),
                  ntaxa(tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank4")),
                  ntaxa(tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank5")),
                  ntaxa(tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank6")),
                  ntaxa(tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank7")),
                  ntaxa(rare_inland_physeq_1in3_rel))

## MICHIGAN
mich_5in10 <- c(ntaxa(tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank1")),
                  ntaxa(tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank2")),
                  ntaxa(tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank3")),
                  ntaxa(tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank4")),
                  ntaxa(tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank5")),
                  ntaxa(tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank6")),
                  ntaxa(tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank7")),
                  ntaxa(rare_michigan_physeq_5in10_rel))

mich_1in3 <- c(ntaxa(tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank1")),
                 ntaxa(tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank2")),
                 ntaxa(tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank3")),
                 ntaxa(tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank4")),
                 ntaxa(tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank5")),
                 ntaxa(tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank6")),
                 ntaxa(tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank7")),
                 ntaxa(rare_michigan_physeq_1in3_rel))

# Combine into dataframe for each dataset
# MUSKEGON
dat_musk_5in10 <- data.frame(ranks, musk_5in10) %>%
  as.data.frame() %>%
  rename(TaxRank = ranks, Num_Unique = musk_5in10) %>%
  mutate(dataset = "5in10", lake = "Muskegon")

dat_musk_1in3 <- data.frame(ranks, musk_1in3) %>%
  as.data.frame() %>%
  rename(TaxRank = ranks, Num_Unique = musk_1in3) %>%
  mutate(dataset = "1in3", lake = "Muskegon")

# INLAND
dat_inland_5in10 <- data.frame(ranks, inland_5in10) %>%
  as.data.frame() %>%
  rename(TaxRank = ranks, Num_Unique = inland_5in10) %>%
  mutate(dataset = "5in10", lake = "Inland")

dat_inland_1in3 <- data.frame(ranks, inland_1in3) %>%
  as.data.frame() %>%
  rename(TaxRank = ranks, Num_Unique = inland_1in3) %>%
  mutate(dataset = "1in3", lake = "Inland")

# MICHIGAN
dat_michigan_5in10 <- data.frame(ranks, mich_5in10) %>%
  as.data.frame() %>%
  rename(TaxRank = ranks, Num_Unique = mich_5in10) %>%
  mutate(dataset = "5in10", lake = "Michigan")

dat_michigan_1in3 <- data.frame(ranks, mich_1in3) %>%
  as.data.frame() %>%
  rename(TaxRank = ranks, Num_Unique = mich_1in3) %>%
  mutate(dataset = "1in3", lake = "Michigan")


# Combine all datasets together into one df
test <- bind_rows(dat_musk_5in10, dat_musk_1in3, dat_inland_5in10, dat_inland_1in3, dat_michigan_5in10, dat_michigan_1in3) %>%
  mutate(TaxRank = factor(TaxRank, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")))

ggplot(test, aes(x = TaxRank, y = Num_Unique, group = lake)) +
  geom_line() + geom_point() +
  labs(x = "Taxonomic Ranking", y = "Number of Unique Values") + 
  facet_grid(dataset~lake, scales = "free_y") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

ggsave("nums_tax_rankings.jpg", width = 7, height = 5)




############### CREATE NEW PHYLOSEQ OBJECTS
set.seed(777)

rare_muskegon_physeq_5in10_rel
# Relative Abundance
species_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank7")
species_musk_physeq_5in10
# Write files 
write.table(otu_table(species_musk_physeq_5in10), file="5in10/muskegon/species/rel_muskegon_species_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(species_musk_physeq_5in10)[,-1], file="5in10/muskegon/species/rel_muskegon_species_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(species_musk_physeq_5in10), file="5in10/muskegon/species/rel_muskegon_species_taxonomy_5in10.tsv", row.names=TRUE)

# Genus Level: Relative Abundance
genus_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank6")
genus_musk_physeq_5in10
# Write files 
write.table(otu_table(genus_musk_physeq_5in10), file="5in10/muskegon/genus/rel_muskegon_genus_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(genus_musk_physeq_5in10)[,-1], file="5in10/muskegon/genus/rel_muskegon_genus_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(genus_musk_physeq_5in10), file="5in10/muskegon/genus/rel_muskegon_genus_taxonomy_5in10.tsv", row.names=TRUE)

# Family level: Relative Abundance
family_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank5")
family_musk_physeq_5in10
# Write files 
write.table(otu_table(family_musk_physeq_5in10), file="5in10/muskegon/family/rel_muskegon_family_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(family_musk_physeq_5in10)[,-1], file="5in10/muskegon/family/rel_muskegon_family_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(family_musk_physeq_5in10), file="5in10/muskegon/family/rel_muskegon_family_taxonomy_5in10.tsv", row.names=TRUE)

# Order level: Relative Abundance
order_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank4")
order_musk_physeq_5in10
# Write files 
write.table(otu_table(order_musk_physeq_5in10), file="5in10/muskegon/order/rel_muskegon_order_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(order_musk_physeq_5in10)[,-1], file="5in10/muskegon/order/rel_muskegon_order_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(order_musk_physeq_5in10), file="5in10/muskegon/order/rel_muskegon_order_taxonomy_5in10.tsv", row.names=TRUE)

# class level: Relative Abundance
class_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank3")
class_musk_physeq_5in10
# Write files 
write.table(otu_table(class_musk_physeq_5in10), file="5in10/muskegon/class/rel_muskegon_class_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(class_musk_physeq_5in10)[,-1], file="5in10/muskegon/class/rel_muskegon_class_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(class_musk_physeq_5in10), file="5in10/muskegon/class/rel_muskegon_class_taxonomy_5in10.tsv", row.names=TRUE)

# Phylum level: Relative Abundance
phylum_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank2")
phylum_musk_physeq_5in10
# Write files 
write.table(otu_table(phylum_musk_physeq_5in10), file="5in10/muskegon/phylum/rel_muskegon_phylum_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(phylum_musk_physeq_5in10)[,-1], file="5in10/muskegon/phylum/rel_muskegon_phylum_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(phylum_musk_physeq_5in10), file="5in10/muskegon/phylum/rel_muskegon_phylum_taxonomy_5in10.tsv", row.names=TRUE)


### SAVE ALL THE PHYSEQS
save(list=c("species_musk_physeq_5in10", "genus_musk_physeq_5in10","family_musk_physeq_5in10", 
            "order_musk_physeq_5in10", "class_musk_physeq_5in10","phylum_musk_physeq_5in10"), 
     file=paste0("5in10/muskegon/rel_tax_collapse_muskegon_5in10_physeqs.RData"))

##########
#########
##########
### 1in3: Write the tsv files
# Species: Relative Abundance
species_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank7")
species_musk_physeq_1in3
# Write files 
write.table(otu_table(species_musk_physeq_1in3), file="1in3/muskegon/species/rel_muskegon_species_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(species_musk_physeq_1in3)[,-1], file="1in3/muskegon/species/rel_muskegon_species_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(species_musk_physeq_1in3), file="1in3/muskegon/species/rel_muskegon_species_taxonomy_1in3.tsv", row.names=TRUE)

# Genus: Relative Abundance
genus_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank6")
genus_musk_physeq_1in3
# Write files 
write.table(otu_table(genus_musk_physeq_1in3), file="1in3/muskegon/genus/rel_muskegon_genus_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(genus_musk_physeq_1in3)[,-1], file="1in3/muskegon/genus/rel_muskegon_genus_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(genus_musk_physeq_1in3), file="1in3/muskegon/genus/rel_muskegon_genus_taxonomy_1in3.tsv", row.names=TRUE)

# Family level
# Relative Abundance
family_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank5")
family_musk_physeq_1in3
# Write files 
write.table(otu_table(family_musk_physeq_1in3), file="1in3/muskegon/family/rel_muskegon_family_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(family_musk_physeq_1in3)[,-1], file="1in3/muskegon/family/rel_muskegon_family_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(family_musk_physeq_1in3), file="1in3/muskegon/family/rel_muskegon_family_taxonomy_1in3.tsv", row.names=TRUE)


# Order level: Relative Abundance
order_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank4")
order_musk_physeq_1in3
# Write files 
write.table(otu_table(order_musk_physeq_1in3), file="1in3/muskegon/order/rel_muskegon_order_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(order_musk_physeq_1in3)[,-1], file="1in3/muskegon/order/rel_muskegon_order_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(order_musk_physeq_1in3), file="1in3/muskegon/order/rel_muskegon_order_taxonomy_1in3.tsv", row.names=TRUE)


# Class level: Relative Abundance
class_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank3")
class_musk_physeq_1in3
# Write files 
write.table(otu_table(class_musk_physeq_1in3), file="1in3/muskegon/class/rel_muskegon_class_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(class_musk_physeq_1in3)[,-1], file="1in3/muskegon/class/rel_muskegon_class_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(class_musk_physeq_1in3), file="1in3/muskegon/class/rel_muskegon_class_taxonomy_1in3.tsv", row.names=TRUE)

#Phylum level
# Relative Abundance
phylum_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank2")
phylum_musk_physeq_1in3
# Write files 
write.table(otu_table(phylum_musk_physeq_1in3), file="1in3/muskegon/phylum/rel_muskegon_phylum_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(phylum_musk_physeq_1in3)[,-1], file="1in3/muskegon/phylum/rel_muskegon_phylum_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(phylum_musk_physeq_1in3), file="1in3/muskegon/phylum/rel_muskegon_phylum_taxonomy_1in3.tsv", row.names=TRUE)


### SAVE ALL THE PHYSEQS
save(list=c("species_musk_physeq_1in3", "genus_musk_physeq_1in3","family_musk_physeq_1in3", 
            "order_musk_physeq_1in3", "class_musk_physeq_1in3","phylum_musk_physeq_1in3"), 
     file=paste0("1in3/muskegon/rel_tax_collapse_muskegon_1in3_physeqs.RData"))







#################################### INLAND #################################### 
#################################### INLAND #################################### 
#################################### 5in10 #################################### 
rare_inland_physeq_5in10_rel
# Species Level: Relative Abundance 
species_inland_physeq_5in10 <- tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank7")
species_inland_physeq_5in10
# Write files 
write.table(otu_table(species_inland_physeq_5in10), file="5in10/inland/species/rel_inland_species_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(species_inland_physeq_5in10)[,-1], file="5in10/inland/species/rel_inland_species_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(species_inland_physeq_5in10), file="5in10/inland/species/rel_inland_species_taxonomy_5in10.tsv", row.names=TRUE)

# Genus Level: Relative Abundance
genus_inland_physeq_5in10 <- tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank6")
genus_inland_physeq_5in10
# Write files 
write.table(otu_table(genus_inland_physeq_5in10), file="5in10/inland/genus/rel_inland_genus_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(genus_inland_physeq_5in10)[,-1], file="5in10/inland/genus/rel_inland_genus_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(genus_inland_physeq_5in10), file="5in10/inland/genus/rel_inland_genus_taxonomy_5in10.tsv", row.names=TRUE)

# Family level: Relative Abundance
family_inland_physeq_5in10 <- tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank5")
family_inland_physeq_5in10
# Write files 
write.table(otu_table(family_inland_physeq_5in10), file="5in10/inland/family/rel_inland_family_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(family_inland_physeq_5in10)[,-1], file="5in10/inland/family/rel_inland_family_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(family_inland_physeq_5in10), file="5in10/inland/family/rel_inland_family_taxonomy_5in10.tsv", row.names=TRUE)

# Order level: Relative Abundance
order_inland_physeq_5in10 <- tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank4")
order_inland_physeq_5in10
# Write files 
write.table(otu_table(order_inland_physeq_5in10), file="5in10/inland/order/rel_inland_order_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(order_inland_physeq_5in10)[,-1], file="5in10/inland/order/rel_inland_order_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(order_inland_physeq_5in10), file="5in10/inland/order/rel_inland_order_taxonomy_5in10.tsv", row.names=TRUE)

# class level: Relative Abundance
class_inland_physeq_5in10 <- tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank3")
class_inland_physeq_5in10
# Write files 
write.table(otu_table(class_inland_physeq_5in10), file="5in10/inland/class/rel_inland_class_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(class_inland_physeq_5in10)[,-1], file="5in10/inland/class/rel_inland_class_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(class_inland_physeq_5in10), file="5in10/inland/class/rel_inland_class_taxonomy_5in10.tsv", row.names=TRUE)

# Phylum level: Relative Abundance
phylum_inland_physeq_5in10 <- tax_glom(rare_inland_physeq_5in10_rel, taxrank = "Rank2")
phylum_inland_physeq_5in10
# Write files 
write.table(otu_table(phylum_inland_physeq_5in10), file="5in10/inland/phylum/rel_inland_phylum_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(phylum_inland_physeq_5in10)[,-1], file="5in10/inland/phylum/rel_inland_phylum_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(phylum_inland_physeq_5in10), file="5in10/inland/phylum/rel_inland_phylum_taxonomy_5in10.tsv", row.names=TRUE)



### SAVE ALL THE PHYSEQS
save(list=c("species_inland_physeq_5in10", "genus_inland_physeq_5in10","family_inland_physeq_5in10", 
            "order_inland_physeq_5in10", "class_inland_physeq_5in10","phylum_inland_physeq_5in10"), 
     file=paste0("5in10/inland/rel_tax_collapse_inland_5in10_physeqs.RData"))


#################################### 1in3 ####################################
rare_inland_physeq_1in3_rel
# Species Level: Relative Abundance 
species_inland_physeq_1in3 <- tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank7")
species_inland_physeq_1in3
# Write files 
write.table(otu_table(species_inland_physeq_1in3), file="1in3/inland/species/rel_inland_species_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(species_inland_physeq_1in3)[,-1], file="1in3/inland/species/rel_inland_species_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(species_inland_physeq_1in3), file="1in3/inland/species/rel_inland_species_taxonomy_1in3.tsv", row.names=TRUE)

# Genus Level: Relative Abundance
genus_inland_physeq_1in3 <- tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank6")
genus_inland_physeq_1in3
# Write files 
write.table(otu_table(genus_inland_physeq_1in3), file="1in3/inland/genus/rel_inland_genus_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(genus_inland_physeq_1in3)[,-1], file="1in3/inland/genus/rel_inland_genus_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(genus_inland_physeq_1in3), file="1in3/inland/genus/rel_inland_genus_taxonomy_1in3.tsv", row.names=TRUE)

# Family level: Relative Abundance
family_inland_physeq_1in3 <- tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank5")
family_inland_physeq_1in3
# Write files 
write.table(otu_table(family_inland_physeq_1in3), file="1in3/inland/family/rel_inland_family_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(family_inland_physeq_1in3)[,-1], file="1in3/inland/family/rel_inland_family_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(family_inland_physeq_1in3), file="1in3/inland/family/rel_inland_family_taxonomy_1in3.tsv", row.names=TRUE)

# Order level: Relative Abundance
order_inland_physeq_1in3 <- tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank4")
order_inland_physeq_1in3
# Write files 
write.table(otu_table(order_inland_physeq_1in3), file="1in3/inland/order/rel_inland_order_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(order_inland_physeq_1in3)[,-1], file="1in3/inland/order/rel_inland_order_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(order_inland_physeq_1in3), file="1in3/inland/order/rel_inland_order_taxonomy_1in3.tsv", row.names=TRUE)

# class level: Relative Abundance
class_inland_physeq_1in3 <- tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank3")
class_inland_physeq_1in3
# Write files 
write.table(otu_table(class_inland_physeq_1in3), file="1in3/inland/class/rel_inland_class_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(class_inland_physeq_1in3)[,-1], file="1in3/inland/class/rel_inland_class_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(class_inland_physeq_1in3), file="1in3/inland/class/rel_inland_class_taxonomy_1in3.tsv", row.names=TRUE)

# Phylum level: Relative Abundance
phylum_inland_physeq_1in3 <- tax_glom(rare_inland_physeq_1in3_rel, taxrank = "Rank2")
phylum_inland_physeq_1in3
# Write files 
write.table(otu_table(phylum_inland_physeq_1in3), file="1in3/inland/phylum/rel_inland_phylum_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(phylum_inland_physeq_1in3)[,-1], file="1in3/inland/phylum/rel_inland_phylum_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(phylum_inland_physeq_1in3), file="1in3/inland/phylum/rel_inland_phylum_taxonomy_1in3.tsv", row.names=TRUE)



### SAVE ALL THE PHYSEQS
save(list=c("species_inland_physeq_1in3", "genus_inland_physeq_1in3","family_inland_physeq_1in3", 
            "order_inland_physeq_1in3", "class_inland_physeq_1in3","phylum_inland_physeq_1in3"), 
     file=paste0("1in3/inland/rel_tax_collapse_inland_1in3_physeqs.RData"))






#################################### michigan #################################### 
#################################### michigan #################################### 
#################################### 5in10 #################################### 
rare_michigan_physeq_5in10_rel
# Species Level: Relative Abundance 
species_michigan_physeq_5in10 <- tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank7")
species_michigan_physeq_5in10
# Write files 
write.table(otu_table(species_michigan_physeq_5in10), file="5in10/michigan/species/rel_michigan_species_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(species_michigan_physeq_5in10)[,-1], file="5in10/michigan/species/rel_michigan_species_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(species_michigan_physeq_5in10), file="5in10/michigan/species/rel_michigan_species_taxonomy_5in10.tsv", row.names=TRUE)

# Genus Level: Relative Abundance
genus_michigan_physeq_5in10 <- tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank6")
genus_michigan_physeq_5in10
# Write files 
write.table(otu_table(genus_michigan_physeq_5in10), file="5in10/michigan/genus/rel_michigan_genus_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(genus_michigan_physeq_5in10)[,-1], file="5in10/michigan/genus/rel_michigan_genus_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(genus_michigan_physeq_5in10), file="5in10/michigan/genus/rel_michigan_genus_taxonomy_5in10.tsv", row.names=TRUE)

# Family level: Relative Abundance
family_michigan_physeq_5in10 <- tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank5")
family_michigan_physeq_5in10
# Write files 
write.table(otu_table(family_michigan_physeq_5in10), file="5in10/michigan/family/rel_michigan_family_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(family_michigan_physeq_5in10)[,-1], file="5in10/michigan/family/rel_michigan_family_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(family_michigan_physeq_5in10), file="5in10/michigan/family/rel_michigan_family_taxonomy_5in10.tsv", row.names=TRUE)

# Order level: Relative Abundance
order_michigan_physeq_5in10 <- tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank4")
order_michigan_physeq_5in10
# Write files 
write.table(otu_table(order_michigan_physeq_5in10), file="5in10/michigan/order/rel_michigan_order_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(order_michigan_physeq_5in10)[,-1], file="5in10/michigan/order/rel_michigan_order_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(order_michigan_physeq_5in10), file="5in10/michigan/order/rel_michigan_order_taxonomy_5in10.tsv", row.names=TRUE)

# class level: Relative Abundance
class_michigan_physeq_5in10 <- tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank3")
class_michigan_physeq_5in10
# Write files 
write.table(otu_table(class_michigan_physeq_5in10), file="5in10/michigan/class/rel_michigan_class_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(class_michigan_physeq_5in10)[,-1], file="5in10/michigan/class/rel_michigan_class_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(class_michigan_physeq_5in10), file="5in10/michigan/class/rel_michigan_class_taxonomy_5in10.tsv", row.names=TRUE)

# Phylum level: Relative Abundance
phylum_michigan_physeq_5in10 <- tax_glom(rare_michigan_physeq_5in10_rel, taxrank = "Rank2")
phylum_michigan_physeq_5in10
# Write files 
write.table(otu_table(phylum_michigan_physeq_5in10), file="5in10/michigan/phylum/rel_michigan_phylum_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(phylum_michigan_physeq_5in10)[,-1], file="5in10/michigan/phylum/rel_michigan_phylum_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(phylum_michigan_physeq_5in10), file="5in10/michigan/phylum/rel_michigan_phylum_taxonomy_5in10.tsv", row.names=TRUE)



### SAVE ALL THE PHYSEQS
save(list=c("species_michigan_physeq_5in10", "genus_michigan_physeq_5in10","family_michigan_physeq_5in10", 
            "order_michigan_physeq_5in10", "class_michigan_physeq_5in10","phylum_michigan_physeq_5in10"), 
     file=paste0("5in10/michigan/rel_tax_collapse_michigan_5in10_physeqs.RData"))


#################################### 1in3 ####################################
rare_michigan_physeq_1in3_rel
# Species Level: Relative Abundance 
species_michigan_physeq_1in3 <- tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank7")
species_michigan_physeq_1in3
# Write files 
write.table(otu_table(species_michigan_physeq_1in3), file="1in3/michigan/species/rel_michigan_species_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(species_michigan_physeq_1in3)[,-1], file="1in3/michigan/species/rel_michigan_species_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(species_michigan_physeq_1in3), file="1in3/michigan/species/rel_michigan_species_taxonomy_1in3.tsv", row.names=TRUE)

# Genus Level: Relative Abundance
genus_michigan_physeq_1in3 <- tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank6")
genus_michigan_physeq_1in3
# Write files 
write.table(otu_table(genus_michigan_physeq_1in3), file="1in3/michigan/genus/rel_michigan_genus_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(genus_michigan_physeq_1in3)[,-1], file="1in3/michigan/genus/rel_michigan_genus_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(genus_michigan_physeq_1in3), file="1in3/michigan/genus/rel_michigan_genus_taxonomy_1in3.tsv", row.names=TRUE)

# Family level: Relative Abundance
family_michigan_physeq_1in3 <- tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank5")
family_michigan_physeq_1in3
# Write files 
write.table(otu_table(family_michigan_physeq_1in3), file="1in3/michigan/family/rel_michigan_family_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(family_michigan_physeq_1in3)[,-1], file="1in3/michigan/family/rel_michigan_family_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(family_michigan_physeq_1in3), file="1in3/michigan/family/rel_michigan_family_taxonomy_1in3.tsv", row.names=TRUE)

# Order level: Relative Abundance
order_michigan_physeq_1in3 <- tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank4")
order_michigan_physeq_1in3
# Write files 
write.table(otu_table(order_michigan_physeq_1in3), file="1in3/michigan/order/rel_michigan_order_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(order_michigan_physeq_1in3)[,-1], file="1in3/michigan/order/rel_michigan_order_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(order_michigan_physeq_1in3), file="1in3/michigan/order/rel_michigan_order_taxonomy_1in3.tsv", row.names=TRUE)

# class level: Relative Abundance
class_michigan_physeq_1in3 <- tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank3")
class_michigan_physeq_1in3
# Write files 
write.table(otu_table(class_michigan_physeq_1in3), file="1in3/michigan/class/rel_michigan_class_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(class_michigan_physeq_1in3)[,-1], file="1in3/michigan/class/rel_michigan_class_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(class_michigan_physeq_1in3), file="1in3/michigan/class/rel_michigan_class_taxonomy_1in3.tsv", row.names=TRUE)

# Phylum level: Relative Abundance
phylum_michigan_physeq_1in3 <- tax_glom(rare_michigan_physeq_1in3_rel, taxrank = "Rank2")
phylum_michigan_physeq_1in3
# Write files 
write.table(otu_table(phylum_michigan_physeq_1in3), file="1in3/michigan/phylum/rel_michigan_phylum_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(phylum_michigan_physeq_1in3)[,-1], file="1in3/michigan/phylum/rel_michigan_phylum_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(phylum_michigan_physeq_1in3), file="1in3/michigan/phylum/rel_michigan_phylum_taxonomy_1in3.tsv", row.names=TRUE)



### SAVE ALL THE PHYSEQS
save(list=c("species_michigan_physeq_1in3", "genus_michigan_physeq_1in3","family_michigan_physeq_1in3", 
            "order_michigan_physeq_1in3", "class_michigan_physeq_1in3","phylum_michigan_physeq_1in3"), 
     file=paste0("1in3/michigan/rel_tax_collapse_michigan_1in3_physeqs.RData"))