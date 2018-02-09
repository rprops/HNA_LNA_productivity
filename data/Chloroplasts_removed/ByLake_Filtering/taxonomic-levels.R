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
rare_muskegon_physeq_5in10_rel
# Relative Abundance
species_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank7")
species_musk_physeq_5in10
# Write files 
write.table(otu_table(species_musk_physeq_5in10), file="5in10/muskegon/species/rel_muskegon_species_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(species_musk_physeq_5in10)[,-1], file="5in10/muskegon/species/rel_muskegon_species_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(species_musk_physeq_5in10), file="5in10/muskegon/species/rel_muskegon_species_taxonomy_5in10.tsv", row.names=TRUE)

# Absolute abundance
abs_species_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_abs, taxrank = "Rank7")
abs_species_musk_physeq_5in10
# Write files
write.table(otu_table(abs_species_musk_physeq_5in10), file="5in10/muskegon/species/abs_muskegon_species_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(abs_species_musk_physeq_5in10)[,-1], file="5in10/muskegon/species/abs_muskegon_species_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(abs_species_musk_physeq_5in10), file="5in10/muskegon/species/abs_muskegon_species_taxonomy_5in10.tsv", row.names=TRUE)

# Family level
# Relative Abundance
family_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank5")
family_musk_physeq_5in10
# Write files 
write.table(otu_table(family_musk_physeq_5in10), file="5in10/muskegon/family/rel_muskegon_family_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(family_musk_physeq_5in10)[,-1], file="5in10/muskegon/family/rel_muskegon_family_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(family_musk_physeq_5in10), file="5in10/muskegon/family/rel_muskegon_family_taxonomy_5in10.tsv", row.names=TRUE)

# Absolute abundance
abs_family_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_abs, taxrank = "Rank5")
abs_family_musk_physeq_5in10
# Write files
write.table(otu_table(abs_family_musk_physeq_5in10), file="5in10/muskegon/family/abs_muskegon_family_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(abs_family_musk_physeq_5in10)[,-1], file="5in10/muskegon/family/abs_muskegon_family_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(abs_family_musk_physeq_5in10), file="5in10/muskegon/family/abs_muskegon_family_taxonomy_5in10.tsv", row.names=TRUE)


#Phylum level
# Relative Abundance
phylum_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_rel, taxrank = "Rank2")
phylum_musk_physeq_5in10
# Write files 
write.table(otu_table(phylum_musk_physeq_5in10), file="5in10/muskegon/phylum/rel_muskegon_phylum_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(phylum_musk_physeq_5in10)[,-1], file="5in10/muskegon/phylum/rel_muskegon_phylum_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(phylum_musk_physeq_5in10), file="5in10/muskegon/phylum/rel_muskegon_phylum_taxonomy_5in10.tsv", row.names=TRUE)

# Absolute abundance
abs_phylum_musk_physeq_5in10 <- tax_glom(rare_muskegon_physeq_5in10_abs, taxrank = "Rank2")
abs_phylum_musk_physeq_5in10
# Write files
write.table(otu_table(abs_phylum_musk_physeq_5in10), file="5in10/muskegon/phylum/abs_muskegon_phylum_5in10_otu_5in10.tsv", row.names=TRUE)
write.table(sample_data(abs_phylum_musk_physeq_5in10)[,-1], file="5in10/muskegon/phylum/abs_muskegon_phylum_sampledata_5in10.tsv", row.names=TRUE)
write.table(tax_table(abs_phylum_musk_physeq_5in10), file="5in10/muskegon/phylum/abs_muskegon_phylum_taxonomy_5in10.tsv", row.names=TRUE)

### SAVE ALL THE PHYSEQS
save(list=c("species_musk_physeq_5in10","family_musk_physeq_5in10", "phylum_musk_physeq_5in10"), file=paste0("5in10/muskegon/rel_tax_collapse_muskegon_5in10_physeqs.RData"))
save(list=c("abs_species_musk_physeq_5in10","abs_family_musk_physeq_5in10", "abs_phylum_musk_physeq_5in10"), file=paste0("5in10/muskegon/abs_tax_collapse_muskegon_5in10_physeqs.RData"))



##########
#########
##########
### 1in3: Write the tsv files
# Species
# Relative Abundance
species_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank7")
species_musk_physeq_1in3
# Write files 
write.table(otu_table(species_musk_physeq_1in3), file="1in3/muskegon/species/rel_muskegon_species_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(species_musk_physeq_1in3)[,-1], file="1in3/muskegon/species/rel_muskegon_species_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(species_musk_physeq_1in3), file="1in3/muskegon/species/rel_muskegon_species_taxonomy_1in3.tsv", row.names=TRUE)

# Absolute abundance
abs_species_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_abs, taxrank = "Rank7")
abs_species_musk_physeq_1in3
# Write files
write.table(otu_table(abs_species_musk_physeq_1in3), file="1in3/muskegon/species/abs_muskegon_species_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(abs_species_musk_physeq_1in3)[,-1], file="1in3/muskegon/species/abs_muskegon_species_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(abs_species_musk_physeq_1in3), file="1in3/muskegon/species/abs_muskegon_species_taxonomy_1in3.tsv", row.names=TRUE)



# Family level
# Relative Abundance
family_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank5")
family_musk_physeq_1in3
# Write files 
write.table(otu_table(family_musk_physeq_1in3), file="1in3/muskegon/family/rel_muskegon_family_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(family_musk_physeq_1in3)[,-1], file="1in3/muskegon/family/rel_muskegon_family_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(family_musk_physeq_1in3), file="1in3/muskegon/family/rel_muskegon_family_taxonomy_1in3.tsv", row.names=TRUE)

# Absolute abundance
abs_family_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_abs, taxrank = "Rank5")
abs_family_musk_physeq_1in3
# Write files
write.table(otu_table(abs_family_musk_physeq_1in3), file="1in3/muskegon/family/abs_muskegon_family_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(abs_family_musk_physeq_1in3)[,-1], file="1in3/muskegon/family/abs_muskegon_family_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(abs_family_musk_physeq_1in3), file="1in3/muskegon/family/abs_muskegon_family_taxonomy_1in3.tsv", row.names=TRUE)


#Phylum level
# Relative Abundance
phylum_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_rel, taxrank = "Rank2")
phylum_musk_physeq_1in3
# Write files 
write.table(otu_table(phylum_musk_physeq_1in3), file="1in3/muskegon/phylum/rel_muskegon_phylum_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(phylum_musk_physeq_1in3)[,-1], file="1in3/muskegon/phylum/rel_muskegon_phylum_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(phylum_musk_physeq_1in3), file="1in3/muskegon/phylum/rel_muskegon_phylum_taxonomy_1in3.tsv", row.names=TRUE)

# Absolute abundance
abs_phylum_musk_physeq_1in3 <- tax_glom(rare_muskegon_physeq_1in3_abs, taxrank = "Rank2")
abs_phylum_musk_physeq_1in3
# Write files
write.table(otu_table(abs_phylum_musk_physeq_1in3), file="1in3/muskegon/phylum/abs_muskegon_phylum_1in3_otu_1in3.tsv", row.names=TRUE)
write.table(sample_data(abs_phylum_musk_physeq_1in3)[,-1], file="1in3/muskegon/phylum/abs_muskegon_phylum_sampledata_1in3.tsv", row.names=TRUE)
write.table(tax_table(abs_phylum_musk_physeq_1in3), file="1in3/muskegon/phylum/abs_muskegon_phylum_taxonomy_1in3.tsv", row.names=TRUE)

### SAVE ALL THE PHYSEQS
save(list=c("species_musk_physeq_1in3","family_musk_physeq_1in3", "phylum_musk_physeq_1in3"), file=paste0("1in3/muskegon/rel_tax_collapse_muskegon_1in3_physeqs.RData"))
save(list=c("abs_species_musk_physeq_1in3","abs_family_musk_physeq_1in3", "abs_phylum_musk_physeq_1in3"), file=paste0("1in3/muskegon/abs_tax_collapse_muskegon_1in3_physeqs.RData"))

