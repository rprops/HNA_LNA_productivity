# Figure 4
# Feb 28, 2018

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Set HNA and LNA discrete colors
fcm_colors <- c(
  "HNA" = "deepskyblue4",
  "LNA" = "darkgoldenrod1",
  "Total" = "black")

# Read in Data
HNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "HNA", 
         RL.ranking = 1/RL.ranking)

LNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "LNA",
         RL.ranking = -1/RL.ranking)

scores_df <- bind_rows(HNA, LNA) %>%
  dplyr::select(OTU, RL.ranking, FCM_type)

ggplot(scores_df, aes(FCM_type, OTU)) + 
  geom_tile(aes(fill = RL.ranking),colour = "white") + 
  scale_fill_gradient2(low = "darkgoldenrod1", high = "deepskyblue4", mid = "white", 
                      midpoint = 0, na.value = "black") +
  theme(axis.title.x = element_blank())


#### Attempt 2
ggplot(dfscores, aes(y=RL.ranking, x=OTU, fill=FCM_type)) + 
  geom_bar(stat="identity", position="identity") + coord_flip() + #ggtitle("Summed OTUs") +
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  theme_classic() +
  ylab("Inverse RL Ranking") +
  scale_fill_manual(values = fcm_colors) +
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.95))
ggsave("heatmap_figs/muskegon_bar_0.15.jpg", width = 4, height = 8)



#### Attempt 3
# Read in Data
dfHNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = 1/RL.ranking)%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = -1/RL.ranking) %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

dfscores <- 
  full_join(dfHNA, dfLNA, by = "OTU") %>%
  rename(HNA = RL.ranking.HNA, LNA = RL.ranking.LNA) %>%
  gather(key = FCM_type, value = RL.ranking, HNA:LNA)
    

ggplot(dfscores, aes(FCM_type, OTU)) + 
  geom_tile(aes(fill = RL.ranking),colour = "white") + 
  scale_fill_gradient2(low = "darkgoldenrod1", high = "deepskyblue4", mid = "white", 
                       midpoint = 0, na.value = "black") +
  theme(axis.title.x = element_blank())
ggsave("heatmap_figs/muskegon_heat_0.15.jpg", width = 4, height = 10)


#######
####
#
#
#
#
#
# Read in Data
# Muskegon
HNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "HNA", Lake = "Muskegon",
         RL.ranking = 1/RL.ranking)

LNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "LNA", Lake = "Muskegon",
         RL.ranking = -1/RL.ranking)

# Inland
HNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "HNA", Lake = "Inland",
         RL.ranking = 1/RL.ranking)

LNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "LNA", Lake = "Inland",
         RL.ranking = -1/RL.ranking)


# Michigan
HNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "HNA", Lake = "Michigan",
         RL.ranking = 1/RL.ranking)

LNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(FCM_type = "LNA", Lake = "Michigan",
         RL.ranking = -1/RL.ranking)


lakes_scores_df <- bind_rows(HNA_musk, LNA_musk, HNA_inland, LNA_inland, HNA_michigan, LNA_michigan) %>%
  dplyr::select(OTU, RL.ranking, FCM_type, Lake)


ggplot(lakes_scores_df, aes(y=RL.ranking, x=OTU, fill=FCM_type)) + 
  geom_bar(stat="identity", position="identity") + coord_flip() + #ggtitle("Summed OTUs") +
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  theme_classic() + 
  facet_grid(. ~ Lake) +
  ylab("Inverse RL Ranking") +
  scale_fill_manual(values = fcm_colors) +
  theme(legend.title = element_blank(), legend.position = c(0.65, 0.98))
ggsave("heatmap_figs/all_lakes_bar_0.15_byLake.jpg", width = 8, height = 25)






### HEATMAP ATTEMPT
#### Attempt 3
# Read in Data
dfHNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = 1/RL.ranking, Lake = "Muskegon")%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = -1/RL.ranking, Lake = "Muskegon") %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

dfHNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = 1/RL.ranking, Lake = "Inland")%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = -1/RL.ranking, Lake = "Inland") %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

dfHNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = 1/RL.ranking, Lake = "Michigan")%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.15) %>%
  mutate(RL.ranking = -1/RL.ranking, Lake = "Michigan") %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

muskegon_df <- full_join(dfHNA_musk, dfLNA_musk, by = "OTU") %>%
  mutate(Lake = "Muskegon")
inland_df <- full_join(dfHNA_inland, dfLNA_inland, by = "OTU") %>%
  mutate(Lake = "Inland")
michigan_df <- full_join(dfHNA_michigan, dfLNA_michigan, by = "OTU") %>%
  mutate(Lake = "Michigan")

lake_dfscores <-
  bind_rows(muskegon_df, inland_df, michigan_df) %>%
  rename(HNA = RL.ranking.HNA, LNA = RL.ranking.LNA) %>%
  gather(key = FCM_type, value = RL.ranking, HNA:LNA)


ggplot(lake_dfscores, aes(FCM_type, OTU)) + 
  geom_tile(aes(fill = RL.ranking),colour = "white") + 
  scale_fill_gradient2(low = "darkgoldenrod1", high = "deepskyblue4", mid = "white", 
                       midpoint = 0, na.value = "black") +
  facet_grid(~Lake) +
  theme(axis.title.x = element_blank())
ggsave("heatmap_figs/all_lakes_heat_0.15_byLake.jpg", width = 8, height = 25)


ggplot(lake_dfscores, aes(Lake, OTU)) + 
  geom_tile(aes(fill = RL.ranking),colour = "white") + 
  scale_fill_gradient2(low = "darkgoldenrod1", high = "deepskyblue4", mid = "white", 
                       midpoint = 0, na.value = "black") +
  facet_grid(~FCM_type) +
  theme(axis.title.x = element_blank())
ggsave("heatmap_figs/all_lakes_heat_0.15_byFCM_type.jpg", width = 8, height = 25)



### CLUSTERING
library(gplots)
library(tibble)
library(d3heatmap)

?heatmap.2

# Reformat data to be matrix style
musk_dat <- 
  muskegon_df %>%
  dplyr::select(-Lake) %>%
  dplyr::filter(RL.ranking.HNA > 0.15 | RL.ranking.LNA < -0.15) %>%
  rename(musk_HNA = RL.ranking.HNA, 
         musk_LNA = RL.ranking.LNA)%>%
  mutate(musk_LNA = musk_LNA * -1) 


inland_dat <- 
  inland_df %>%
  dplyr::select(-Lake) %>%
  dplyr::filter(RL.ranking.HNA > 0.15 | RL.ranking.LNA < -0.15) %>%
  rename(inland_HNA = RL.ranking.HNA, 
         inland_LNA = RL.ranking.LNA) %>%
  mutate(inland_LNA = inland_LNA * -1)

michigan_dat <- 
  michigan_df %>%
  dplyr::select(-Lake) %>%
  dplyr::filter(RL.ranking.HNA > 0.15 | RL.ranking.LNA < -0.15) %>%
  rename(michigan_HNA = RL.ranking.HNA, 
         michigan_LNA = RL.ranking.LNA) %>%
  mutate(michigan_LNA = michigan_LNA * -1)

matrix_scores <- 
  full_join(musk_dat, inland_dat, by = "OTU") %>%
  full_join(michigan_dat, by = "OTU") %>%
  tibble::column_to_rownames(var = "OTU") %>%
  as.matrix()

matrix_scores[is.na(matrix_scores)] <- 0

breakers <- seq(min(matrix_scores, na.rm = T), max(matrix_scores, na.rm = T), length.out = 21)
my_palette <- colorRampPalette(c("grey","black", "red")) (n=20)

col <- c("deepskyblue4", "darkgoldenrod1", "deepskyblue4", 
         "darkgoldenrod1", "deepskyblue4", "darkgoldenrod1")

colz <- c("pink", "pink", "green", "green", "orange", "orange")

jpeg(file = "heatmap_figs/clustering_heatmap_all.jpg", 
     units = "in", width = 8, height = 10, res = 300)
heatmap.2(matrix_scores, col = my_palette, breaks=breakers, trace="none",
          key=TRUE, symkey=FALSE, density.info="none", srtCol=45,
          margins=c(8,7), cexRow=1.25,cexCol=1.25,
          ColSideColors=colz)
dev.off()  # 318 bytes file in current directory



# Separate plots for HNA and LNA
HNA_matrix <- 
  full_join(musk_dat, inland_dat, by = "OTU") %>%
  full_join(michigan_dat, by = "OTU") %>%
  dplyr::select(OTU, musk_HNA, inland_HNA, michigan_HNA) %>%
  dplyr::filter(musk_HNA > 0.15 | inland_HNA > 0.15 | michigan_HNA > 0.15) %>%
  rename(Muskegon = musk_HNA, 
         Inland = inland_HNA,
         Michigan = michigan_HNA) %>%
  tibble::column_to_rownames(var = "OTU") %>%
  as.matrix()
  
HNA_matrix[is.na(HNA_matrix)] <- 0

hna_breaks <- seq(min(HNA_matrix, na.rm = T), max(HNA_matrix, na.rm = T), length.out = 11)
hna_palette <- colorRampPalette(c("grey","black", "deepskyblue4")) (n=10)


jpeg(file = "heatmap_figs/clustering_HNA.jpg", 
     units = "in", width = 7, height = 6.5, res = 300)
heatmap.2(HNA_matrix, col = hna_palette, breaks=hna_breaks, trace="none",
          key=TRUE, symkey=FALSE, density.info="none", srtCol=0,
          margins=c(3,7), cexRow=1.25,cexCol=1.25, key.xlab="Lasso Score", key.title = NA)
dev.off() 

# Separate plots for HNA and LNA
LNA_matrix <- 
  full_join(musk_dat, inland_dat, by = "OTU") %>%
  full_join(michigan_dat, by = "OTU") %>%
  dplyr::select(OTU, musk_LNA, inland_LNA, michigan_LNA) %>%
  dplyr::filter(musk_LNA > 0.15 | inland_LNA > 0.15 | michigan_LNA > 0.15) %>%
  rename(Muskegon = musk_LNA, 
         Inland = inland_LNA,
         Michigan = michigan_LNA) %>%
  tibble::column_to_rownames(var = "OTU") %>%
  as.matrix()

LNA_matrix[is.na(LNA_matrix)] <- 0

lna_breaks <- seq(min(LNA_matrix, na.rm = T), max(LNA_matrix, na.rm = T), length.out = 11)
lna_palette <- colorRampPalette(c("grey","black", "darkgoldenrod1")) (n=10)

jpeg(file = "heatmap_figs/clustering_LNA.jpg", 
     units = "in", width = 7, height = 6.5, res = 300)
heatmap.2(LNA_matrix, col = lna_palette, breaks=lna_breaks, trace="none",
          key=TRUE, symkey=FALSE, density.info="none", srtCol=0,
          margins=c(3,7), cexRow=1.25,cexCol=1.25, key.xlab="Lasso Score", key.title = NA)
dev.off() 


# Taxonomic analysis 
load("data/Chloroplasts_removed/ByLake_Filtering/5in10/muskegon/muskegon_5in10_physeqs.RData")

muskegonOTUs <- 
  matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, musk_HNA, musk_LNA) %>%
  dplyr::filter(musk_HNA > 0.15 | musk_LNA > 0.15) 

list_muskegon_otus <- muskegonOTUs$OTU

muskegon_taxonomy <- rare_muskegon_physeq_5in10_abs %>%
  subset_taxa(., taxa_names(rare_muskegon_physeq_5in10_abs) %in% list_muskegon_otus)

as.data.frame(tax_table(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(muskegonOTUs, by = "OTU") %>%
  rename(HNA = musk_HNA, LNA = musk_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.15) %>%
  ggplot(aes(y = RL_score, x = Rank5, fill = Rank5)) +
  geom_bar(stat = "identity", position = "stack", color = "black") + facet_wrap(~fcm_group) +
  scale_fill_brewer(palette="Set3") + xlab("Family") +
  ggtitle("Muskegon Lake") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))

as.data.frame(sample_data(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU")




load("data/Chloroplasts_removed/ByLake_Filtering/5in10/inland/inland_5in10_physeqs.RData")
rare_inland_physeq_5in10_abs

inlandOTUs <- 
  matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, inland_HNA, inland_LNA) %>%
  dplyr::filter(inland_HNA > 0.15 | inland_LNA > 0.15) 

list_inland_otus <- inlandOTUs$OTU

inland_taxonomy <- rare_inland_physeq_5in10_abs %>%
  subset_taxa(., taxa_names(rare_inland_physeq_5in10_abs) %in% list_inland_otus)

as.data.frame(tax_table(inland_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(inlandOTUs, by = "OTU") %>%
  rename(HNA = inland_HNA, LNA = inland_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.15) %>%
  ggplot(aes(y = RL_score, x = Rank4, fill = Rank4)) +
  geom_bar(stat = "identity") + facet_wrap(~fcm_group, scale = "free_x") +
  scale_fill_brewer(palette="Set1") + xlab("Order") +
  ggtitle("Inland Lakes") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))



load("data/Chloroplasts_removed/ByLake_Filtering/5in10/michigan/michigan_5in10_physeqs.RData")
rare_michigan_physeq_5in10_abs

michiganOTUs <- 
  matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, michigan_HNA, michigan_LNA) %>%
  dplyr::filter(michigan_HNA > 0.15 | michigan_LNA > 0.15) 

list_michigan_otus <- michiganOTUs$OTU

michigan_taxonomy <- rare_michigan_physeq_5in10_abs %>%
  subset_taxa(., taxa_names(rare_michigan_physeq_5in10_abs) %in% list_michigan_otus)

as.data.frame(tax_table(michigan_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(michiganOTUs, by = "OTU") %>%
  rename(HNA = michigan_HNA, LNA = michigan_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.15) %>%
  ggplot(aes(y = RL_score, x = Rank4, fill = Rank4)) +
  geom_bar(stat = "identity", color = "black") + facet_wrap(~fcm_group, scale = "free_x") +
  scale_fill_brewer(palette="Set3") + xlab("Order") +
  ggtitle("Lake Michigan") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))


##### THE TOP 6
ones <- c("Otu000173", "Otu000029", "Otu000369", "Otu000555", "Otu000025", "Otu000168")


top_inland_otu_tax <- 
  as.data.frame(tax_table(inland_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(inlandOTUs, by = "OTU") %>%
  dplyr::filter(OTU %in% ones) %>%
  rename(HNA = inland_HNA, LNA = inland_LNA)

top_musk_otu_tax <- 
  as.data.frame(tax_table(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(muskegonOTUs, by = "OTU") %>%
  dplyr::filter(OTU %in% ones) %>%
  rename(HNA = musk_HNA, LNA = musk_LNA)

top_michigan_otu_tax <- 
  as.data.frame(tax_table(michigan_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(michiganOTUs, by = "OTU") %>%
  dplyr::filter(OTU %in% ones) %>%
  rename(HNA = michigan_HNA, LNA = michigan_LNA)

df <- bind_rows(top_inland_otu_tax, top_musk_otu_tax, top_michigan_otu_tax)




#### PHYLOGENETIC ANALYSIS
load("data/fasttree/ALL-physeq-for-phylo.RData")
ALL_rel_physeq_5in10

otu_scores_df <- matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU")

vector_of_otus <- as.vector(otu_scores_df$OTU)
write(vector_of_otus, file = "heatmap_figs/OTUnames_based_on_RLscores.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")

### Make a phyloseq object 
otu_scores_df <- matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "OTU") 

library(phytools)
HNALNA_otu_tree <- read.newick(file="data/fasttree/newick_tree_HNALNA_rmN.tre")

## Prepare the tax table 
michigan_otu_tax <- 
  as.data.frame(tax_table(michigan_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(michiganOTUs, by = "OTU") %>%
  rename(HNA = michigan_HNA, LNA = michigan_LNA) %>%
  filter(OTU != "Otu000242")%>%
  dplyr::select(-c(HNA:LNA))

inland_otu_tax <- 
  as.data.frame(tax_table(inland_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(inlandOTUs, by = "OTU") %>%
  rename(HNA = inland_HNA, LNA = inland_LNA)%>%
  dplyr::select(-c(HNA:LNA))

muskegon_otu_tax <- 
  as.data.frame(tax_table(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(muskegonOTUs, by = "OTU") %>%
  rename(HNA = musk_HNA, LNA = musk_LNA) %>%
  dplyr::select(-c(HNA:LNA))

taxtable <- bind_rows(michigan_otu_tax, inland_otu_tax, muskegon_otu_tax) %>%
  rename(Domain = Rank1, Phylum = Rank2, Class = Rank3, 
         Order = Rank4, Family = Rank5, Genus = Rank6, Species = Rank7) %>%
  tibble::column_to_rownames(var = "OTU") 


merge_phyloseq(HNALNA_otu_tree, otu_table(otu_scores_df, taxa_are_rows = TRUE), tax_table(taxtable))


phylo_otus <- ALL_rel_physeq_5in10 %>%
  subset_taxa(., taxa_names(ALL_rel_physeq_5in10) %in% otu_scores_df$OTU)  

scratch_otu <- otu_table(phylo_otus)

PHYTREE <- merge_phyloseq(phylo_otus, phy_tree(HNALNA_otu_tree))



plot_tree(phylo_otus, color="Rank2", ladderize="left") 



