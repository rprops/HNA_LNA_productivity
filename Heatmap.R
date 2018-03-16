# Figure 4
# Feb 28, 2018

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gplots)
library(tibble)
library(d3heatmap)

# Set HNA and LNA discrete colors
fcm_colors <- c(
  "HNA" = "deepskyblue4",
  "LNA" = "darkgoldenrod1",
  "Total" = "black")

# Optimal Thresholds from March 14th
# HNA/LNA: Inland: 0.13/0.108,  Michigan: 0.248/0.286,  Muskegon: 0.09/0.09

### HEATMAP 
# Read in Data
dfHNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.09) %>%
  mutate(RL.ranking = 1/RL.ranking, Lake = "Muskegon")%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.09) %>%
  mutate(RL.ranking = -1/RL.ranking, Lake = "Muskegon") %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

dfHNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.13) %>%
  mutate(RL.ranking = 1/RL.ranking, Lake = "Inland")%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.108) %>%
  mutate(RL.ranking = -1/RL.ranking, Lake = "Inland") %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

dfHNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.248) %>%
  mutate(RL.ranking = 1/RL.ranking, Lake = "Michigan")%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.286) %>%
  mutate(RL.ranking = -1/RL.ranking, Lake = "Michigan") %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

# Combine HNA and LNA datasets together
muskegon_df <- full_join(dfHNA_musk, dfLNA_musk, by = "OTU") %>%   mutate(Lake = "Muskegon")
inland_df <- full_join(dfHNA_inland, dfLNA_inland, by = "OTU") %>%  mutate(Lake = "Inland")
michigan_df <- full_join(dfHNA_michigan, dfLNA_michigan, by = "OTU") %>%  mutate(Lake = "Michigan")

# Combine all of the three lakes together into one dataframe! 
lake_dfscores <-
  bind_rows(muskegon_df, inland_df, michigan_df) %>%
  rename(HNA = RL.ranking.HNA, LNA = RL.ranking.LNA) %>%
  gather(key = FCM_type, value = RL.ranking, HNA:LNA)



# Reformat data to be matrix style for heatmapping
# HNA/LNA: Inland: 0.13/0.108,  Michigan: 0.248/0.286,  Muskegon: 0.09/0.09
musk_dat <- 
  muskegon_df %>%
  dplyr::select(-Lake) %>%
  dplyr::filter(RL.ranking.HNA > 0.09 | RL.ranking.LNA < -0.09) %>%
  rename(Muskegon_HNA = RL.ranking.HNA, 
         Muskegon_LNA = RL.ranking.LNA)%>%
  mutate(Muskegon_LNA = Muskegon_LNA * -1) 

inland_dat <- 
  inland_df %>%
  dplyr::select(-Lake) %>%
  dplyr::filter(RL.ranking.HNA > 0.13 | RL.ranking.LNA < -0.108) %>%
  rename(Inland_HNA = RL.ranking.HNA, 
         Inland_LNA = RL.ranking.LNA) %>%
  mutate(Inland_LNA = Inland_LNA * -1)

michigan_dat <- 
  michigan_df %>%
  dplyr::select(-Lake) %>%
  dplyr::filter(RL.ranking.HNA > 0.248 | RL.ranking.LNA < -0.286) %>%
  rename(Michigan_HNA = RL.ranking.HNA, 
         Michigan_LNA = RL.ranking.LNA) %>%
  mutate(Michigan_LNA = Michigan_LNA * -1)

matrix_scores <- 
  full_join(musk_dat, inland_dat, by = "OTU") %>%
  full_join(michigan_dat, by = "OTU") %>%
  tibble::column_to_rownames(var = "OTU") %>%
  as.matrix()

matrix_scores[is.na(matrix_scores)] <- 0

breakers <- seq(min(matrix_scores, na.rm = T), max(matrix_scores, na.rm = T), length.out = 21)
my_palette <- colorRampPalette(c("white","gray65", "red")) (n=20)

colz <- c("pink", "pink", "green", "green", "orange", "orange")

jpeg(file = "heatmap_figs/clustering_heatmap_all.jpg", 
     units = "in", width = 8, height = 10, res = 300)
heatmap.2(matrix_scores, col = my_palette, breaks=breakers, trace="none",
          key=TRUE, symkey=FALSE, density.info="none", srtCol=30,
          margins=c(5.5,7), cexRow=1.25,cexCol=1.25,
          ColSideColors=colz)
dev.off() 


# Separate plots for HNA and LNA 
#   HNA/LNA --> Inland: 0.13/0.108,  Michigan: 0.248/0.286,  Muskegon: 0.09/0.09
HNA_matrix <- 
  full_join(musk_dat, inland_dat, by = "OTU") %>%
  full_join(michigan_dat, by = "OTU") %>%
  dplyr::select(OTU, Muskegon_HNA, Inland_HNA, Michigan_HNA) %>%
  dplyr::filter(Muskegon_HNA > 0.09 | Inland_HNA > 0.13 | Michigan_HNA > 0.248) %>%
  rename(Muskegon = Muskegon_HNA, 
         Inland = Inland_HNA,
         Michigan = Michigan_HNA) %>%
  tibble::column_to_rownames(var = "OTU") %>%
  as.matrix()
  
HNA_matrix[is.na(HNA_matrix)] <- 0

hna_breaks <- seq(min(HNA_matrix, na.rm = T), max(HNA_matrix, na.rm = T), length.out = 13)
hna_palette <- colorRampPalette(c("white","grey65", "deepskyblue4")) (n=12)

# Plot only the HNA
jpeg(file = "heatmap_figs/clustering_HNA.jpg", 
     units = "in", width = 7, height = 6.5, res = 300)
heatmap.2(HNA_matrix, col = hna_palette, breaks=hna_breaks, trace="none",
          key=TRUE, symkey=FALSE, density.info="none", srtCol=30,
          margins=c(4,7), cexRow=1.25,cexCol=1.25, key.xlab="Lasso Score", key.title = NA)
dev.off() 

# Separate plots for HNA and LNA
#   HNA/LNA --> Inland: 0.13/0.108,  Michigan: 0.248/0.286,  Muskegon: 0.09/0.09
LNA_matrix <- 
  full_join(musk_dat, inland_dat, by = "OTU") %>%
  full_join(michigan_dat, by = "OTU") %>%
  dplyr::select(OTU, Muskegon_LNA, Inland_LNA, Michigan_LNA) %>%
  dplyr::filter(Muskegon_LNA > 0.09 | Inland_LNA > 0.108 | Michigan_LNA > 0.286) %>%
  rename(Muskegon = Muskegon_LNA, 
         Inland = Inland_LNA,
         Michigan = Michigan_LNA) %>%
  tibble::column_to_rownames(var = "OTU") %>%
  as.matrix()

LNA_matrix[is.na(LNA_matrix)] <- 0

lna_breaks <- seq(min(LNA_matrix, na.rm = T), max(LNA_matrix, na.rm = T), length.out = 13)
lna_palette <- colorRampPalette(c("white","grey65", "darkgoldenrod1")) (n=12)

# Plot only the LNA
jpeg(file = "heatmap_figs/clustering_LNA.jpg", 
     units = "in", width = 7, height = 6.5, res = 300)
heatmap.2(LNA_matrix, col = lna_palette, breaks=lna_breaks, trace="none",
          key=TRUE, symkey=FALSE, density.info="none", srtCol=30,
          margins=c(4,7), cexRow=1.25,cexCol=1.25, key.xlab="Lasso Score", key.title = NA)
dev.off() 


# Taxonomic analysis
#   HNA/LNA --> Inland: 0.13/0.108,  Michigan: 0.248/0.286,  Muskegon: 0.09/0.09
load("data/Chloroplasts_removed/ByLake_Filtering/5in10/muskegon/muskegon_5in10_physeqs.RData")

muskegonOTUs <- 
  matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, Muskegon_HNA, Muskegon_LNA) %>%
  dplyr::filter(Muskegon_HNA > 0.09 | Muskegon_LNA > 0.09) 

list_muskegon_otus <- muskegonOTUs$OTU
muskegon_taxonomy <- rare_muskegon_physeq_5in10_abs %>%
  subset_taxa(., taxa_names(rare_muskegon_physeq_5in10_abs) %in% list_muskegon_otus)

as.data.frame(tax_table(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(muskegonOTUs, by = "OTU") %>%
  rename(HNA = Muskegon_HNA, LNA = Muskegon_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.09) %>%
  ggplot(aes(y = RL_score, x = Rank5, fill = Rank2)) + xlab("Family") +
  geom_bar(stat = "identity", position = "stack",color = "black") + 
  facet_wrap(~fcm_group, scale = "free") +  ggtitle("Muskegon Lake") +  
  scale_fill_brewer(palette="Set3") + xlab("Family") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.3)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))

as.data.frame(sample_data(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU")


#### Load Inland lake data!
load("data/Chloroplasts_removed/ByLake_Filtering/5in10/inland/inland_5in10_physeqs.RData")
rare_inland_physeq_5in10_abs
#   HNA/LNA --> Inland: 0.13/0.108,  Michigan: 0.248/0.286,  Muskegon: 0.09/0.09

inlandOTUs <- 
  matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, Inland_HNA, Inland_LNA) %>%
  dplyr::filter(Inland_HNA > 0.13 | Inland_LNA > 0.108) 

list_inland_otus <- inlandOTUs$OTU

inland_taxonomy <- rare_inland_physeq_5in10_abs %>%
  subset_taxa(., taxa_names(rare_inland_physeq_5in10_abs) %in% list_inland_otus)

as.data.frame(tax_table(inland_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(inlandOTUs, by = "OTU") %>%
  dplyr::filter(Inland_HNA > 0.13 | Inland_LNA > 0.108) %>%
  rename(HNA = Inland_HNA, LNA = Inland_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.15) %>%
  ggplot(aes(y = RL_score, x = Rank4, fill = Rank2)) +
  geom_bar(stat = "identity") + facet_wrap(~fcm_group, scale = "free_x") +
  scale_fill_brewer(palette="Set1") + xlab("Order") +
  ggtitle("Inland Lakes") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))

##### Load Lake Michigan Data!
load("data/Chloroplasts_removed/ByLake_Filtering/5in10/michigan/michigan_5in10_physeqs.RData")
rare_michigan_physeq_5in10_abs
#   HNA/LNA --> Inland: 0.13/0.108,  Michigan: 0.248/0.286,  Muskegon: 0.09/0.09

michiganOTUs <- 
  matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, Michigan_HNA, Michigan_LNA) %>%
  dplyr::filter(Michigan_HNA > 0.248 | Michigan_LNA > 0.286) 

list_michigan_otus <- michiganOTUs$OTU

michigan_taxonomy <- rare_michigan_physeq_5in10_abs %>%
  subset_taxa(., taxa_names(rare_michigan_physeq_5in10_abs) %in% list_michigan_otus)

as.data.frame(tax_table(michigan_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(michiganOTUs, by = "OTU") %>%
  dplyr::filter(Michigan_HNA > 0.248 | Michigan_LNA > 0.286) %>%
  rename(HNA = Michigan_HNA, LNA = Michigan_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.15) %>%
  ggplot(aes(y = RL_score, x = Rank5, fill = Rank2)) +
  geom_bar(stat = "identity", color = "black") + facet_wrap(~fcm_group, scale = "free_x") +
  scale_fill_brewer(palette="Set3") + xlab("Family") +
  ggtitle("Lake Michigan") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))





#### PHYLOGENETIC ANALYSIS
load("/data/phyloseq.RData")
physeq.otu

otu_scores_df <- matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("OTU")

vector_of_otus <- as.vector(otu_scores_df$OTU)

physeq <- physeq.otu %>%
  subset_taxa(., taxa_names(physeq.otu) %in% vector_of_otus)  

# Which OTU is missing?
setdiff(sort(vector_of_otus), sort(taxa_names(physeq)))
setdiff(sort(taxa_names(physeq)), sort(vector_of_otus))

write(vector_of_otus, file = "data/fasttree/OTUnames_based_on_RLscores.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")

### Make a phyloseq object 
otu_scores_df <- matrix_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "OTU") 

library(phytools)
HNALNA_otu_tree <- read.tree(file="data/fasttree/RAxML_bipartitions.newick_tree_HNALNA_rmN.tre")

tree_tip_order <- data.frame(HNALNA_otu_tree$tip.label) %>%
  rename(OTU = HNALNA_otu_tree.tip.label)


final_physeq <- merge_phyloseq(physeq, phy_tree(HNALNA_otu_tree))

# Fix the taxonomy names 
colnames(tax_table(final_physeq)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

###################################################################### ADD THE PROTEOBACTERIA TO THE PHYLA
phy <- data.frame(tax_table(final_physeq))
Phylum <- as.character(phy$Phylum)
Class <- as.character(phy$Class)

for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    Phylum[i] <- Class[i]
  } 
}

phy$Phylum <- Phylum # Add the new phylum level data back to phy
#phy$OTU <- row.names(tax_table(final_physeq)) # Make a column for OTU

t <- tax_table(as.matrix(phy))

tax_table(final_physeq) <- t

# Set global colors for the different taxonomic phyla
phylum_colors <- c( 
  Acidobacteria = "navy", 
  Actinobacteria = "blue", 
  Alphaproteobacteria = "tomato", 
  Aminicenantes = "cornflowerblue",
  Armatimonadetes = "wheat", 
  Bacteria_unclassified = "#508578", 
  Bacteroidetes = "gold", 
  Betaproteobacteria = "plum1", 
  "Candidate_division_OP3" = "slategray3",
  Chlamydiae = "#A20E42",
  Chlorobi="magenta", 
  Chloroflexi="black", 
  Cyanobacteria = "limegreen",
  "Deinococcus-Thermus" = "black",
  Deltaproteobacteria = "olivedrab", 
  Firmicutes = "black",
  Gammaproteobacteria = "cyan",
  Gemmatimonadetes = "yellow",
  Gracilibacteria = "#FD823F",
  JTB23 = "#B5D6AA",
  Latescibacteria = "salmon4",
  Lentisphaerae = "palevioletred1",
  Nitrospirae = "forestgreen",
  Omnitrophica = "violet",
  Parcubacteria = "#531A4D",
  Planctomycetes = "darkorange", 
  Proteobacteria_unclassified = "greenyellow",
  Spirochaetae = "royalblue",
  TA06 = "peachpuff",
  Omnitrophica = "burlywood", 
  unknown_unclassified = "grey",
  Verrucomicrobia = "purple4",
  Proteobacteria_unclassified = "green",
  Proteobacteria = "red")


plot_tree(final_physeq, "treeonly", nodelabf = nodeplotboot(), ladderize = "left", 
          label.tips="taxa_names")
  

## Test on a smaller tree
test <- final_physeq %>%
  subset_taxa(., taxa_names(final_physeq) %in% vector_of_otus[1:10])  


plot_tree(test, "treeonly", nodelabf = nodeplotboot(), ladderize = "left", 
          label.tips="taxa_names", color = "Phylum") 


library(ggtree)

tax <- data.frame(tax_table(final_physeq)) %>%
  tibble::rownames_to_column(var = "OTU")

df <- read.csv("data/fasttree/OTUnames_based_on_RLscores_MANUAL.csv") %>%
  left_join(., tax, by = "OTU") %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::select(Lake, fcm_type)

p <- ggplot(HNALNA_otu_tree, aes(x, y)) + geom_tree() + theme_tree() +
  geom_tiplab(size=3, align=TRUE, linesize=.5) 
gheatmap(p, df, offset = 0.5, width=0.5, font.size=3, colnames_angle=0, hjust=0.5) +
  scale_fill_brewer(palette = "Set1")
ggsave("heatmap_figs/phylogenetic_tree_fcm_lake.jpg", width = 8, height = 8)


df2 <- read.csv("data/fasttree/OTUnames_based_on_RLscores_MANUAL.csv") %>%
  left_join(., tax, by = "OTU") %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::select(Phylum)

p2 <- ggplot(HNALNA_otu_tree, aes(x, y)) + geom_tree() + theme_tree() +
  geom_tiplab(size=3, align=TRUE, linesize=.5) 
gheatmap(p2, df2, offset = 0.5, width=0.5, font.size=3, colnames_angle=0, hjust=0.5) +
  scale_fill_manual(values = phylum_colors)
ggsave("heatmap_figs/phylogenetic_tree_phylum.jpg", width = 8, height = 8)






############################################################
############################################################
######## EXTRA CODE ########################################
######## OLDER CODE FOR geom_tile and geom_bar plots #######
############################################################
############################################################

# Plot by Lake
ggplot(lake_dfscores, aes(FCM_type, OTU)) + 
  geom_tile(aes(fill = RL.ranking),colour = "white") + 
  scale_fill_gradient2(low = "darkgoldenrod1", high = "deepskyblue4", mid = "white", 
                       midpoint = 0, na.value = "black") +
  facet_grid(~Lake) + theme(axis.title.x = element_blank())
ggsave("heatmap_figs/all_lakes_heat_optimalthreshold_byLake.jpg", width = 8, height = 25)

# Plot by FCM Type
ggplot(lake_dfscores, aes(Lake, OTU)) + 
  geom_tile(aes(fill = RL.ranking),colour = "white") + 
  scale_fill_gradient2(low = "darkgoldenrod1", high = "deepskyblue4", mid = "white", 
                       midpoint = 0, na.value = "black") +
  facet_grid(~FCM_type) + theme(axis.title.x = element_blank())
ggsave("heatmap_figs/all_lakes_heat_optimalthreshold_byFCM_type.jpg", width = 8, height = 25)



# Read in Data
#   Thresholds for HNA/LNA Muskegon: 0.09/0.09
HNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.09) %>%
  mutate(FCM_type = "HNA", 
         RL.ranking = 1/RL.ranking)

LNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.09) %>%
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
#   Thresholds for HNA/LNA Muskegon: 0.09/0.09
# Read in Data
dfHNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.09) %>%
  mutate(RL.ranking = 1/RL.ranking)%>%
  rename(OTU = X, RL.ranking.HNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.HNA)

dfLNA <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  dplyr::filter(RL.score > 0.09) %>%
  mutate(RL.ranking = -1/RL.ranking) %>%
  rename(OTU = X, RL.ranking.LNA = RL.ranking) %>%
  dplyr::select(OTU, RL.ranking.LNA)

dfscores <- 
  full_join(dfHNA, dfLNA, by = "OTU") %>%
  rename(HNA = RL.ranking.HNA, LNA = RL.ranking.LNA) %>%
  gather(key = FCM_type, value = RL.ranking, HNA:LNA)

# Bar Plot 
ggplot(dfscores, aes(y=RL.ranking, x=OTU, fill=FCM_type)) + 
  geom_bar(stat="identity", position="identity") + coord_flip() + #ggtitle("Summed OTUs") +
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  theme_classic() +
  ylab("Inverse RL Ranking") +
  scale_fill_manual(values = fcm_colors) +
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.95))
ggsave("heatmap_figs/muskegon_bar_0.09.jpg", width = 4, height = 8)



#### Attempt 3
ggplot(dfscores, aes(FCM_type, OTU)) + 
  geom_tile(aes(fill = RL.ranking),colour = "white") + 
  scale_fill_gradient2(low = "darkgoldenrod1", high = "deepskyblue4", mid = "white", 
                       midpoint = 0, na.value = "black") +
  theme(axis.title.x = element_blank())
ggsave("heatmap_figs/muskegon_heat_0.09.jpg", width = 4, height = 10)


######## Read in Data for all of the lakes 
# Optimal thresholds for HNA/LNA: 
#   Inland: 0.13/0.108
#   Michigan: 0.248/0.286
#   Muskegon: 0.09/0.09

HNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.09) %>%
  mutate(FCM_type = "HNA", Lake = "Muskegon",
         RL.ranking = 1/RL.ranking)
LNA_musk <- read.csv("Final/FS_Scores/Muskegon_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.09) %>%
  mutate(FCM_type = "LNA", Lake = "Muskegon",
         RL.ranking = -1/RL.ranking)

# Inland
HNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.13) %>%
  mutate(FCM_type = "HNA", Lake = "Inland",
         RL.ranking = 1/RL.ranking)
LNA_inland <- read.csv("Final/FS_Scores/Inland_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.108) %>%
  mutate(FCM_type = "LNA", Lake = "Inland",
         RL.ranking = -1/RL.ranking)

# Michigan
HNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_HNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.248) %>%
  mutate(FCM_type = "HNA", Lake = "Michigan",
         RL.ranking = 1/RL.ranking)
LNA_michigan <- read.csv("Final/FS_Scores/Michigan_fs_scores_LNA_5seq10.csv") %>%
  rename(OTU = X) %>%
  dplyr::filter(RL.score > 0.286) %>%
  mutate(FCM_type = "LNA", Lake = "Michigan",
         RL.ranking = -1/RL.ranking)

# Combine into one dataframe
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
ggsave("heatmap_figs/all_lakes_bar_optimalthresholds_byLake.jpg", width = 8, height = 25)





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






####
## Prepare the tax table 
michigan_otu_tax <- 
  as.data.frame(tax_table(michigan_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(michiganOTUs, by = "OTU") %>%
  rename(HNA = Michigan_HNA, LNA = Michigan_LNA) %>%
  filter(OTU != "Otu000242")%>%
  dplyr::select(-c(HNA:LNA))

inland_otu_tax <- 
  as.data.frame(tax_table(inland_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(inlandOTUs, by = "OTU") %>%
  rename(HNA = Inland_HNA, LNA = Inland_LNA)%>%
  dplyr::select(-c(HNA:LNA))

muskegon_otu_tax <- 
  as.data.frame(tax_table(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(muskegonOTUs, by = "OTU") %>%
  rename(HNA = Muskegon_HNA, LNA = Muskegon_LNA) %>%
  dplyr::select(-c(HNA:LNA))

taxtable <- bind_rows(michigan_otu_tax, inland_otu_tax, muskegon_otu_tax) %>%
  rename(Domain = Rank1, Phylum = Rank2, Class = Rank3, 
         Order = Rank4, Family = Rank5, Genus = Rank6, Species = Rank7)

ordered_tax_table <- left_join(tree_tip_order, taxtable, by = "OTU") %>%
  tibble::column_to_rownames(var = "OTU")

ordered_tax_table$OTU <- row.names(ordered_tax_table)

