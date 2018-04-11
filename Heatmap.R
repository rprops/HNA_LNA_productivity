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
library(phytools)
library(ggtree)

# Set HNA and LNA discrete colors
fcm_colors <- c(
  "HNA" = "deepskyblue4",
  "LNA" = "darkgoldenrod1",
  "Total" = "black")

# Set global colors for the different taxonomic phyla
phylum_colors <- c( 
  Acidobacteria = "navy", 
  Actinobacteria = "blue", 
  Alphaproteobacteria = "orangered", 
  Aminicenantes = "",
  Armatimonadetes = "wheat", 
  Bacteria_unclassified = "grey47", 
  Bacteroidetes = "cornflowerblue", 
  Betaproteobacteria = "plum1", 
  "Candidate_division_OP3" = "slategray3",
  Chlamydiae = "#A20E42",
  Chlorobi="magenta", 
  Chloroflexi="black", 
  Cyanobacteria = "limegreen",
  "Deinococcus-Thermus" = "black",
  Deltaproteobacteria = "olivedrab", 
  Firmicutes = "navy",
  Gammaproteobacteria = "cyan",
  Gemmatimonadetes = "yellow",
  Gracilibacteria = "#FD823F",
  JTB23 = "#B5D6AA",
  Latescibacteria = "salmon4",
  Lentisphaerae = "palevioletred1",
  Nitrospirae = "forestgreen",
  Omnitrophica = "red4",
  Parcubacteria = "#531A4D",
  Planctomycetes = "darkorange", 
  Proteobacteria_unclassified = "greenyellow",
  Spirochaetae = "royalblue",
  TA06 = "peachpuff",
  Omnitrophica = "burlywood", 
  unknown_unclassified = "grey88",
  Verrucomicrobia = "purple",
  Proteobacteria_unclassified = "green",
  Proteobacteria = "red", 
  HNA =  "deepskyblue4",
  LNA = "darkgoldenrod1",
  Both = "black" )

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
heatmap.2(matrix_scores, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x,method = "complete"),
          col = my_palette, breaks=breakers, trace="none",
          key=TRUE, symkey=FALSE, density.info="none", srtCol=30,
          margins=c(5.5,7), cexRow=1.25,cexCol=1.25,
          key.xlab="1/RL Score",
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

musk_plot <- 
  as.data.frame(tax_table(muskegon_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(muskegonOTUs, by = "OTU") %>%
  rename(HNA = Muskegon_HNA, LNA = Muskegon_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.09) %>%
  ggplot(aes(y = RL_score, x = Rank5, fill = Rank2)) + 
  geom_bar(stat = "identity", position = "stack",color = "black") + 
  facet_wrap(~fcm_group, scale = "free_x") +  ggtitle("Muskegon Lake") +  
  scale_fill_manual(values=phylum_colors, name = "Phylum") + xlab("Family") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.3)) +
  theme(legend.position = "right", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))

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

inland_plot <- as.data.frame(tax_table(inland_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(inlandOTUs, by = "OTU") %>%
  dplyr::filter(Inland_HNA > 0.13 | Inland_LNA > 0.108) %>%
  rename(HNA = Inland_HNA, LNA = Inland_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.15) %>%
  ggplot(aes(y = RL_score, x = Rank4, fill = Rank2)) +
  geom_bar(stat = "identity", color = "black") + facet_wrap(~fcm_group, scale = "free_x") +
  scale_fill_manual(values=phylum_colors, name = "Phylum") + xlab("Order") +
  ggtitle("Inland Lakes") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  theme(legend.position = "right", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))

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

michi_plot <- as.data.frame(tax_table(michigan_taxonomy)) %>%
  tibble::rownames_to_column(var = "OTU") %>%
  full_join(michiganOTUs, by = "OTU") %>%
  dplyr::filter(Michigan_HNA > 0.248 | Michigan_LNA > 0.286) %>%
  rename(HNA = Michigan_HNA, LNA = Michigan_LNA) %>%
  gather(key = fcm_group, value = RL_score, HNA:LNA) %>%
  dplyr::filter(RL_score > 0.15) %>%
  ggplot(aes(y = RL_score, x = Rank5, fill = Rank2)) +
  geom_bar(stat = "identity", color = "black") + facet_wrap(~fcm_group, scale = "free_x") +
  scale_fill_manual(values=phylum_colors, name = "Phylum") + xlab("Family") +
  ggtitle("Lake Michigan") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  theme(legend.position = "right", axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))

plot_grid(musk_plot, inland_plot, michi_plot,
          nrow = 3, ncol = 1, labels = c("A", "B", "C"),
          align = "v")
ggsave("heatmap_figs/HNALNA_scores.jpg", width = 8, height = 12)



#### PHYLOGENETIC ANALYSIS
load("data/phyloseq.RData")
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

#write(vector_of_otus, file = "data/fasttree/OTUnames_based_on_RLscores.txt",
#      ncolumns = 1,
#      append = FALSE, sep = "\n")

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



tax <- data.frame(tax_table(final_physeq)) %>%
  tibble::rownames_to_column(var = "OTU")

fcm_groups_df <- read.csv("data/fasttree/OTUnames_based_on_RLscores_MANUAL.csv") %>%
  left_join(., tax, by = "OTU") %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::select(fcm_type)

p_fcm <- ggplot(HNALNA_otu_tree, aes(x, y)) + geom_tree() + theme_tree() +
  geom_tiplab(size=3, align=TRUE, linesize=.5) 
gheatmap(p_fcm, fcm_groups_df, offset = 0.2, width=0.15, font.size=0, colnames_angle=0, hjust=0.5) +
  scale_fill_manual(values = c("black", "deepskyblue4", "darkgoldenrod1"))
ggsave("heatmap_figs/phylogenetic_tree_fcm_only.jpg", width = 8, height = 8)




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



### Is there phylogentic signal?
physig_df <- data.frame(HNALNA_otu_tree$tip.label) %>%
  rename(OTU = HNALNA_otu_tree.tip.label) %>%
  left_join(., otu_scores_df, by = "OTU") %>%
  left_join(., tax, by = "OTU")

stopifnot(physig_df$OTU == HNALNA_otu_tree$tip.label)

# null hypothesis of no signal we just type:
phylosig(HNALNA_otu_tree, physig_df$Muskegon_HNA, method="K",test=TRUE)
phylosig(HNALNA_otu_tree, physig_df$Muskegon_LNA, method="K",test=TRUE)

phylosig(HNALNA_otu_tree, physig_df$Inland_HNA, method="K",test=TRUE)
phylosig(HNALNA_otu_tree, physig_df$Inland_LNA, method="K",test=TRUE)

phylosig(HNALNA_otu_tree, physig_df$Michigan_HNA, method="K",test=TRUE)
phylosig(HNALNA_otu_tree, physig_df$Michigan_LNA, method="K",test=TRUE)

phylosig(HNALNA_otu_tree, physig_df$Phylum, method="K",test=TRUE)
phylosig(HNALNA_otu_tree, physig_df$Michigan_LNA, method="K",test=TRUE)


######################################### OUTGROUP TREE

outgroup_tree <- read.tree(file="data/fasttree/outgroup_tree/RAxML_bipartitions.newick_tree_HNALNA_rmN_outgroup.tre")

outgroup_tree_tip_order <- data.frame(outgroup_tree$tip.label) %>%
  rename(OTU = outgroup_tree.tip.label)


outgroup_physeq <- merge_phyloseq(physeq, phy_tree(outgroup_tree))

# Fix the taxonomy names 
colnames(tax_table(outgroup_physeq)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

###################################################################### ADD THE PROTEOBACTERIA TO THE PHYLA
phy <- data.frame(tax_table(outgroup_physeq))
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

tax_table(outgroup_physeq) <- t

library(ggtree)

outgroup_tax <- data.frame(tax_table(outgroup_physeq)) %>%
  tibble::rownames_to_column(var = "OTU")

outgroup_fcm_groups_df <- read.csv("data/fasttree/OTUnames_based_on_RLscores_MANUAL.csv") %>%
  left_join(., outgroup_tax, by = "OTU") %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::select(fcm_type)

outgroup_p <- ggplot(outgroup_tree, aes(x, y)) + geom_tree() + theme_tree() +
  geom_tiplab(size=3, align=TRUE, linesize=.5) 
gheatmap(outgroup_p, outgroup_fcm_groups_df, offset = 0.2, width=0.15, font.size=0, colnames_angle=0, hjust=0.5) +
  scale_fill_manual(values = c("black", "deepskyblue4", "darkgoldenrod1"))
#ggsave("heatmap_figs/outgroup_phylogenetic_tree_fcm_only.jpg", width = 8, height = 8)



outgroup_df2 <- read.csv("data/fasttree/OTUnames_based_on_RLscores_MANUAL.csv") %>%
  left_join(., outgroup_tax, by = "OTU") %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::select(Phylum)

outgroup_p2 <- ggplot(outgroup_tree, aes(x, y)) + geom_tree() + theme_tree() +
  geom_tiplab(size=3, align=TRUE, linesize=.5) 
gheatmap(outgroup_p2, outgroup_df2, offset = 0.5, width=0.5, font.size=3, colnames_angle=0, hjust=0.5) +
  scale_fill_manual(values = phylum_colors)
#ggsave("heatmap_figs/outgroup_phylogenetic_tree_phylum.jpg", width = 8, height = 8)





######################################### FASTTREE
library(ape)
fast_tree <- read.tree(file="data/fasttree/fasttree_newick_tree_HNALNA_rmN.tre")

fast_tree_tip_order <- data.frame(fast_tree$tip.label) %>%
  rename(OTU = fast_tree.tip.label)


fasttree_physeq <- merge_phyloseq(physeq, phy_tree(fast_tree))

# Fix the taxonomy names 
colnames(tax_table(fasttree_physeq)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

###################################################################### ADD THE PROTEOBACTERIA TO THE PHYLA
phy <- data.frame(tax_table(fasttree_physeq))
Phylum <- as.character(phy$Phylum)
Class <- as.character(phy$Class)

for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    Phylum[i] <- Class[i]  } }

phy$Phylum <- Phylum # Add the new phylum level data back to phy

t <- tax_table(as.matrix(phy))
tax_table(fasttree_physeq) <- t

fasttree_tax <- data.frame(tax_table(fasttree_physeq)) %>%
  tibble::rownames_to_column(var = "OTU")

fasttree_df2 <- read.csv("data/fasttree/OTUnames_based_on_RLscores_MANUAL.csv") %>%
  left_join(., fasttree_tax, by = "OTU") %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::select(Phylum)


## Let's root the tree
is.rooted(fast_tree)
test_tree <- root(fast_tree, outgroup = "Otu001533", resolve.root = TRUE)
is.rooted(test_tree)


phyfcm_fasttree_df3 <- read.csv("data/fasttree/OTUnames_based_on_RLscores_MANUAL.csv") %>%
  left_join(., fasttree_tax, by = "OTU") %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::rename(FCM = fcm_type) %>%
  dplyr::select(Phylum, FCM)

rooted_tree <- 
  ggplot(test_tree, aes(x, y)) + geom_tree() + theme_tree() +
  geom_tiplab(size=3, align=TRUE, linesize=.5) #+
  #geom_nodelab(vjust=-.5, size=3) 

  
gheatmap(rooted_tree, phyfcm_fasttree_df3, offset = 0.1, width=0.3, font.size=3, colnames_angle=0, hjust=0.5) +
  scale_fill_manual(values = phylum_colors,  
                    breaks = c("Actinobacteria", "Alphaproteobacteria", "Bacteria_unclassified", "Bacteroidetes", "Betaproteobacteria", "Cyanobacteria",
                               "Deltaproteobacteria", "Firmicutes", "Gammaproteobacteria","Omnitrophica", "Planctomycetes", "Proteobacteria_unclassified",
                               "unknown_unclassified", "Verrucomicrobia", "Both", "HNA", "LNA"))
ggsave("heatmap_figs/rooted_fasttree_phylumFCM.jpg", width = 8, height = 8)






