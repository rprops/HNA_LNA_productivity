# Name: Marian L. Schmidt
# Date: November 26th, 2017
# Purpose: Phylogenetic diversity Analysis

###################### Load Libraries
library(phyloseq)
library(picante)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggtree)


###################### Load Data
# Read in the data with 2 data objects 
    # 1. ALL_rel_physeq_1in3
    # 2. ALL_rel_physeq_5in10
load("./data/Chloroplasts_removed/Nov2017_Filtering/ALL-physeq-for-phylo.RData")

########### 5in10
mpd_5in10 <- read.table("./data/Chloroplasts_removed/Nov2017_Filtering/unweighted_MPD_5seqs_in_10perc.tsv") %>%
  tibble::rownames_to_column("Sample_16S")

metadata_5in10_raw <- read.table("./data/Chloroplasts_removed/Nov2017_Filtering/5seqs_in_10percent_samples/nochloro_HNA_LNA_5in10percent.tsv") %>%
  mutate(type = "ALL")
row.names(metadata_5in10) <- NULL

metadata_5in10_LNA <- metadata_5in10_raw %>%
  mutate(Sample_16S = paste(Sample_16S, "_LNA", sep=""),
         type = "LNA")

metadata_5in10_HNA <- metadata_5in10_raw %>%
  mutate(Sample_16S = paste(Sample_16S, "_HNA", sep=""),
         type = "HNA")

metadata_5in10 <- bind_rows(metadata_5in10_raw, metadata_5in10_LNA, metadata_5in10_HNA) 

mpd_5in10_data <-left_join(mpd_5in10, metadata_5in10, by = "Sample_16S") %>%
  mutate(type = factor(type, levels = c("ALL", "LNA", "HNA")),
         threshold = "5in10")


########### 1in3
mpd_1in3 <- read.table("./data/Chloroplasts_removed/Nov2017_Filtering/unweighted_MPD_1seqs_in_3samps.tsv") %>%
  tibble::rownames_to_column("Sample_16S")
  

metadata_1in3_raw <- read.table("./data/Chloroplasts_removed/Nov2017_Filtering/1seq_in_3samples/nochloro_HNA_LNA_1seqin3samps.tsv") %>%
  mutate(type = "ALL")
row.names(metadata_1in3_raw) <- NULL

metadata_1in3_LNA <- metadata_1in3_raw %>%
  mutate(Sample_16S = paste(Sample_16S, "_LNA", sep=""),
         type = "LNA")

metadata_1in3_HNA <- metadata_1in3_raw %>%
  mutate(Sample_16S = paste(Sample_16S, "_HNA", sep=""),
         type = "HNA")

metadata_1in3 <- bind_rows(metadata_1in3_raw, metadata_1in3_LNA, metadata_1in3_HNA)

mpd_1in3_data <-left_join(mpd_1in3, metadata_1in3, by = "Sample_16S") %>%
  mutate(type = factor(type, levels = c("ALL", "LNA", "HNA")),
         threshold = "1in3")

all_mpd_data <- bind_rows(mpd_5in10_data, mpd_1in3_data)


ggplot(all_mpd_data, aes(y = mpd.obs.z, x = type, color = type, fill = type)) +
  theme_bw() + ylab("Phylogenetic Diversity \n Negative = Clustered; Positive = Overdispersed") + 
  facet_grid(threshold~Lake, scales = "free") + 
  geom_jitter(alpha = 0.9) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#F8BB04", "#ED8137", "#FF4E50")) +
  scale_color_manual(values = c("#F8BB04", "#ED8137", "#FF4E50")) +
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank(), legend.position = "bottom")
ggsave("./data/Chloroplasts_removed/Nov2017_Filtering/phylo_analysis_figures/PD_by_lake.png",  width = 6, height = 5, units = "in")

ggplot(all_mpd_data, aes(y = ntaxa, x = type, color = type, fill = type)) +
  theme_bw() + ylab("Number of Taxa") + 
  facet_grid(threshold~Lake, scales = "free") + 
  geom_jitter(alpha = 0.9) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#F8BB04", "#ED8137", "#FF4E50")) +
  scale_color_manual(values = c("#F8BB04", "#ED8137", "#FF4E50")) +
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank(), legend.position = "bottom")
ggsave("./data/Chloroplasts_removed/Nov2017_Filtering/phylo_analysis_figures/ntaxa_by_lake.png", width = 6, height = 5, units = "in")


## Muskegon Lake Only 
ggplot(filter(all_mpd_data, Lake == "Muskegon" & Year == "2015" & Depth == "Surface"), aes(y = mpd.obs.z, x = type, color = type, fill = type)) +
  theme_bw() + ylab("Phylogenetic Diversity \n Negative = Clustered; Positive = Overdispersed") + 
  facet_grid(.~threshold) + 
  geom_jitter(alpha = 0.9) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#F8BB04", "#ED8137", "#FF4E50")) +
  scale_color_manual(values = c("#F8BB04", "#ED8137", "#FF4E50")) +
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank(), legend.position = "bottom")
ggsave("./data/Chloroplasts_removed/Nov2017_Filtering/phylo_analysis_figures/PD_musk2015_surface.png", width = 5, height = 4, units = "in")



##### PHYLO TREE
## for 5in10 

# Which OTUs are only found in HNA, LNA, or both?
tree_5in10 <- read.tree(file = "./data/Chloroplasts_removed/Nov2017_Filtering/5seqs_in_10percent_samples/newick_tree_5seqs_in_10perc_rmN.tre")

load("./data/Chloroplasts_removed/Nov2017_Filtering/HNA-LNA-physeq_5in10.RData")
# 2 objects named HNA_physeq_5in10 and LNA_physeq_5in10

load("./data/Chloroplasts_removed/Nov2017_Filtering/ALL-physeq-for-phylo.RData")
# ALL_rel_physeq_5in10


all_otus_5in10 <- taxa_names(ALL_rel_physeq_5in10)

df_5in10 <- data.frame(matrix(ncol = 1, nrow = length(all_otus_5in10)))
colnames(df_5in10)[1] <- "OTU"
df_5in10 <- mutate(df_5in10, OTU = all_otus_5in10)


hna_otus_5in10 <- data.frame(taxa_names(HNA_physeq_5in10)) %>% 
  mutate(HNA = "YES")
colnames(hna_otus_5in10)[1] <- "OTU"


lna_otus_5in10 <- data.frame(taxa_names(LNA_physeq_5in10)) %>% 
  mutate(LNA = "YES")
colnames(lna_otus_5in10)[1] <- "OTU"


otus_5in10 <- left_join(df_5in10, hna_otus_5in10, by = "OTU") %>%
  left_join(lna_otus_5in10, by = "OTU") %>%
  mutate(HNA = ifelse(is.na(HNA), "NO", "YES"),
         LNA = ifelse(is.na(LNA), "NO", "YES")) %>%
  mutate(Type = ifelse(HNA == "YES" & LNA == "YES", "BOTH", 
                       ifelse(HNA == "YES" & LNA == "NO", "HNA", 
                              ifelse(HNA == "NO" & LNA == "YES", "LNA", 
                                     ifelse(HNA == "NO" & LNA == "NO", "NEITHER",
                                            "PANIC"))))) 

plot_percentOTUs_5in10 <- otus_5in10 %>%
  group_by(Type) %>%
  count() %>%
  mutate(Total_OTU_Percent = n/486) %>%
  ggplot(aes(x = Type, y = Total_OTU_Percent, fill = Type)) + theme_bw() + 
  ylab("Percent of Total OTUs") + ggtitle("5in10") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_fill_manual(values = c("#383245", "#41505E", "#569492", "#A4BAA2")) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  theme(legend.position = "bottom", legend.title = element_blank())





## for 1in3 

# Which OTUs are only found in HNA, LNA, or both? 
load("./data/Chloroplasts_removed/Nov2017_Filtering/HNA-LNA-physeq_1in3.RData")
# 2 objects named HNA_physeq_1in3 and LNA_physeq_1in3
# ALL_rel_physeq_1in3


all_otus_1in3 <- taxa_names(ALL_rel_physeq_1in3)

df_1in3 <- data.frame(matrix(ncol = 1, nrow = length(all_otus_1in3)))
colnames(df_1in3)[1] <- "OTU"
df_1in3 <- mutate(df_1in3, OTU = all_otus_1in3)


hna_otus_1in3 <- data.frame(taxa_names(HNA_physeq_1in3)) %>% 
  mutate(HNA = "YES")
colnames(hna_otus_1in3)[1] <- "OTU"


lna_otus_1in3 <- data.frame(taxa_names(LNA_physeq_1in3)) %>% 
  mutate(LNA = "YES")
colnames(lna_otus_1in3)[1] <- "OTU"


otus_1in3 <- left_join(df_1in3, hna_otus_1in3, by = "OTU") %>%
  left_join(lna_otus_1in3, by = "OTU") %>%
  mutate(HNA = ifelse(is.na(HNA), "NO", "YES"),
         LNA = ifelse(is.na(LNA), "NO", "YES")) %>%
  mutate(Type = ifelse(HNA == "YES" & LNA == "YES", "BOTH", 
                       ifelse(HNA == "YES" & LNA == "NO", "HNA", 
                              ifelse(HNA == "NO" & LNA == "YES", "LNA", 
                                     ifelse(HNA == "NO" & LNA == "NO", "NEITHER",
                                            "PANIC"))))) 

plot_percentOTUs_1in3 <- otus_1in3 %>%
  group_by(Type) %>%
  count() %>%
  mutate(Total_OTU_Percent = n/2951) %>%
  ggplot(aes(x = Type, y = Total_OTU_Percent, fill = Type)) + theme_bw() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  ylab("Percent of Total OTUs") + ggtitle("1in3") +
  scale_fill_manual(values = c("#383245", "#41505E", "#569492", "#A4BAA2")) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  theme(legend.position = "bottom", legend.title = element_blank())


library(cowplot)
toprow <- plot_grid(plot_percentOTUs_5in10 + theme(legend.position = "none"), 
                    plot_percentOTUs_1in3 + theme(legend.position = "none"), 
                    labels=c('A', 'B'), align="h")

type_legend <- get_legend(plot_percentOTUs_1in3)

plot_grid(toprow, type_legend, nrow = 2, ncol = 1, rel_heights = c(1, 0.1))



detach("package:phyloseq", unload=TRUE)

################## GGTREE
tree_5in10 <- read.tree(file = "./data/Chloroplasts_removed/Nov2017_Filtering/5seqs_in_10percent_samples/newick_tree_5seqs_in_10perc_rmN.tre")
hna_5in10_scores <- read.csv("./Scores/hnascores_otus_tuned_thr_0.36_5seq10_rel.csv") %>%
  dplyr::rename(OTU = X)
lna_5in10_scores <- read.csv("./Scores/lnascores_otus_tuned_thr_0.21_5seq10_rel.csv")


otus_5in10_rename <- otus_5in10 %>%
  dplyr::rename(label = OTU)

tree_5in10_df <- fortify(tree_5in10) %>%
  left_join(otus_5in10_rename, by = "label") %>%
  mutate(Type = factor(Type, levels = c("BOTH", "HNA", "LNA","NEITHER", NA)))

Type_colors <- c("BOTH" = "#531A4D", "HNA" = "#FD823F", "LNA" = "cornflowerblue", "NEITHER" = "#2F3440", "NA" = "black")

ggtree(tree_5in10_df, aes(color = Type)) +
  theme_tree2() + geom_tiplab(size=2) +
#  scale_color_manual(values = Type_colors)
  theme(legend.position = "right") 

otus_5in10_plotting <- otus_5in10 %>%
  tibble::column_to_rownames("OTU") %>%
  dplyr::select(Type)


ggtree(tree_5in10_df, branch.length = 'none') +
  geom_tiplab(size=2, aes(color = Type)) +
  scale_color_manual(breaks=c("BOTH", "HNA", "LNA","NEITHER"), 
                    values=c("firebrick","gold","darkslategray1", "white"))


p <- ggtree(tree_5in10_df, layout="circular") +
  theme_tree2() + geom_tiplab(size=2) +
  theme(legend.position = "right") 
gheatmap(p, otus_5in10_plotting, offset = 5, width=0.5, font.size=3, colnames_angle=-90, hjust=0) +
  scale_fill_manual(breaks=c("BOTH", "HNA", "LNA","NEITHER"), 
                    values=c("firebrick","gold","darkslategray1", "white"))
ggsave("./data/Chloroplasts_removed/Nov2017_Filtering/phylo_analysis_figures/phylo_tree_type_cir.png",  
       width = 10, height = 10, units = "in")

#scale_fill_manual(breaks=c("NO", "YES", "BOTH", "HNA", "LNA","NEITHER"), 
#                  values=c("gold","lawngreen", "black","gray89","darkslategray1", "firebrick"))



############### 1in3
tree_1in3 <- read.tree(file = "./data/Chloroplasts_removed/Nov2017_Filtering/1seq_in_3samples/newick_tree_1seqs_in_3samps_rmN.tre")
hna_1in3_scores <- read.csv("./Scores/hnascores_otus_tuned_thr_0.15_1seq3_rel.csv") %>%
  dplyr::rename(OTU = X)
lna_1in3_scores <- read.csv("./Scores/lnascores_otus_tuned_thr_0.18_1seq3_rel.csv")



otus_1in3_rename <- otus_1in3 %>%
  dplyr::rename(label = OTU)

tree_1in3_df <- fortify(tree_1in3) %>%
  left_join(otus_1in3_rename, by = "label") %>%
  mutate(Type = factor(Type, levels = c("BOTH", "HNA", "LNA","NEITHER", NA)))

ggtree(tree_1in3_df, aes(color = Type)) +
  theme_tree2() + geom_tiplab(size=2) +
  theme(legend.position = "right") 

