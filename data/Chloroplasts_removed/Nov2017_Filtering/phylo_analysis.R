# Name: Marian L. Schmidt
# Date: November 26th, 2017
# Purpose: Phylogenetic diversity Analysis

###################### Load Libraries
library(phyloseq)
library(picante)
library(tidyr)
library(dplyr)


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
## TBD


