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



scores_df <- bind_rows(HNA, LNA) %>%
  dplyr::select(OTU, RL.ranking, FCM_type)
