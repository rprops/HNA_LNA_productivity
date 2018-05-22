# February 22th, 2018 
# Marian L Schmidt 

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

fcm_colors <- c(
  "HNA" = "deepskyblue4",
  "LNA" = "darkgoldenrod1",
  "Total" = "black")

# Read in the data 
updated_data <- read.table(file="data/Chloroplasts_removed/ByLake_Filtering/5in10/muskegon/muskegon_sampledata_5in10.tsv", header = TRUE) %>%
  mutate(norep_filter_name = paste(substr(Sample_16S,1,4), substr(Sample_16S,6,9), sep = "")) %>%
  arrange(norep_filter_name)


# Subset out only the Muskegon and Surface samples 
muskegon_data <- updated_data %>%
  dplyr::filter(Lake == "Muskegon" & Depth == "Surface") %>%
  dplyr::select(norep_filter_name)


# Load in the productivity data
production <- read.csv(file="data/production_data.csv", header = TRUE) %>%    
  dplyr::filter(fraction == "Free" & limnion == "Surface") %>%                         # Select only rows that are free-living
  dplyr::select(names, tot_bacprod, SD_tot_bacprod) %>%         # Select relevant columns for total productivity
  mutate(tot_bacprod = round(tot_bacprod, digits = 2),          # Round to 2 decimals
         SD_tot_bacprod = round(SD_tot_bacprod, digits = 2) ) %>%      
  dplyr::rename(norep_filter_name = names) %>%                  # Rename to match other data frame
  arrange(norep_filter_name)

# Stop if the names do not match
stopifnot(muskegon_data$norep_filter_name == production$norep_filter_name)

# Combine the two muskegon data frames into one
combined_data <- left_join(muskegon_data, production, by = "norep_filter_name")

# Merge the combined data back into the original data frame (data)
data <- left_join(updated_data, combined_data, by = "norep_filter_name")

# Write out the file 
#write.table(data, file="data/Chloroplasts_removed/productivity_data.tsv", row.names=TRUE)

###### SUBSET MUSKEGON ONLY
muskegon <- dplyr::filter(data, Lake == "Muskegon" & Depth == "Surface") %>%
  mutate(Site = factor(Site, levels = c("MOT", "MDP", "MBR", "MIN"))) %>%
  filter(tot_bacprod < 90)

#############
#############
# Number of Cells in Muskegon
musk_cells <- muskegon %>%
  dplyr::select(1:4, Season:Depth) %>%
  rename(Total = Total.cells, HNA = HNA.cells, LNA = LNA.cells) %>%
  gather(key = FCM_type, value = num_cells, Total:LNA) 

df_cells <- data %>%
  dplyr::select(1:4, Lake, Season:Depth) %>%
  rename(Total = Total.cells, HNA = HNA.cells, LNA = LNA.cells) %>%
  gather(key = FCM_type, value = num_cells, Total:LNA)
  
stats_cells <- df_cells %>%
  group_by(Lake, FCM_type) %>%
  summarise(count = n(),
            mean_cells = mean(num_cells),
            median_cells = median(num_cells))

proportionHNA_plot <- 
  df_cells %>%
  dplyr::select(samples, Lake, FCM_type, num_cells) %>%
  spread(FCM_type, num_cells) %>%
  mutate(prop_HNA = HNA/Total * 100,
         prop_LNA = LNA/Total * 100) %>%
  dplyr::select(Lake, prop_HNA, prop_LNA) %>%
  rename(HNA = prop_HNA, LNA = prop_LNA) %>%
  gather(key = fcm_type, value = Percentage, HNA:LNA) %>%
  ggplot(aes(x = Lake, y = Percentage, color = fcm_type, fill = fcm_type)) +
    facet_grid(~fcm_type) + ylab("Percentage of Total Cells") +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(alpha = 0.9) +
    scale_color_manual(values = fcm_colors) +
    scale_fill_manual(values = fcm_colors) + 
    theme(legend.position = "none")
#ggsave(filename = "data/Chloroplasts_removed/HNA-LNA-proportions.png", 
#       width = 6, height = 3.5, units = "in", dpi = 500)



df_cells %>%
  dplyr::select(samples, Lake, FCM_type, num_cells) %>%
  spread(FCM_type, num_cells) %>%
  mutate(prop_HNA = HNA/Total * 100,
         prop_LNA = LNA/Total * 100) %>%
  dplyr::select(Lake, prop_HNA, prop_LNA) %>%
  rename(HNA = prop_HNA, LNA = prop_LNA) %>%
  group_by(Lake) %>%
  summarize(min_HNA = min(HNA), 
            max_HNA = max(HNA), 
            mean_HNA = mean(HNA),
            min_LNA = min(LNA), 
            max_LNA = max(LNA), 
            mean_LNA = mean(LNA))

df_cells %>%
  dplyr::select(samples, Lake, FCM_type, num_cells) %>%
  spread(FCM_type, num_cells) %>%
  mutate(prop_HNA = HNA/Total * 100,
         prop_LNA = LNA/Total * 100) %>%
  dplyr::select(Lake, prop_HNA, prop_LNA) %>%
  rename(HNA = prop_HNA, LNA = prop_LNA) %>%
  summarize(mean_HNA = mean(HNA), mean_LNA = mean(LNA))



df_cells %>%
  dplyr::select(samples, Lake, FCM_type, num_cells) %>%
  spread(FCM_type, num_cells) %>%
  mutate(prop_HNA = HNA/Total * 100,
         prop_LNA = LNA/Total * 100) %>%
  dplyr::select(Lake, prop_HNA, prop_LNA) %>%
  group_by(Lake) %>%
  summarize(mean_prop_HNA = mean(prop_HNA),
            mean_prop_LNA = mean(prop_LNA))

# DESCRIPTIVE STATS
totcells_df <- df_cells %>%
  filter(FCM_type == "Total")
# Compute the analysis of variance
totcells_aov <- aov(num_cells ~ Lake, data = totcells_df)
summary(totcells_aov)
# Which samples are significant from each other?
TukeyHSD(totcells_aov)


cells_boxplot <- 
  ggplot(df_cells, 
         aes(x = FCM_type, y = num_cells, fill = FCM_type, color = FCM_type, shape = Lake)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = NA, show.legend = FALSE) +
  geom_point(size = 1.5, position = position_jitterdodge()) + 
  scale_fill_manual(values = fcm_colors, guide = "none") + 
  scale_color_manual(values = fcm_colors, guide = "none") + 
  labs(y = "Number of Cells \n (cells/mL)", x = "Cell Type") +
  theme(legend.position = c(0.04, 0.85))

ggplot(df_cells, 
       aes(x = Lake, y = num_cells, fill = FCM_type, color = FCM_type)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = NA, show.legend = FALSE) +
  geom_point(size = 1, position = position_jitterdodge()) + 
  #facet_grid(FCM_type~., scales = "free") + 
  scale_fill_manual(values = fcm_colors, guide = "none") + 
  scale_color_manual(values = fcm_colors) + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(y = "Number of Cells \n (cells/mL)", x = "Lake") +
  theme(legend.position = "top", legend.title = element_blank())

#ggsave(filename = "data/Chloroplasts_removed/HNA-LNA-lakes.png", 
#       width = 4, height = 3.5, units = "in", dpi = 500)

####################################################################################
####################################################################################
########################  Analysis of HNA & LNA Correlations

# 1. Run the linear model 
lm_NA_corr <- lm(LNA.cells ~ HNA.cells, data = muskegon)

## 2. Extract the R2 and p-value from the linear model: 
lm_NA_corr_stats <- paste("atop(R^2 ==", round(summary(lm_NA_corr)$adj.r.squared, digits = 3), ",",
                      "p ==", round(unname(summary(lm_NA_corr)$coefficients[,4][2]), digits = 3), ")")

# 3. Plot it
ggplot(muskegon, aes(x = HNA.cells, y = LNA.cells)) + 
  geom_errorbarh(aes(xmin = HNA.cells - HNA.sd, xmax = HNA.cells + HNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = LNA.cells - LNA.sd, max = LNA.cells + LNA.sd), color = "grey") + 
  geom_point(size = 3, shape = 22, fill = "black") + 
  geom_smooth(method = "lm", color = "black") + 
  labs(y = "LNA Cell Density (cells/mL)", x = "HNA Cell Density (cells/mL)") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), 
                     breaks = c(2e+06, 3e+06)) +
  annotate("text", x=3e+06, y=3e+06, label=lm_NA_corr_stats, parse = TRUE, color = "black", size = 4) 

ggsave(filename = "data/Chloroplasts_removed/HNA-LNA-correlation.png", 
       width = 4, height = 3.5, units = "in", dpi = 500)


####################################################################################
####################################################################################
########################  Analysis of HNA/LNA/Total Cells vs Total Productivity
# 1. Run the linear model 
lm_HNA <- lm(tot_bacprod ~ HNA.cells, data = muskegon)

## 2. Extract the R2 and p-value from the linear model: 
lm_HNA_stats <- paste("atop(R^2 ==", round(summary(lm_HNA)$adj.r.squared, digits = 2), ",",
                                    "p ==", round(unname(summary(lm_HNA)$coefficients[,4][2]), digits = 5), ")")

# 3. Plot it
HNA_vs_prod <- ggplot(muskegon, aes(x = HNA.cells, y = tot_bacprod)) + 
  geom_errorbarh(aes(xmin = HNA.cells - HNA.sd, xmax = HNA.cells + HNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = tot_bacprod - SD_tot_bacprod, max = tot_bacprod + SD_tot_bacprod), color = "grey") + 
  geom_point(size = 3, shape = 22, fill = "deepskyblue4") + 
  geom_smooth(method = "lm", color = "deepskyblue4") + 
  labs(y = "Total Bacterial Production \n (μg C/L/day)", x = "HNA Cell Density (cells/mL)") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), 
                     breaks = c(2e+06, 3e+06)) +
  annotate("text", x=1.5e+06, y=60, label=lm_HNA_stats, parse = TRUE, color = "black", size = 4) 


# 1. Run the linear model 
lm_LNA <- lm(tot_bacprod ~ LNA.cells, data = muskegon)

## 2. Extract the R2 and p-value from the linear model: 
lm_LNA_stats <- paste("atop(R^2 ==", round(summary(lm_LNA)$adj.r.squared, digits = 3), ",",
                      "p ==", round(unname(summary(lm_LNA)$coefficients[,4][2]), digits = 2), ")")

# 3. Plot it
LNA_vs_prod <- ggplot(muskegon, aes(x = LNA.cells, y = tot_bacprod)) + 
  geom_errorbarh(aes(xmin = LNA.cells - LNA.sd, xmax = LNA.cells + LNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = tot_bacprod - SD_tot_bacprod, max = tot_bacprod + SD_tot_bacprod), color = "grey") + 
  geom_point(size = 3, shape = 22, fill = "darkgoldenrod1") + 
  labs(y = "Total Bacterial Production \n (μg C/L/day)", x = "LNA Cell Density (cells/mL)") +
  geom_smooth(method = "lm", se = FALSE, linetype = "longdash", color = "darkgoldenrod1") + 
  annotate("text", x=2.5e+06, y=60, label=lm_LNA_stats, parse = TRUE, color = "red", size = 4) 



# 1. Run the linear model 
lm_total <- lm(tot_bacprod ~ Total.cells, data = muskegon)

## 2. Extract the R2 and p-value from the linear model: 
lm_total_stats <- paste("atop(R^2 ==", round(summary(lm_total)$adj.r.squared, digits = 2), ",",
                      "p ==", round(unname(summary(lm_total)$coefficients[,4][2]), digits = 2), ")")

# 3. Plot it 
Total_vs_prod <- ggplot(muskegon, aes(x = Total.cells, y = tot_bacprod)) + 
  geom_errorbarh(aes(xmin = Total.cells - Total.count.sd, xmax = Total.cells + Total.count.sd), color = "grey") + 
  geom_errorbar(aes(ymin = tot_bacprod - SD_tot_bacprod, max = tot_bacprod + SD_tot_bacprod), color = "grey") + 
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_point(size = 3, shape = 22, fill = "black") + 
  labs(y = "Total Bacterial Production \n (μg C/L/day)", x = "Cell Density (cells/mL)") +
  geom_smooth(method = "lm", color = "black") + 
  #geom_smooth(method = "lm", se = FALSE, linetype = "longdash", color = "red") + 
  annotate("text", x=5e+06, y=60, label=lm_total_stats, parse = TRUE, color = "black", size = 4) 


# Put all three plots together into one 
totprod_plots <- plot_grid(cells_boxplot + theme(legend.text = element_text(size = 10),
                                                 legend.title = element_text(size = 11, face = "bold")),
                           HNA_vs_prod + theme(legend.position = "none"), 
                           LNA_vs_prod + theme(axis.title.y = element_blank(), legend.position = "none"), 
                           Total_vs_prod + theme(axis.title.y = element_blank(), legend.position = "none"), 
          labels = c("A", "B", "C", "D"), 
          ncol = 4, rel_widths = c(1, 0.9, 0.8, 0.8))

totprod_plots

## Plot the fraction of HNA
fmusk <- muskegon %>% mutate(fHNA = HNA.cells/Total.cells)

lm_fHNA <- lm(tot_bacprod ~ fHNA, data = fmusk)

## 2. Extract the R2 and p-value from the linear model: 
lm_fHNA_stats <- paste("atop(R^2 ==", round(summary(lm_fHNA)$adj.r.squared, digits = 2), ",",
                      "p ==", round(unname(summary(lm_fHNA)$coefficients[,4][2]), digits = 3), ")")


fHNA_vs_prod <- ggplot(fmusk, aes(x = fHNA, y = tot_bacprod)) + 
  geom_errorbar(aes(ymin = tot_bacprod - SD_tot_bacprod, max = tot_bacprod + SD_tot_bacprod), color = "grey") + 
  scale_shape_manual(values = c(21, 22, 23, 24)) + 
  geom_point(size = 3, aes(shape = Site), fill = "black") + 
  geom_smooth(method = "lm") + 
  ylab("Bacterial Production") + xlab("Fraction HNA") +
  annotate("text", x= 0.22, y=60, label=lm_fHNA_stats, parse = TRUE, color = "black", size = 4) 

# Save the plot to a file to call in the README
ggsave(plot = totprod_plots, 
       filename = "data/Chloroplasts_removed/HNA_vs_productivity.png", 
       width = 14, height = 4, units = "in", dpi = 500)


fHNA_comparison <- plot_grid(HNA_vs_prod + theme(legend.position = "none"), 
          fHNA_vs_prod + theme(legend.position = "none"),
          labels = c("A", "B"), ncol = 2) 

ggsave(plot = fHNA_comparison, 
       filename = "data/Chloroplasts_removed/fHNA_vs_productivity.png", 
       width = 7, height = 3.5, units = "in", dpi = 500)



####################################################################################
####################################################################################
######################## CORRELATED HNA AND LNA

musk_all_df <- read.table(file="data/Chloroplasts_removed/ByLake_Filtering/5in10/muskegon/muskegon_sampledata_5in10.tsv", header = TRUE) %>%
  mutate(norep_filter_name = paste(substr(Sample_16S,1,4), substr(Sample_16S,6,9), sep = "")) %>%
  arrange(norep_filter_name)

mich_all_df <- read.table(file="data/Chloroplasts_removed/ByLake_Filtering/5in10/michigan/michigan_sampledata_5in10.tsv", header = TRUE) %>%
  mutate(norep_filter_name = paste(substr(Sample_16S,1,4), substr(Sample_16S,6,9), sep = "")) %>%
  arrange(norep_filter_name)

inla_all_df <- read.table(file="data/Chloroplasts_removed/ByLake_Filtering/5in10/inland/inland_sampledata_5in10.tsv", header = TRUE) %>%
  mutate(norep_filter_name = paste(substr(Sample_16S,1,4), substr(Sample_16S,6,9), sep = "")) %>%
  arrange(norep_filter_name)

stopifnot(colnames(musk_all_df) == colnames(mich_all_df))
stopifnot(colnames(musk_all_df) == colnames(inla_all_df))

lakes <- bind_rows(musk_all_df, mich_all_df, inla_all_df)

# 1. Run the linear model 
lm_allNA_corr <- lm(LNA.cells ~ HNA.cells, data = lakes)
summary(lm(LNA.cells ~ HNA.cells * Lake, data = lakes))

## 2. Extract the R2 and p-value from the linear model: 
lm_allNA_corr_stats <- paste("atop(R^2 ==", round(summary(lm_allNA_corr)$adj.r.squared, digits = 2), ",",
                          "p ==", round(unname(summary(lm_allNA_corr)$coefficients[,4][2]), digits = 24), ")")

lake_colors <- c(
  "Muskegon" = "#FF933F",   #"#1AB58A",
  "Michigan" =  "#EC4863", #"#FFC543",
  "Inland" =  "#5C2849")  #"#FF2151")

mytheme <- theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8, face = "bold"), 
                 plot.title = element_text(size = 10),
                 axis.title = element_text(size = 10, face = "bold"), axis.text = element_text(size = 8),
                 legend.key.width=unit(0.1,"line"), legend.key.height=unit(0.1,"line")) 



# 3. Plot it
ggplot(lakes, aes(x = HNA.cells, y = LNA.cells)) + 
  geom_errorbarh(aes(xmin = HNA.cells - HNA.sd, xmax = HNA.cells + HNA.sd), color = "grey", alpha = 0.8) + 
  geom_errorbar(aes(ymin = LNA.cells - LNA.sd, max = LNA.cells + LNA.sd), color = "grey", alpha = 0.8) + 
  geom_point(size = 2.5, shape = 22, alpha = 0.8, aes(fill = Lake)) + 
  geom_smooth(method = "lm", color = "black") + 
  scale_fill_manual(values = lake_colors) +
  labs(y = "LNA Cell Density (cells/mL)", x = "HNA Cell Density (cells/mL)") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 6.1e+06)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 1e+07)) +
  annotate("text", x=5e+06, y=0.8e+06, label=lm_allNA_corr_stats, parse = TRUE, color = "black", size = 3) +
  mytheme + theme(legend.position = c(0.01, 0.9))

ggsave(filename = "data/Chloroplasts_removed/AllLakes-HNA-LNA-correlation.png", 
       width = 4, height = 3.5, units = "in", dpi = 500)



### MUSKEGON
# 1. Run the linear model 
lm_muskNA_corr <- lm(LNA.cells ~ HNA.cells, data = musk_all_df)

## 2. Extract the R2 and p-value from the linear model: 
lm_muskNA_corr_stats <- paste("atop(R^2 ==", round(summary(lm_muskNA_corr)$adj.r.squared, digits = 2), ",",
                             "p ==", round(unname(summary(lm_muskNA_corr)$coefficients[,4][2]), digits = 9), ")")

musk_corr_plot <- ggplot(musk_all_df, aes(x = HNA.cells, y = LNA.cells)) + 
  geom_errorbarh(aes(xmin = HNA.cells - HNA.sd, xmax = HNA.cells + HNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = LNA.cells - LNA.sd, max = LNA.cells + LNA.sd), color = "grey") + 
  geom_point(size = 3, shape = 22, aes(fill = Lake)) + 
  geom_smooth(method = "lm", color = "black") + 
  ggtitle("Muskegon Lake") + scale_fill_manual(values = lake_colors) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 6.1e+06)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 1e+07)) +
  labs(y = "LNA Cell Density (cells/mL)", x = "HNA Cell Density (cells/mL)") +
  annotate("text", x=5e+06, y=0.8e+06, label=lm_muskNA_corr_stats, parse = TRUE, color = "black", size = 3) +
  mytheme

### INLAND
# 1. Run the linear model 
lm_inlaNA_corr <- lm(LNA.cells ~ HNA.cells, data = filter(inla_all_df, Sample_16S != "Z14003F"))

## 2. Extract the R2 and p-value from the linear model: 
lm_inlaNA_corr_stats <- paste("atop(R^2 ==", round(summary(lm_inlaNA_corr)$adj.r.squared, digits = 2), ",",
                              "p ==", round(unname(summary(lm_inlaNA_corr)$coefficients[,4][2]), digits = 2), ")")

inla_corr_plot <- ggplot(filter(inla_all_df, Sample_16S != "Z14003F"), aes(x = HNA.cells, y = LNA.cells)) + 
  geom_errorbarh(aes(xmin = HNA.cells - HNA.sd, xmax = HNA.cells + HNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = LNA.cells - LNA.sd, max = LNA.cells + LNA.sd), color = "grey") + 
  geom_point(size = 3, shape = 22, aes(fill = Lake)) + 
  geom_smooth(method = "lm", color = "black") + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 6.1e+06)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 1e+07)) +
  ggtitle("Inland Lakes") + scale_fill_manual(values = lake_colors) +
  labs(y = "LNA Cell Density (cells/mL)", x = "HNA Cell Density (cells/mL)") +
  annotate("text", x=5e+06, y=0.8e+06, label=lm_inlaNA_corr_stats, parse = TRUE, color = "black", size = 3) +
  mytheme

### MICHIGAN
# 1. Run the linear model 
lm_michNA_corr <- lm(LNA.cells ~ HNA.cells, data = mich_all_df)

## 2. Extract the R2 and p-value from the linear model: 
lm_michNA_corr_stats <- paste("atop(R^2 ==", round(summary(lm_michNA_corr)$adj.r.squared, digits = 2), ",",
                              "p ==", round(unname(summary(lm_michNA_corr)$coefficients[,4][2]), digits = 11), ")")

mich_corr_plot <- ggplot(mich_all_df, aes(x = HNA.cells, y = LNA.cells)) + 
  geom_errorbarh(aes(xmin = HNA.cells - HNA.sd, xmax = HNA.cells + HNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = LNA.cells - LNA.sd, max = LNA.cells + LNA.sd), color = "grey") + 
  geom_point(size = 3, shape = 22, aes(fill = Lake)) + 
  geom_smooth(method = "lm", color = "black") + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 6.1e+06)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 1e+07)) +
  ggtitle("Lake Michigan") + scale_fill_manual(values = lake_colors) +
  labs(y = "LNA Cell Density (cells/mL)", x = "HNA Cell Density (cells/mL)") +
  annotate("text", x=5e+06, y=0.8e+06, label=lm_michNA_corr_stats, parse = TRUE, color = "black", size = 3) +
  mytheme

plot_grid(mich_corr_plot, inla_corr_plot, musk_corr_plot, 
          labels = c("A", "B", "C"), 
          nrow = 1, ncol = 3)

ggsave(filename = "data/Chloroplasts_removed/AllLakes-separate-HNA-LNA-correlation.png", 
       width = 10, height = 3.5, units = "in", dpi = 500)