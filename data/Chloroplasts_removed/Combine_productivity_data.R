# March 4th, 2017 
# Marian L Schmidt 

# Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

# Read in the data 
raw_data <- read.table(file="data/Chloroplasts_removed/nochloro_HNA_LNA.tsv", header = TRUE)  %>%
  # Create matching column to go with Muskegon productivity data
  mutate(norep_filter_name = paste(substr(Sample_16S,1,4), substr(Sample_16S,6,9), sep = "")) %>%
  arrange(norep_filter_name)

# Subset out only the Muskegon and Surface samples 
muskegon_data <- raw_data %>%
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
data <- left_join(raw_data, combined_data, by = "norep_filter_name")

# Write out the file 
#write.table(data, file="data/Chloroplasts_removed/productivity_data.tsv", row.names=TRUE)




####################################################################################
####################################################################################
########################  Analysis of HNA/LNA/Total Cells vs Total Productivity
muskegon <- dplyr::filter(data, Lake == "Muskegon" & Depth == "Surface") %>%
  mutate(Site = factor(Site, levels = c("MOT", "MDP", "MBR", "MIN"))) 


lm_HNA <- lm(tot_bacprod ~ HNA.cells, data = muskegon)

HNA_vs_prod <- ggplot(muskegon, aes(x = HNA.cells, y = tot_bacprod)) + 
  geom_errorbarh(aes(xmin = HNA.cells - HNA.sd, xmax = HNA.cells + HNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = tot_bacprod - SD_tot_bacprod, max = tot_bacprod + SD_tot_bacprod), color = "grey") + 
  scale_shape_manual(values = c(21, 22, 23, 24)) + 
  geom_point(size = 3, aes(shape = Site), fill = "black") + 
  geom_smooth(method = "lm") + 
  ylab("Bacterial Production") +
  annotate("text", x = 2e+06, y=75, color = "black", fontface = "bold", size = 3.5,
           label = paste("Adj R2 =", round(summary(lm_HNA)$adj.r.squared, digits = 3), "\n", 
                         "p =", round(unname(summary(lm_HNA)$coefficients[,4][2]), digits = 8)))

lm_LNA <- lm(tot_bacprod ~ LNA.cells, data = muskegon)

LNA_vs_prod <- ggplot(muskegon, aes(x = LNA.cells, y = tot_bacprod)) + 
  geom_errorbarh(aes(xmin = LNA.cells - LNA.sd, xmax = LNA.cells + LNA.sd), color = "grey") + 
  geom_errorbar(aes(ymin = tot_bacprod - SD_tot_bacprod, max = tot_bacprod + SD_tot_bacprod), color = "grey") + 
  scale_shape_manual(values = c(21, 22, 23, 24)) + 
  geom_point(size = 3, aes(shape = Site), fill = "black") + 
  ylab("Bacterial Production") +
  geom_smooth(method = "lm", se = FALSE, linetype = "longdash", color = "red") + 
  annotate("text", x = 4e+06, y=75, color = "red", fontface = "bold", size = 3.5,
           label = paste("Adj R2 =", round(summary(lm_LNA)$adj.r.squared, digits = 3), "\n", 
                         "p =", round(unname(summary(lm_LNA)$coefficients[,4][2]), digits = 3)))

lm_total <- lm(tot_bacprod ~ Total.cells, data = muskegon)

Total_vs_prod <- ggplot(muskegon, aes(x = Total.cells, y = tot_bacprod)) + 
  geom_errorbarh(aes(xmin = Total.cells - Total.count.sd, xmax = Total.cells + Total.count.sd), color = "grey") + 
  geom_errorbar(aes(ymin = tot_bacprod - SD_tot_bacprod, max = tot_bacprod + SD_tot_bacprod), color = "grey") + 
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_point(size = 3, aes(shape = Site), fill = "black") + 
  ylab("Bacterial Production") +
  geom_smooth(method = "lm", se = FALSE, linetype = "longdash", color = "red") + 
  annotate("text", x = 7e+06, y=75, color = "red", fontface = "bold", size = 3.5,
           label = paste("Adj R2 =", round(summary(lm_total)$adj.r.squared, digits = 3), "\n", 
                         "p =", round(unname(summary(lm_total)$coefficients[,4][2]), digits = 3)))

# Put all three plots together into one 
totprod_plots <- plot_grid(HNA_vs_prod + theme(legend.position = c(0.84, 0.18),
                              legend.text = element_text(size = 10),
                              legend.title = element_text(size = 11, face = "bold")), 
          LNA_vs_prod + theme(legend.position = "none"), 
          Total_vs_prod + theme(legend.position = "none"),
          labels = c("A", "B", "C"), 
          ncol = 3)

# Save the plot to a file to call in the README
ggsave(plot = totprod_plots, 
       filename = "data/Chloroplasts_removed/HNA_vs_productivity.png", 
       width = 10, height = 4, units = "in", dpi = 500)


