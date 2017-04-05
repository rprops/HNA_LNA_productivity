# Purpose:  Plot the OTU abundances from Peter's Models
# Date:  April 4th, 2017 
# Author:  Marian L Schmidt 

#############################################################################################
################################# LOAD LIBRARIES ############################################
#############################################################################################
library(tidyverse)
library(cowplot)

# Set Phylum colors for plotting
colors <- c("#547980",  "#FFC543", "#A73E5C", "forestgreen", "#FF2151", "black", "#FF9E9D")


#############################################################################################
#################################### LOAD DATA ##############################################
#############################################################################################
# Read in the absolute abundance data 
absolute_otu <- read.table(file="data/Chloroplasts_removed/nochloro_absolute_otu.tsv", header = TRUE) # Absolute OTU abundance table 

# Read in the taxonomy data 
tax <- read.table(file="data/Chloroplasts_removed/nochloro_taxonomy_otu.tsv", header = TRUE) %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(Kingdom = Rank1,
         Phylum = Rank2, 
         Class = Rank3,
         Order = Rank4,
         Family = Rank5,
         Genus = Rank6,
         Species = Rank7,
         OTU = rowname) # Fix the Taxonomy


# Replace the phylum Proteobacteria with the class level
Phylum <- as.character(tax$Phylum)
Class <- as.character(tax$Class)

for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    Phylum[i] <- Class[i]
  } 
}

# Overwrite the Phylum level with the new phylum classification
tax$Phylum <- Phylum # Add the new phylum level data back to phy


# Read in the productivity and flow cytometry data 
productivity <- read.table(file = "data/Chloroplasts_removed/productivity_data.tsv", header = TRUE) # Metadata file


#############################################################################################
################################## RNADOMMIZED LASSO ########################################
########################### COMBINE DATA + PLOT OF 12 OTUS ##################################
#############################################################################################

# Vector of 12 OTUs pulled out by Peter's randomized lasso model
otus <- c("Otu000123", "Otu000027", "Otu000043","Otu000057", "Otu000176", "Otu000005", "Otu000048", "Otu000040", "Otu000058", "Otu000067", "Otu000029", "Otu000025")

# Put all the data together into one dataframe with only the important OTUs
data <- absolute_otu %>%
  dplyr::select(one_of(otus)) %>%  ### ONLY SELECT THE 12 OTUs FROM PETER'S MODEL
  tibble::rownames_to_column() %>%
  dplyr::rename(Sample_16S = rowname) %>%
  dplyr::left_join(productivity, by = "Sample_16S") %>%
  dplyr::filter(Lake == "Muskegon" & Depth == "Surface") %>%
  dplyr::select(-c(Platform, samples, Lake)) %>%
  mutate(Site = factor(Site, levels = c("MOT", "MDP", "MBR", "MIN"))) %>%
  gather("OTU", "Abs_Abund", 2:13) %>%
  dplyr::left_join(tax, by = "OTU") %>%
  filter(!is.na(Sample_16S)) %>%
  mutate(OTU_fraction_HNA = Abs_Abund/HNA.cells,
         OTU = factor(OTU, levels = OTU[order(Phylum)]))




# Plot the absolute abundance data 
plot_absabund <- ggplot(data, aes(x = reorder(OTU, Phylum), y = Abs_Abund, fill = Phylum, color = Phylum)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3.5e+5)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab("log10(Abundance)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# Plot the fraction of the HNA pool that each OTU is
plot_fracHNA <- ggplot(data, aes(x = reorder(OTU, Phylum), y = OTU_fraction_HNA, fill = Phylum, color = Phylum)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 0.35)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab("\n Abundance/HNA.cells ") +
  guides(fill = guide_legend(ncol=2),
         color = guide_legend(ncol=2)) +
  theme(legend.position = c(0.33, 0.85),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) 

# Combine the plots together 
plot_both <- plot_grid(plot_absabund, plot_fracHNA,
          labels = c("A", "B"),
          nrow = 2, ncol = 1)

# Save the plot to a .jpeg file 
#ggsave(plot_both, filename = "OTU_abundance_plot.jpeg", dpi = 400, width = 10, height = 7)



#############################################################################################
########################## CORRFILTERING RNADOMMIZED LASSO ##################################
########################### COMBINE DATA + PLOT OF 3 OTUS ###################################
#############################################################################################

# Vector of 3 OTUs pulled out by Peter's randomized lasso model WITH CORRFILTERING
three_otus <- c("Otu000060", "Otu000124", "Otu000357")

# Put all the data together into one dataframe with only the important OTUs
corrfilt_data <-  absolute_otu %>%
  dplyr::select(one_of(three_otus)) %>%     ### Use only 3 OTUs that Peter pulled out with corrfilt
  tibble::rownames_to_column() %>%
  dplyr::rename(Sample_16S = rowname) %>%  
  dplyr::left_join(productivity, by = "Sample_16S") %>%
  dplyr::filter(Lake == "Muskegon" & Depth == "Surface") %>%
  dplyr::select(-c(Platform, samples, Lake)) %>%
  mutate(Site = factor(Site, levels = c("MOT", "MDP", "MBR", "MIN"))) %>%
  gather("OTU", "Abs_Abund", 2:4) %>%     #### Gather only columns 2:4 
  dplyr::left_join(tax, by = "OTU") %>%
  filter(!is.na(Sample_16S)) %>%
  mutate(OTU_fraction_HNA = Abs_Abund/HNA.cells,
         OTU = factor(OTU, levels = OTU[order(Phylum)]))



# Plot the absolute abundance data 
plot_corrfilt_absabund <-   ggplot(corrfilt_data, aes(x = reorder(OTU, Phylum), y = Abs_Abund, fill = Phylum, color = Phylum)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2.5e+5)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab("log10(Abundance)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# Plot the fraction of the HNA pool that each OTU is
plot_corrfilt_fracHNA <-   ggplot(corrfilt_data, aes(x = reorder(OTU, Phylum), y = OTU_fraction_HNA, fill = Phylum, color = Phylum)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 0.065)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab("\n Abundance/HNA.cells ") +
  guides(fill = guide_legend(ncol=2),
         color = guide_legend(ncol=2)) +
  theme(legend.position = "none", #c(0.22, 0.85),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) 

# Combine the plots together 
plot_corrfilt_both <- plot_grid(plot_corrfilt_absabund, plot_corrfilt_fracHNA,
                                labels = c("C", "D"),
                                nrow = 2, ncol = 1)

plot_four <- plot_grid(plot_both, plot_corrfilt_both,
                       rel_widths = c(2.5,1),
                       nrow = 1, ncol = 2)
plot_four


# Save the plot to a .jpeg file 
ggsave(plot_four, filename = "OTU_abundance_plot_corrfilt.jpeg", dpi = 400, width = 10, height = 7)

