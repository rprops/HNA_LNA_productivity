-   [Load the necessary libraries and set colors](#load-the-necessary-libraries-and-set-colors)
-   [Load in the data](#load-in-the-data)
-   [Based on 1):](#based-on-1)
-   [Proportion of HNA Pool?](#proportion-of-hna-pool)
    -   [WHAT PROPORTION OF THE HNA POOL IS MADE UP BY THE 26 OTUs?](#what-proportion-of-the-hna-pool-is-made-up-by-the-26-otus)
-   [Sum OTUs vs HNA](#sum-otus-vs-hna)

### Load the necessary libraries and set colors

``` r
################################# LOAD LIBRARIES ############################################
library(tidyverse)
library(cowplot)

# Set Phylum colors for plotting
colors <- c(
  Actinobacteria = "skyblue",
  Alphaproteobacteria = "#547980",  
  Bacteroidetes = "#FFC543", 
  Betaproteobacteria = "#A73E5C", 
  Chlorobi = "#BEDB39",
  Cyanobacteria = "forestgreen", 
  Deltaproteobacteria = "#FF2151", 
  Gammaproteobacteria = "black", 
  Planctomycetes = "#FD7400",
  Proteobacteria_unclassified = "#019879",
  Verrucomicrobia = "#562258")
```

### Load in the data

``` r
#################################### LOAD DATA ##############################################

# Read in the absolute abundance data 
absolute_otu <- read.table(file="data/Chloroplasts_removed/nochloro_absolute_otu.tsv", header = TRUE) # Absolute OTU abundance table 
```

    ## Error in file(file, "rt"): cannot open the connection

``` r
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
```

    ## Error in file(file, "rt"): cannot open the connection

``` r
# Replace the phylum Proteobacteria with the class level
Phylum <- as.character(tax$Phylum)
```

    ## Error in eval(expr, envir, enclos): object 'tax' not found

``` r
Class <- as.character(tax$Class)
```

    ## Error in eval(expr, envir, enclos): object 'tax' not found

``` r
for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    Phylum[i] <- Class[i]
  } 
}
```

    ## Error in eval(expr, envir, enclos): object 'Phylum' not found

``` r
# Overwrite the Phylum level with the new phylum classification
tax$Phylum <- Phylum # Add the new phylum level data back to phy
```

    ## Error in eval(expr, envir, enclos): object 'Phylum' not found

``` r
# Read in the productivity and flow cytometry data 
productivity <- read.table(file = "data/Chloroplasts_removed/productivity_data.tsv", header = TRUE) # Metadata file
```

    ## Error in file(file, "rt"): cannot open the connection

Based on 1):
============

@prubbens performed an analysis where he removed three outliers from the productivity data on April 5th (*see `Final/OutliersRemoved/analysis_final_prod_outliersremoved.ipynb`*), which pulled out **19** otus, **3** OTUs which matched both pipelines.

``` r
# Vector of 3 OTUs pulled out by Peter's model in the above file 
OTUs_26 <- read.csv("final/OutliersRemoved/HNAscores_prod_outliersremoved_abun0.0075_Thr0.63.csv", header = FALSE) %>%
  dplyr::rename(OTU = V1,
                Corr = V2)
```

    ## Error in file(file, "rt"): cannot open the connection

``` r
otu_names_26 <- as.character(OTUs_26$OTU)
```

    ## Error in eval(expr, envir, enclos): object 'OTUs_26' not found

``` r
# What is the taxonomy of these 3 OTUs?
tax %>%
  dplyr::filter(OTU %in% otu_names_26)
```

    ## Error in eval(expr, envir, enclos): object 'tax' not found

``` r
# Put all the data together into one dataframe with only the important OTUs
AbsAbund_OTUs_26 <-  absolute_otu %>%
  dplyr::select(one_of(otu_names_26)) %>%     ### Use only 26 OTUs
  tibble::rownames_to_column() %>%
  dplyr::rename(Sample_16S = rowname) %>%  
  dplyr::left_join(productivity, by = "Sample_16S") %>%
  dplyr::filter(Lake == "Muskegon" & Depth == "Surface") %>%
  dplyr::select(-c(Platform, samples, Lake)) %>%
  mutate(Site = factor(Site, levels = c("MOT", "MDP", "MBR", "MIN"))) %>%
  gather("OTU", "Abs_Abund", 2:27) %>%     #### Gather only columns 2:12
  dplyr::left_join(tax, by = "OTU") %>%
  filter(!is.na(Sample_16S)) %>%
  mutate(OTU_fraction_HNA = Abs_Abund/HNA.cells,
         OTU = factor(OTU, levels = OTU[order(Phylum)])) %>%
  dplyr::filter(tot_bacprod < 90) # REMOVE OUTLIERS TO MATCH PETER'S ANALYSIS
```

    ## Error in eval(expr, envir, enclos): object 'absolute_otu' not found

``` r
# Plot the absolute abundance data 
ggplot(abs_abund_OTUs_26, aes(x = reorder(OTU, Phylum), y = Abs_Abund, fill = Phylum, color = Phylum)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, 2.5e+5)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab("log10(Abundance)") +
  ggtitle("26 HNA OTUs from HNA Prediction") + 
  guides(fill = guide_legend(ncol=2),
         color = guide_legend(ncol=2)) +
  theme(legend.position = c(0.33, 0.85),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

    ## Error in ggplot(abs_abund_OTUs_26, aes(x = reorder(OTU, Phylum), y = Abs_Abund, : object 'abs_abund_OTUs_26' not found

``` r
ggplot(abs_abund_OTUs_26, aes(x = reorder(OTU, Phylum), y = OTU_fraction_HNA, fill = Phylum, color = Phylum)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, 2.5e+5)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab("Fraction of HNA Cells") +
  ggtitle("26 HNA OTUs from HNA Prediction") + 
  guides(fill = guide_legend(ncol=2),
         color = guide_legend(ncol=2)) +
  theme(legend.position = c(0.33, 0.85),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

    ## Error in ggplot(abs_abund_OTUs_26, aes(x = reorder(OTU, Phylum), y = OTU_fraction_HNA, : object 'abs_abund_OTUs_26' not found

Proportion of HNA Pool?
=======================

### WHAT PROPORTION OF THE HNA POOL IS MADE UP BY THE 26 OTUs?

``` r
# Calculate the sum of the HNA pool, the max of the HNA pool, and the median/mean
frac_HNA_stats_AbsAbund_OTUs_26 <- AbsAbund_OTUs_26 %>%
  dplyr::select(Sample_16S, OTU, OTU_fraction_HNA, Abs_Abund, HNA.cells) %>%
  group_by(Sample_16S) %>%
  summarise(sum_fracHNA = sum(OTU_fraction_HNA), 
            max_fracHNA = max(OTU_fraction_HNA), 
            median_fracHNA = median(OTU_fraction_HNA), 
            mean_fracHNA = mean(OTU_fraction_HNA),
            sum_abs_abund = sum(Abs_Abund)) %>%
  mutate(All_Samples = "AllSamps")
```

    ## Error in eval(expr, envir, enclos): object 'AbsAbund_OTUs_26' not found

``` r
# Plot it
ggplot(frac_HNA_stats_AbsAbund_OTUs_26, aes(y = sum_fracHNA, x = All_Samples, 
                                            color = "All_Samples", fill = "All_Samples")) +
  geom_boxplot(alpha = 0.5) +   geom_point(size = 3, position = position_jitterdodge()) +
  ggtitle("26 HNA OTUs from HNA Prediction") + 
  xlab("All Samples") + scale_color_manual(values = "black") +
  scale_fill_manual(values = "black") +
  scale_y_continuous(expand = c(0,0),limits = c(0, 1.3), 
                     breaks = seq(0, 1.2, by = 0.1)) +
  ylab("\n Sum(Abundance/HNA.cells)") + xlab("Sample") +
  theme(legend.position = "none", axis.text.x = element_blank())
```

    ## Error in ggplot(frac_HNA_stats_AbsAbund_OTUs_26, aes(y = sum_fracHNA, : object 'frac_HNA_stats_AbsAbund_OTUs_26' not found

Sum OTUs vs HNA
===============

``` r
all_data <- inner_join(frac_HNA_stats_AbsAbund_OTUs_26, productivity, by = "Sample_16S") %>%
  mutate(pred_totHNA_counts = sum_fracHNA*HNA.cells)
```

    ## Error in inner_join(frac_HNA_stats_AbsAbund_OTUs_26, productivity, by = "Sample_16S"): object 'frac_HNA_stats_AbsAbund_OTUs_26' not found

``` r
dplyr::select(all_data, Sample_16S, HNA.cells, pred_totHNA_counts, sum_fracHNA)
```

    ## Error in select_(.data, .dots = lazyeval::lazy_dots(...)): object 'all_data' not found

``` r
ggplot(all_data, aes(x = sum_abs_abund, y= HNA.cells)) +
  geom_point(size = 3) + ylab("HNA Cell Count") + 
  xlab("Sum(Abs_Abund of 26 OTUs)") + 
  geom_abline(intercept = 0, slope = 1)
```

    ## Error in ggplot(all_data, aes(x = sum_abs_abund, y = HNA.cells)): object 'all_data' not found

``` r
plot_OTU25 <- ggplot(dplyr::filter(threeOTU_data, OTU == "Otu000025"), 
       aes(y = tot_bacprod,x = log10(Abs_Abund), fill = Phylum, color = Phylum)) +
  geom_point(size = 3) + ggtitle("Otu000025") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  xlab("log10(Abundance)") +
  geom_smooth(method = "lm", color = "#FFC543") +
  ylab("Total Production (ug C/L/hr)") + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

    ## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'threeOTU_data' not found

``` r
plot_OTU41 <- ggplot(dplyr::filter(threeOTU_data, OTU == "Otu000041"), 
       aes(y = tot_bacprod,x = log10(Abs_Abund), fill = Phylum, color = Phylum)) +
  geom_point(size = 3) +  ggtitle("Otu000041") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  xlab("log10(Abundance)") +
  geom_smooth(method = "lm", color = "#562258") +
  ylab("Total Production (ug C/L/hr)") + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

    ## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'threeOTU_data' not found

``` r
plot_OTU176 <- ggplot(dplyr::filter(threeOTU_data, OTU == "Otu000176"), 
       aes(y = tot_bacprod,x = log10(Abs_Abund), fill = Phylum, color = Phylum)) +
  geom_point(size = 3) + ggtitle("Otu000176") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  xlab("log10(Abundance)") +
  geom_smooth(method = "lm", color = "#FF2151") +
  ylab("Total Production (ug C/L/hr)") + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

    ## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'threeOTU_data' not found

``` r
plot_grid(plot_OTU25, plot_OTU176, plot_OTU41, nrow = 1, ncol = 3,
          labels = c("A", "B", "C"))
```

    ## Error in plot_grid(plot_OTU25, plot_OTU176, plot_OTU41, nrow = 1, ncol = 3, : object 'plot_OTU25' not found
