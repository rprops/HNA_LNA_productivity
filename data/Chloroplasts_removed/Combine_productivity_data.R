# March 4th, 2017 
# Marian L Schmidt 

# Load libraries
library(dplyr)

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
final_data <- left_join(raw_data, combined_data, by = "norep_filter_name")

# Write out the file 
write.table(final_data, file="data/Chloroplasts_removed/productivity_data.tsv", row.names=FALSE)



