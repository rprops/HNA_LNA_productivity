## Purpose:  Helpful functions for Analysis
## Author: Marian L. Schmidt 
## Date: April 28th, 2017



## Short list and description of functions 
# 1. combine_OTU_data = combine absolute OTU abundances for selected OTUs with flow cy and productivity data
# 2. calc_fraction_HNA = For each sample calculate the sum/median/mean of the HNA pool
# 3. 
# 4. 
# 5. 






#################################################################################### 1
#################################################################################### 1
# Combine all many types of data together and make a dataframe 

# Inputs/Arguments: 
    # 1. absolute_otu_table = data.frame of absolute abudance OTU table with taxa as columns and samples as rows 
    # 2. otu_vector_names = vector of OTUs to be selected from the OTU table 
    # 3. productivity_fcm_data = data.frame of productivity and flow cytometry and other meta data

# Output: 
    # A long dataframe with flow cytometry, productivity, metadata, taxonomy, AND OTU absolute abundance data! 
    # Will include 2 columns that are key: 
              # 1. Abs_Abund = Absolute Abundance of that OTU
              # 2. OTU_fraction_HNA = What fraction of the total HNA pool is that OTU?

combine_OTU_data <- function(absolute_otu_table, otu_vector_names, productivity_fcm_data, taxonomy_table){
  
  # Select correct column number for step 9 below 
  select_cols <- length(otu_vector_names)+1
  
  absolute_otu_table %>%                                                   # 1. Start with the OTU table with absolute abundance counts 
    dplyr::select(one_of(otu_vector_names)) %>%                              # 2. Pull out only the relevant 26 OTUs from the above OTU table
    tibble::rownames_to_column() %>%                                         # 3. Change the sample names to be a column
    dplyr::rename(Sample_16S = rowname) %>%                                  # 4. Rename the sample names column to match up with other data frames
    dplyr::left_join(productivity_fcm_data, by = "Sample_16S") %>%           # 5. Join 26 OTU absolute abundance counts with rest of metadata 
    dplyr::select(-c(Platform, samples)) %>%                                 # 6. Remove unnecessary columns 
    gather("OTU", "Abs_Abund", 2:select_cols) %>%                            # 7. Gather only relevant columns, which represent OTU abs abundance counts, and put it in *long* format
    mutate(OTU_fraction_HNA = Abs_Abund/HNA.cells) %>%                       # 8. Calculate the fraction that each individual OTU takes up within the HNA pool for each sample
    dplyr::left_join(taxonomy_table, by = "OTU")                             # 9. Add the taxonomic information for each OTU
  
}

#################################################################################### 1
#################################################################################### 1






#################################################################################### 2
#################################################################################### 2
# For each sample calculate the sum/median/mean of the HNA pool
# Use the output from the first function in this file 

# Input/Argument:
    # 1. AbsAbund_OTUs = dataframe output from combine_OTU_data function 

calc_fraction_HNA <- function(AbsAbund_OTUs){
  
  AbsAbund_OTUs %>%                        # Take the dataframe from above
  dplyr::select(Sample_16S, OTU, OTU_fraction_HNA, Abs_Abund, HNA.cells) %>%   # Select only relevant columns
  group_by(Sample_16S) %>%                                                     # Make a data frame for each of the individual samples
  summarise(sum_fracHNA = sum(OTU_fraction_HNA),                               # Calculate the sum of the fraction of each OTU from the HNA pool (total HNA pool represented by each OTU)
            median_fracHNA = median(OTU_fraction_HNA),                         # Calculate the median of 2 lines above
            mean_fracHNA = mean(OTU_fraction_HNA),                             # Calculate the mean of 3 lines above
            sum_abs_abund = sum(Abs_Abund)) %>%                                # What's the sum of the absolute abundance of the 26 OTUs within each sample? 
  mutate(All_Samples = "AllSamps_26_OTUs")
}
#################################################################################### 2
#################################################################################### 2