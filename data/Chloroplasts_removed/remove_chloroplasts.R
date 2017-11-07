# Script to remove chloroplasts from data

library("ggplot2")
library("dplyr")
library("phyloseq")


# Load in the raw data 
load("data/phyloseq.RData")
physeq.otu@otu_table <- t(physeq.otu@otu_table)

# Clean up the metadata file
metadata <- read.csv2("data/metadata.csv", stringsAsFactors = FALSE)
metadata <- metadata[metadata$Platform == "Accuri",]
metadata$Sample_fcm <- gsub(metadata$Sample_fcm, pattern="_rep.*", replacement="")
metadata <- do.call(rbind,by(metadata, INDICES = factor(metadata$Sample_fcm), 
                             FUN = unique))
metadata$Sample_16S[metadata$Lake=="Inland"] <- gsub(metadata$Sample_16S[metadata$Lake=="Inland"], pattern="-", replacement="")

# Import counts 
counts <- read.csv2("data/count_total_HNALNA.csv", stringsAsFactors = FALSE)
counts$samples[125:188] <- paste0(counts$samples[125:188],"-1")

# Merge metadata 
counts.total <- inner_join(counts, metadata, by=c("samples"="Sample_fcm"))

# replace unnecessary "-"
sample_names(physeq.otu) <- gsub(sample_names(physeq.otu),pattern="-",replacement="")

#remove cDNA
otu_table(physeq.otu) <- otu_table(physeq.otu)[grep(x=sample_names(physeq.otu),pattern="cD.",invert=TRUE), ]

# remove "renamed" tag
sample_names(physeq.otu) <- gsub(sample_names(physeq.otu),pattern=".renamed",replacement="")

### Select the metadata for which you have sequencing data (876 samples at this point)
temp1 <- counts.total 
temp2 <- data.frame(Sample=sample_names(physeq.otu))
temp3 <- semi_join(temp1,temp2,by=c("Sample_16S"="Sample"))
rownames(temp3) <- temp3$Sample
sample_data(physeq.otu) <- temp3

#################################################################################################
############## REMOVE CHLOROPLASTS
no_chloro_physeq <- subset_taxa(physeq.otu, Rank3 != "Chloroplast")
no_chloro_physeq_pruned <- prune_taxa(taxa_sums(no_chloro_physeq) > 0, no_chloro_physeq) 

# Calculate the sequencing depth of each sample
sums <- data.frame(Sample_16S=as.character(names(sample_sums(no_chloro_physeq_pruned))),
                   Total_Sequences=sample_sums(no_chloro_physeq_pruned), row.names=NULL)

sums_data <- left_join(data.frame(sample_data(no_chloro_physeq_pruned)), sums, by = "Sample_16S")
row.names(sums_data) <- sums_data$Sample_16S
sample_data(no_chloro_physeq_pruned) <- sums_data

# Remove samples with low sequencing depth 
no_chloro_physeq_pruned_seqs <- subset_samples(no_chloro_physeq_pruned, Total_Sequences > 4700)
# Now we have removed samples that have 4700 sequences or less 

# Remove OTUs with counts 0 and less
no_chloro_physeq_pruned_seqs_rm0 <- prune_taxa(taxa_sums(no_chloro_physeq_pruned_seqs) > 0, no_chloro_physeq_pruned_seqs) 
  
#################################################################################################
############## PREVALENCE FILTERING
## On November 6th we decided to run the analysis 2 ways:
     # 1. Liberal Cutoff: Remove OTUs that have only one sequence in three samples 
     # 2. Conservative Cutoff: Remove OTUs that have less than 5 sequences in 10% of samples
### WILL HAVE 2 PHYLOSEQ OBJECTS FROM NOW ON!

# Remove OTUs that have only one sequence in three samples 
filter_val_1 <- 3/173
no_chloro_physeq_pruned_seqs_rm_1in3 <- filter_taxa(no_chloro_physeq_pruned_seqs_rm0, function(x) sum(x > 1) > (filter_val_1*length(x)), TRUE)


# Remove OTUs that have less than 5 sequences in 10% of samples
filter_val_2 <- 0.10
no_chloro_physeq_pruned_seqs_rm5_in_10percent <- filter_taxa(no_chloro_physeq_pruned_seqs_rm0, function(x) sum(x > 5) > (filter_val_2*length(x)), TRUE)




#################################################################################################
############## RELATIVE ABUNDANCES VIA RAREFY-ING
## RAREFY EVEN DEPTH 
min_lib_rm1in3 <- min(sample_sums(no_chloro_physeq_pruned_seqs_rm_1in3)) - 1
min_lib_rm_5in10percent <- min(sample_sums(no_chloro_physeq_pruned_seqs_rm5_in_10percent)) - 1
min_lib_rm1in3 # 4726 is the smallest library size
min_lib_rm_5in10percent # 4482 is the smallest library size 


#################################################################################################
#################################################################################################
################################################################################################# YOU ARE HERE MARIAN 
#################################################################################################
#################################################################################################

# Set the seed for the randomization of rarefy-ing
set.seed(777)

# Rarefy to the even depth 
rare_nochloro_rm0 <- rarefy_even_depth(no_chloro_physeq_pruned_seqs_rm0, sample.size = min_lib_rm0,
                                        verbose = FALSE, replace = TRUE)

# Sanity Check - Samples all have sample counts of 4760
apply(otu_table(rare_nochloro_rm0), 1, sum)


### New phyloseq objects with relative abundances
rare_nochloro_rm0_rel <- transform_sample_counts(rare_nochloro_rm0, function(x) x/sum(x))

# Sanity Check
apply(otu_table(rare_nochloro_rm0_rel), 1, sum)

### Replacing NA by "Unclassified"
tax_table(rare_nochloro_rm0_rel)[is.na(tax_table(rare_nochloro_rm0_rel))] <- "Unclassified"

# Final Relative Abundance Phyloseq Object 
rare_nochloro_rm0_rel


#################################################################################################
############## ABSOLUTE ABUNDANCES
### New phyloseq object with absolute abundances in cells/mL
rare_nochloro_rm0_abs <- rare_nochloro_rm0_rel
otu_table(rare_nochloro_rm0_abs) <- 1000*otu_table(rare_nochloro_rm0_rel)*sample_data(rare_nochloro_rm0_rel)$Total.cells

### to cells/mL
sample_data(rare_nochloro_rm0_rel)[,c(3:8)] <- sample_data(rare_nochloro_rm0_rel)[,c(3:8)]*1000
sample_data(rare_nochloro_rm0_abs)[,c(3:8)] <- sample_data(rare_nochloro_rm0_abs)[,c(3:8)]*1000

# Sanity Check
apply(otu_table(rare_nochloro_rm0_abs), 1, sum)


### Write the tsv files
write.table(otu_table(rare_nochloro_rm0_abs), file="data/Chloroplasts_removed/nochloro_absolute_otu.tsv", row.names=TRUE)
write.table(otu_table(rare_nochloro_rm0_rel), file="data/Chloroplasts_removed/nochloro_relative_otu.tsv", row.names=TRUE)
write.table(sample_data(rare_nochloro_rm0_rel)[,-1], file="data/Chloroplasts_removed/nochloro_HNA_LNA.tsv", row.names=TRUE)
write.table(tax_table(rare_nochloro_rm0_rel), file="data/Chloroplasts_removed/nochloro_taxonomy_otu.tsv", row.names=TRUE)
