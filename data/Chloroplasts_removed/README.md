## Information on the files in this repo:

1. `remove_chloroplasts.R`:  R script for removing chloroplast OTUs and rarefying. Creates the following files: 
  
	- `nochloro_HNA_LNA.tsv`:  Sample meta- and flow cytometry data. 
	- `nochloro_relative_otu.tsv`: Relative abundance OTU table.    
	- `nochloro_taxonomy_otu.tsv`: OTU taxonomy file.  
	- `nochloro_absolute_otu.tsv`: Absolute abundance OTU table calculated from flow cytometry data. 

2. `Combine_producitivy_data.R`: R script for combining sample metadata with productivity data. Creates the following file:   
	- `productivity_data.tsv`: Sample meta-, flow cytometry, **and productivity data.**  


