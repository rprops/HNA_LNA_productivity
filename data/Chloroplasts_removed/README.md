## Information on the files in this repo:

1. `remove_chloroplasts_bylake.Rmd`: This R script creates 6 folders of data in the [ByLake_Filtering/](ByLake_Filtering/) organized by the two filtering methods [ByLake_Filtering/1in3/](ByLake_Filtering/1in3/) and [ByLake_Filtering/5in10/](ByLake_Filtering/5in10/). *NOTE: This filtering occurs **before** rarefying the data. removing chloroplast OTUs and rarefying.*

> Please see [remove_chloroplasts_bylake.md](remove_chloroplasts_bylake.md) for preparation of the data with some basic analysis
  
2. `Combine_producitivy_data.R`: R script for combining sample metadata with productivity data. Creates the following file:   
	- `productivity_data.tsv`: Sample meta-, flow cytometry, **and productivity data.**  
	- And the figure below:

![Relationships between Total Bacterial Production in ug C/L/Hr and  (A) High-Nucleic Acid Cells,  (B) Low-Nucleic Acid Cells, and (C) Total Cells](Cutoff_Analysis_Figures/HNA_vs_productivity.png)
