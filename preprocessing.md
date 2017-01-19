# Absolute versus relative abundances
Ruben Props and Marian L. Schmidt  
January 2017  


  
  


# Load libraries


```r
library(dplyr)
library(phyloseq)
```


# Load in Phyloseq Data
This phyloseq object will contain all samples from:  

- The inland lakes  
- Lake Michigan  
- Muskegon Lake 


```r
load("data/phyloseq.RData") # The phyloseq object will be physeq.otu

# Rename to something more intuitive
raw_physeq <- physeq.otu
```

