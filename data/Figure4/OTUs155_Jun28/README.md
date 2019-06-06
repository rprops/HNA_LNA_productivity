# Phylogenetic Tree Construction with Fasttree
## Marian Schmidt
## June 28th, 2018


## For All OTUs that pass the RL score threshold: 155 OTUs

1. Create this repo in the fasttree folder.

2. To plot a phylogeny with OTUs with HNA and LNA Lasso Scores to look for phylogenetic signal. The following code was ran in R with the vector of OTUs with RL scores of at least 0.15:

```
vector_of_otus <- as.vector(otu_scores_df$OTU)
write(vector_of_otus, file = "data/fasttree/Figure5/OTUnames_RLscores_155.txt",
      ncolumns = 1, append = FALSE, sep = "\n")
```

3. Next, I ran the following code in the shell: 

```
~/bbmap/filterbyname.sh in=../cut_test.fasta names=OTUnames_RLscores_155.txt out=OTUs155.fasta -include t
```

4. Since the above function puts an "N" where all the dashes are for alignment - let's replace them with "-" by running the following command.

```
sed 's/N/-/g' OTUs155.fasta > OTUs155_rmN.fasta
```



## For the Top 10 OTUs from each system: 41 OTUs

5. Make file with OTU names in R:

```
top10_vector_of_otus <- as.vector(top10_otu_scores_df$OTU)
# Write the file 
write(top10_vector_of_otus, file = "data/fasttree/Figure5/OTUnames_RLscores_top10.txt",
      ncolumns = 1, append = FALSE, sep = "\n")
```

6. Subset the sequences from the top 10 OTUs from the larger fasta file with filterbyname.sh from bbtools.

**NOTE:** If a file already exists - the command will not run because of `overwrite=FALSE`

```
~/bbmap/filterbyname.sh in=../cut_test.fasta names=OTUnames_RLscores_top10.txt out=OTUs_top10.fasta -include t
```

7. Remove the N's for -.

```
sed 's/N/-/g' OTUs_top10.fasta > OTUs_top10_rmN.fasta
```

## Copy Files to flux

8. I will use the transfer host through U of M (requires password). Working on my own personal computer 

```
# All 155 RLscore threshold OTUs
scp OTUs155_rmN.fasta flux-xfer.arc-ts.umich.edu:/scratch/lsa_fluxm/marschmi/HNA_LNA/JUNE_2018

# Top 10 OTUs from each syste
scp OTUs_top10_rmN.fasta flux-xfer.arc-ts.umich.edu:/scratch/lsa_fluxm/marschmi/HNA_LNA/JUNE_2018
```

## Run FastTree

9. Run the fasttree pbs scripts with `qsub fasttree_top10.pbs` and `qsub fasttree_OTUs_155.pbs`. 

10. After running fasttree, transfer those files back local computer. I ran code like (password required):

```
scp flux-xfer.arc-ts.umich.edu:/scratch/lsa_fluxm/marschmi/HNA_LNA/JUNE_2018/FastTree-Top10* .
```




