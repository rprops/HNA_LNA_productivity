# Date: November 10th, 2017 & March 13th, 2018
# Author: Marian L. Schmidt
# This folder was created to fix fasta header names and to run fasttree for phylogenetic diversity analysis. 


# Fix header names in fasta file 
## Copy and Rename original fasta file
cp total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta rep_16s_seqs.fasta

## Remove unnecessary information from the fasta headers 

cut -f1 -d "|" test.fasta | cut -f2 | sed 's/Otu/>Otu/'  > cut_test.fasta  


########### 
On March 13, 2018

### Prepare fasta files for FastTree

1. To plot a phylogeny with OTUs with HNA and LNA Lasso Scores to look for phylogenetic signal. The following code was ran in R with the vector of OTUs with RL scores of at least 0.15:

```
vector_of_otus <- as.vector(otu_scores_df$OTU)
write(vector_of_otus, file = "heatmap_figs/OTUnames_based_on_RLscores.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")
```

2. Next, I ran the following code in the shell: 

```
~/bbmap/filterbyname.sh in=cut_test.fasta names=OTUnames_based_on_RLscores.txt out=HNA_LNA_OTUs.fasta -include t
```

3. Since the above function puts an "N" where all the dashes are for alignment - let's replace them with "-" by running the following command.

```
sed 's/N/-/g' HNA_LNA_OTUs.fasta > HNA_LNA_OTUs_rmN.fasta
```

### Run FastTree

4. Log onto flux  

5. Make a new working directory, for me it was called `/scratch/lsa_fluxm/marschmi/HNA_LNA/March_2018`.

6. Load fasttree: `module load fasttree`

10. On flux run the `fasttree.pbs` script. It's main function is the following line of code:  


```
# Infer a tree with fasttree with the GTR+CAT 
### GTR: General time reversible model 

## Input file to fasttree = rmN_to_dash.fasta
## Output file to fasttree = newick_tree_rmN_to_dash.tre 

FastTree -gtr -nt -fastest  < HNA_LNA_OTUs_rmN.fasta > newick_tree_HNALNA_rmN.tre
```
