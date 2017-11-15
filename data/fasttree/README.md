# Date: November 10th, 2017
# Author: Marian L. Schmidt
# This folder was created to fix fasta header names and to run fasttree for phylogenetic diversity analysis. 


# Fix header names in fasta file 
## Copy and Rename original fasta file
cp total.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta rep_16s_seqs.fasta

## Remove unnecessary information from the fasta headers 

cut -f1 -d "|" test.fasta | cut -f2 | sed 's/Otu/>Otu/'  > cut_test.fasta  

