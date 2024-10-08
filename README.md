
An example to run DISHIC on [GSE80006](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80006) chromosome 19.

scHiCNorm features for other species and resolutions can be generated from their [website](http://dna.cs.miami.edu/scHiCNorm/)

```
# data folder root path
file_path <- "/zhaoy/DISHIC/data"
# group1 data subfolder name
file_name1 <- "chr19-fold"
# group2 data subfolder name
file_name2 <- "chr19-ori"
#DISHIC code folder path
code_path <- "/zhaoy/DISHIC"
#scHiCNorm feature folder path
feature_path <- "/zhaoy/DISHIC/feature"

# chromosomes to be analyzed
chr <- 19
# cell-level covariate matrix, nrows is the cell number and ncols is feature number
cell_feature <- NULL

# number of nodes for parallel computing
cores <- 40
# data binned resolution
bin_size <- 200000
# max genomic distance between analyzed bin-pairs
limit_size <- 10000000
 #if the file is too large, group the files into several groups with group_size bin-pairs and calculate them sequentially.
group_size <- 25000

DISHIC(file_path, feature_path, code_path, chr, cores, bin_size, limit_size, group_size)
```
