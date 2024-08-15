
An example to run DISHIC on GSE80006 chromosome 19.

```
file_path <- "/zhaoy/DISHIC/data" # data folder root path
file_name1 <- "chr19-fold" # group1 data subfolder name
file_name2 <- "chr19-ori" # group2 data subfolder name
code_path <- "/zhaoy/DISHIC" #DISHIC code folder path
feature_path <- "/zhaoy/DISHIC/feature" #schiCNorm feature folder path
chr <- 19  # chromosomes to be analyzed
cell_feature <- NULL # cell-level covariate matrix, nrows is the cell number and ncols is feature number
cores <- 40 # number of nodes for parallel computing
bin_size <- 200000 # data binned resolution
limit_size <- 10000000 # max genomic distance between analyzed bin-pairs
group_size <- 25000 #if the file is too large, group the files into several groups with group_size bin-pairs and calculate them sequentially.

DISHIC(file_path, feature_path, code_path, chr, cores, bin_size, limit_size, group_size)
```
