
Here is an example to run DISHIC on [GSE80006](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80006) chromosome 19.


```
# data folder root path
file_path <- "./data"
# group1 data subfolder name
file_name1 <- "chr19-fold"
# group2 data subfolder name
file_name2 <- "chr19-ori"
#DISHIC code folder path
code_path <- "./"
#scHiCNorm feature folder path
feature_path <- "./feature"

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
This is the detail information.

### 1. Input data
```
# data folder root path
file_path <- "./data"
# group1 data subfolder name
file_name1 <- "chr19-fold"
# group2 data subfolder name
file_name2 <- "chr19-ori"
```
The data folder should contain two subfolders representing the two control groups for analysis of variance: **chr19-fold** and **chr19-ori**.

#### File Structure
In each folder, the scHi-C data files for all samples are included. Each sample is stored in a separate file, and each file has three columns:

- **V1**: Chromosome bin1 start position (divided by the resolution)
- **V2**: Chromosome bin1 start position (divided by the resolution)
- **V3**: Interaction values between V1 and V2 bins

For example, in the file `GSM2109888_1_oocyte_NSN.200kb.txt`, the data is structured as follows:

| V1  | V2  | V3  |
| --- | --- | --- |
| 16  | 17  | 0   |
| 16  | 18  | 0   |
| 17  | 18  | 13  |
| 16  | 19  | 0   |
| 17  | 19  | 1   |
|……|……|……|


scHiCNorm features for other species and resolutions can be generated from their [website](http://dna.cs.miami.edu/scHiCNorm/)
