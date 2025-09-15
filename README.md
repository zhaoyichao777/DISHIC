DISHIC method has been published on DOI:[10.1109/BIBM62325.2024.10821719](https://ieeexplore.ieee.org/document/10821719)

![image](https://github.com/user-attachments/assets/322e8b32-b4ee-41ce-b700-d17d7117a287)


Here is an example to run DISHIC on [GSE80006](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80006) chromosome 19.

```
source("DISHIC.R")
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
# cell-level covariate matrix, nrows is the cell number and ncols is feature number
cell_feature <- NULL

# chromosomes to be analyzed
chr <- 19
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
The data folder `file_path` should contain two subfolders `file_name1` and `file_name2` representing the two control groups for analysis of variance, such as **chr19-fold** and **chr19-ori** here.

#### File Structure
In each folder, the scHi-C data files for all samples are included. Each sample is stored in a separate file, and each file has three columns:

- **V1**: Chromosome bin1 start position (divided by the resolution)
- **V2**: Chromosome bin2 start position (divided by the resolution)
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

### 2. Code
```
code_path <- "./"
```
Folder containing the DISHIC code files:

- `DISHIC.R`: Main code for the method.
- `data_process.R`: Code for preprocessing input scHi-C data and features.
- `model_class.R`: DISHIC model.
- `solve_regression.R`: Regression function used to update model parameters.
- `zinb_initialize.R`: Function to initialize model parameters.
- `zinb_optimize.R`: Function to update model parameters.

### 3. Feature
```
#scHiCNorm feature folder path
feature_path <- "./feature"
# cell-level covariate matrix
cell_feature <- NULL
```
- The `feature_path` folder contains genome-level features. By default, DISHIC uses scHiCNorm features as genome-level features. Features for other species or resolutions not provided can be generated from their [website](http://dna.cs.miami.edu/scHiCNorm/).
- The `cell_feature` is a matrix of cell-level covariates, with the number of cells in the row and the columns as features. By default it is NULL, but users can define it and the method will detect it.

### 4. Other settings
```
# chromosomes to be analyzed
chr <- 19
# number of nodes for parallel computing
cores <- 40
# data binned resolution
bin_size <- 200000
# max genomic distance between analyzed bin-pairs
limit_size <- 10000000
 #if the file is too large, group the files into several groups with group_size bin-pairs and calculate them sequentially.
group_size <- 25000
 #whether to use gaussian process or not
GP <- False
```
- DISHIC analyzes intra-interactions in certain chromosome `chr`.
- DISHIC supports multi-core parallelism on CPU, `cores` represents the number of parallel cores.
- `bin_size` represents the resolution of input scHi-C data. For example, 200000 equals to resolution of 200kb.
- Due to the lack of analytical value in interactions between bin-pairs that are too far from the diagonal, the `limit_size` parameter restricts the maximum genomic distance between the analyzed genome pairs. For example, if we set limit_size=10000000, DISHIC will only analyze bin-pairs within 10Mbp.
- If the scHi-C data is too large to process, DISHIC can divide the files into several groups, each with `group_size` bin-pairs (rows), and calculate them sequentially.

### 5. Run
```
DISHIC(file_path, feature_path, code_path, chr, cores, bin_size, limit_size, group_size, GP)
```
<img width="1440" height="960" alt="image" src="https://github.com/user-attachments/assets/89038bc2-33cc-4ec8-9642-586b63538119" />
