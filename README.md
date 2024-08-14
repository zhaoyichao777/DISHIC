# DISHIC
# library(optparse)
# option_list <- list(
#     make_option(c("--file_path"), type = "character", help = "Path of data files"),
#     make_option(c("--feature_path"), type = "character", help = "Path of scHiCNorm feature folder"),
#     make_option(c("--code_path"), type = "character", help = "Path of the scDI code"),
#     make_option(c("--chr"), type = "character", help = "chromosome"),
#     make_option(c("--cores"), type = "integer", help = "core numbers"),
#     make_option(c("--binsize"), type = "integer", help = "resulotion bin size"),
#     make_option(c("--limitsize"), type = "integer", help = "limit genomic distance")
# )

```
file_path <- "/home/zhaoy/my-scHiCDiff/code/pipeline/DISHIC-git/data"
file_name1 <- "chr19-fold"
file_name2 <- "chr19-ori"
code_path <- "/home/zhaoy/my-scHiCDiff/code/pipeline/DISHIC-git"
feature_path <- "/home/zhaoy/my-scHiCDiff/feature"
chr <- 19
cell_feature <- NULL
cores <- 40
bin_size <- 200000
limit_size <- 10000000
group_size <- 25000
```

DISHIC <- function(file_path, feature_path, code_path, chr, cores, bin_size, limit_size, group_size)