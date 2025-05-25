library(Matrix)
library(dplyr)
library(optparse)
library(parallel)

# Merge data into pseudo bulk data
mergeData <- function(filePaths, noiseProp, chromosome, nCells) {
    print(paste("heree", chromosome))

    print("start merge data")
    file_bin <- paste0(filePaths, "/", chromosome)
    fileNames <- dir(file_bin)
    fileNames <- sample(fileNames, nCells)
    fileLists <- list()
    hicMerge <- NULL
    count <- numeric(length(fileNames))
    chrLengthMax <- 0
    # combine
    for (i in 1:length(fileNames)) {
        data <- read.table(paste(file_bin, fileNames[i], sep = "/")) %>%
            mutate(V1_new = pmin(V1, V2), V2_new = pmax(V1, V2)) %>%
            group_by(V1_new, V2_new) %>%
            summarise(V3 = sum(V3), .groups = "drop") %>%
            ungroup() %>%
            select(V1 = V1_new, V2 = V2_new, V3)
        # data <- data[(data$V2 - data$V1) * opt$bin_size <= opt$limit_size, ]
        fileLists[[i]] <- as.matrix(data)
        count[i] <- sum(fileLists[[i]][, 3])
        chrLengthCurrent <- max(fileLists[[i]][, 1:2])
        chrLengthMax <- max(chrLengthMax, chrLengthCurrent)
        hicSparse <- sparseMatrix(
            i = fileLists[[i]][, 1], j = fileLists[[i]][, 2],
            x = fileLists[[i]][, 3], dims = c(chrLengthMax, chrLengthMax)
        )
        if (is.null(hicMerge)) {
            hicMerge <- hicSparse
        } else {
            if (dim(hicMerge)[1] < chrLengthMax) {
                # 扩展行
                hicMerge <- rbind(hicMerge, Matrix(0, nrow = chrLengthMax - dim(hicMerge)[1], ncol = dim(hicMerge)[1], sparse = TRUE))
                hicMerge <- cbind(hicMerge, Matrix(0, nrow = dim(hicMerge)[1], ncol = chrLengthMax - dim(hicMerge)[2], sparse = TRUE))
            }
            hicMerge <- hicMerge + hicSparse
        }
    }
    hicMerge[lower.tri(hicMerge)] <- 0

    # p.prop
    pProp <- count / sum(hicMerge)

    allPositions <- do.call(rbind, lapply(fileLists, function(x) x[, 1:2]))
    allPositions <- unique(allPositions)
    propMatrix <- matrix(0, nrow = nrow(allPositions), ncol = length(fileLists))
    colnames(propMatrix) <- fileNames
    rownames(propMatrix) <- paste0(allPositions[, 1], "-", allPositions[, 2])
    for (i in 1:length(fileLists)) {
        indexVector <- paste0(fileLists[[i]][, 1], "-", fileLists[[i]][, 2])
        propMatrix[match(indexVector, rownames(propMatrix)), i] <- fileLists[[i]][, 3]
    }
    back <- propMatrix
    # calculate p.prop
    for (i in 1:nrow(propMatrix)) {
        propMatrix[i, ] <- noiseProp * propMatrix[i, ] / sum(propMatrix[i, ]) + (1 - noiseProp) * pProp
    }
    return(list(hicMerge, propMatrix))
}

# Generation of pseudo-bulk HiC matrices with predefined differential interactions
generatePseudoBulkHiCMatrices <- function(mergedData, foldChange, pDiff) {
    print("start generate pseudo data")
    summaryData <- summary(mergedData[[1]])
    positionUsable <- which(summaryData$i != summaryData$j & abs(summaryData$j - summaryData$i) <= opt$limit_size / opt$bin_size)
    # set.seed(123)
    numDiffPos <- floor(pDiff * length(positionUsable))
    diffPos <- sample(positionUsable, numDiffPos)

    numUp <- numDiffPos %/% 2 # Ensure equal number of up and down
    numDown <- numDiffPos - numUp
    diffPosUp <- diffPos[1:numUp]
    diffPosDown <- diffPos[(numUp + 1):numDiffPos]

    summaryData$x2 <- summaryData$x
    summaryData[diffPosUp, "x2"] <- foldChange * summaryData[diffPosUp, "x2"]
    summaryData[diffPosDown, "x2"] <- summaryData[diffPosDown, "x2"] / foldChange
    summaryData$fold <- 0
    summaryData[diffPos, ]$fold <- 1
    summaryData$type <- "n"
    summaryData[diffPosUp, ]$type <- "up"
    summaryData[diffPosDown, ]$type <- "down"
    hicOri <- summaryData[, c(1, 2, 3)]
    hicFold <- summaryData[, c(1, 2, 4)]
    stateMatrix <- summaryData[, c(1, 2, 5, 6)]
    return(list(hicOri, hicFold, stateMatrix))
}

# Downsampling for simulated single-cell samples
downsampleToSingleCell <- function(pseudoBulkMatrices, mergedData, chromosome, foldChange, pDiff, noiseProp) {
    print("start downsample data")

    hicOri <- pseudoBulkMatrices[[1]]
    hicFold <- pseudoBulkMatrices[[2]]

    adjusted_sizes_ori <- round(hicOri$x * runif(nrow(hicOri), min = 0.95, max = 1.05))
    adjusted_sizes_fold <- round(hicFold$x * runif(nrow(hicFold), min = 0.95, max = 1.05))

    stateMatrix <- pseudoBulkMatrices[[3]]
    pProp <- mergedData[[2]][paste0(hicOri$i, "-", hicOri$j), ]
    datOri <- array(0, dim = c(ncol(pProp), nrow(hicOri), 3))
    for (i in 1:nrow(hicOri)) {
        datOri[, i, 1] <- rep(hicOri$i[i], ncol(pProp))
        datOri[, i, 2] <- rep(hicOri$j[i], ncol(pProp))
        datOri[, i, 3] <- rmultinom(1, adjusted_sizes_ori[i], prob = pProp[i, ])
    }
    datFold <- array(0, dim = c(ncol(pProp), nrow(hicFold), 3))
    for (i in 1:nrow(hicFold)) {
        datFold[, i, 1] <- rep(hicFold$i[i], ncol(pProp))
        datFold[, i, 2] <- rep(hicFold$j[i], ncol(pProp))
        datFold[, i, 3] <- rmultinom(1, adjusted_sizes_fold[i], prob = pProp[i, ])
    }
    for (i in 1:ncol(pProp)) {
        name <- colnames(pProp)[i]
        file_write <- paste0(opt$output_dir, "/random-result-fold", foldChange, "-pdiff", pDiff, "-noise", noiseProp, "-10Mb-rep", opt$index)
        file_ori <- paste0(file_write, "/", chromosome, "-ori")
        dir.create(file_ori, recursive = T)
        print(file_ori)
        fileName_ori <- paste0(file_ori, "/", name)
        write.table(datOri[i, , ], fileName_ori, sep = "\t", row.names = F, quote = F)

        file_fold <- paste0(file_write, "/", chromosome, "-fold")
        dir.create(file_fold, recursive = T)
        fileName_fold <- paste0(file_fold, "/", name)
        # print(fileName_fold)
        write.table(datFold[i, , ], fileName_fold, sep = "\t", row.names = F, quote = F)
    }
    write.table(stateMatrix, paste0(file_write, "/stateMatrix_", chromosome, ".txt"), sep = "\t", row.names = F, quote = F)
}

# main
scHiCDiffSim <- function(filePaths, foldChange, pDiff, noiseProp, chromosome) {
    mergedData <- mergeData(filePaths, noiseProp, chromosome, opt$ncells)
    pseudoBulkMatrices <- generatePseudoBulkHiCMatrices(mergedData, foldChange, pDiff)
    singleCellSimulations <- downsampleToSingleCell(pseudoBulkMatrices, mergedData, chromosome, foldChange, pDiff, noiseProp)
}


option_list <- list(
    make_option(c("--filePaths"), type = "character", help = "Path to data files"),
    make_option(c("--noiseProp"), type = "numeric", help = "Noise proportion"),
    make_option(c("--foldChange"), type = "numeric", help = "Fold change"),
    make_option(c("--pDiff"), type = "numeric", help = "Probability difference"),
    make_option(c("--index"), type = "integer"),
    make_option(c("--bin_size"), type = "integer"),
    make_option(c("--limit_size"), type = "integer"),
    make_option(c("--cores"), type = "integer"),
    make_option(c("--ncells"), type = "integer"),
    make_option(c("--output_dir"), type = "integer"),
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

chro_bins <- list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")

mclapply(chro_bins, function(chro) {
    scHiCDiffSim(
        filePaths = opt$filePaths,
        foldChange = opt$foldChange,
        pDiff = opt$pDiff,
        noiseProp = opt$noiseProp,
        chromosome = chro
    )
}, mc.cores = opt$cores)
