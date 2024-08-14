combine_cell <- function(file_path) {
    combine_file <- data.frame()
    files <- list.files(file_path)
    for (file_name in files) {
        file <- read.table(file.path(file_path, file_name), header = T)
        colnames(file) <- c("V1", "V2", "V3")
        file <- file[abs(file$V2 - file$V1) <= limit_size / bin_size, ]
        if (nrow(combine_file) == 0) {
            combine_file <- file
        } else {
            if (nrow(combine_file) == nrow(file)) {
                if (all(combine_file[, c(1, 2)] == file[, c(1, 2)])) {
                    combine_file <- cbind(combine_file, file[, 3])
                } else {
                    combine_file <- merge(combine_file, file, by = c("V1", "V2"), all = TRUE)
                    print("merge")
                }
            } else {
                combine_file <- merge(combine_file, file, by = c("V1", "V2"), all = TRUE)
                print("merge due to different row counts")
            }
        }
        #print(combine_file)
        colnames(combine_file)[ncol(combine_file)] <- sub("\\.txt$", "", file_name)
        colnames(combine_file)[ncol(combine_file)] <- sub("\\.csv$", "", file_name)
    }
    combine_file[is.na(combine_file)]=0
    rownames(combine_file) <- paste0(combine_file$V1, "-", combine_file$V2)
    combine_file <- combine_file[, -c(1, 2)]
    combine_file <- combine_file[rowSums(combine_file) != 0, ]
    return(combine_file)
}


get_scHiCNorm_feature <- function(address) {
    dat_feature <- read.table(address)
    dat_feature$V2 <- dat_feature$V2 / bin_size
    dat_feature <- dat_feature[, -1]
    colnames(dat_feature) <- c("pos", "cut", "gc", "map")
    dat_feature <- dat_feature[-which(dat_feature$cut == 0 | dat_feature$gc == 0 | dat_feature$map == 0), ]

    dat_feature <- as.matrix(dat_feature)
    mat_cut <- as.matrix(log(dat_feature[, 2] %o% dat_feature[, 2]))
    mat_gc <- as.matrix(log(dat_feature[, 3] %o% dat_feature[, 3]))
    mat_map <- as.matrix(log(dat_feature[, 4] %o% dat_feature[, 4]))

    mat_cut <- (mat_cut - mean(c(mat_cut))) / sd(c(mat_cut))
    mat_gc <- (mat_gc - mean(c(mat_gc))) / sd(c(mat_gc))

    upper_tri_indices <- which(upper.tri(mat_map, diag = TRUE), arr.ind = TRUE)
    upper_tri_names <- data.frame(
        row = rownames(mat_map)[upper_tri_indices[, 1]],
        col = colnames(mat_map)[upper_tri_indices[, 2]]
    )

    feature <- data.frame(
        cut = mat_cut[upper_tri_indices],
        gc = mat_gc[upper_tri_indices],
        map = mat_map[upper_tri_indices]
    )
    rownames(feature) <- paste0(upper_tri_names$row, "-", upper_tri_names$col)
    return(feature)
}

prepare_data <- function(file, feature, bin_pair_feature, cell_feature = NULL,group) {
    m <- nrow(file) # bin-pair number
    n <- ncol(file) # cell number

    # A
    A <- as.matrix(cbind(rep(1, m), bin_pair_feature))
    
    #B
    if (is.null(cell_feature)) {
        if (missing(group)) {
            B <- as.matrix(rep(1, n))
        } else {
            B <- as.matrix(cbind(rep(1, n), group))
        }
    } else {
        if (missing(group)) {
            B <- as.matrix(cbind(rep(1, n), cell_feature))
        } else {
            B <- as.matrix(cbind(rep(1, n), cell_feature, group))
        }
    }
    model <- create_model(m = m, n = n, A = A, B = B)
    return(model)
}
