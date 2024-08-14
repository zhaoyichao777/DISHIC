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

DISHIC <- function(file_path, feature_path, code_path, chr, cores, bin_size, limit_size, group_size) {
    options(scipen=999)
    setwd(code_path)
    source("data_process.R", local = TRUE)
    source("model_class.R", local = TRUE)
    source("solveRidgeRegression.R", local = TRUE)
    source("zinbInitialize.R", local = TRUE)
    source("zinbOptimize.R", local = TRUE)

    library(BiocParallel)
    BPPARAM <- MulticoreParam(workers = cores, progressbar = FALSE)

    zinb.loglik.single <- function(Y, mu, theta, logitPi) {
        logPnb <- dnbinom(Y, size = theta, mu = mu, log = TRUE)
        lognorm <- -log1pexp(logitPi)
        log_likelihood <- rowSums(logPnb * (Y > 0)) + rowSums(logPnb * (Y == 0) + log1pexp(logitPi - logPnb) * (Y == 0)) + rowSums(lognorm)
        return(log_likelihood)
    }

    log1pexp <- function(x, c0 = -37, c1 = 18, c2 = 33.3) {
        if (has.na <- any(ina <- is.na(x))) {
            y <- x
            x <- x[ok <- !ina]
        }
        r <- exp(x)
        if (any(i <- c0 < x & (i1 <- x <= c1))) {
            r[i] <- log1p(r[i])
        }
        if (any(i <- !i1 & (i2 <- x <= c2))) {
            r[i] <- x[i] + 1 / r[i]
        }
        if (any(i3 <- !i2)) {
            r[i3] <- x[i3]
        }
        if (has.na) {
            y[ok] <- r
            y
        } else {
            r
        }
    }
    ########################## data prepare #################################
    file_path1 <- file.path(file_path,file_name1)
    file_path2 <- file.path(file_path,file_name2)

    ######### combile all cells to one dataframe ########
    ori_file <- combine_cell(file_path1)
    fold_file <- combine_cell(file_path2)
    group <- c(rep(0, NCOL(ori_file)), rep(1, NCOL(fold_file)))

    file <- merge(ori_file, fold_file, by = "row.names", all = TRUE)
    rownames(file) <- file$Row.names
    file <- file[, -1]
    file[is.na(file)] <- 0
    file <- file[rowSums(file) != 0, ]
    ###
    ###
    if(!is.null(cell_feature) & NROW(cell_feature) != NCOL(file)){
        print("Cell feature size wrong!")
        stop("Stopping function due to incorrect cell feature size.")
    }

    # calculate number of groups
    n_rows <- nrow(file)
    n_groups <- ceiling(n_rows / group_size)
    groups <- rep(1:n_groups, each = group_size)[1:n_rows]
    file_list <- split(file, groups)
    number_of_groups <- length(file_list)
    print(paste0("number_of_groups:", number_of_groups))

    dir.create(paste0(file_path, "/result2/chr", chr), recursive = TRUE)

    results_list <- mapply(function(group_index, file_part) {
        feature_scHiCNorm <- get_scHiCNorm_feature(address = file.path(feature_path, paste0(chr, ".bin_features_", bin_size)))
        common_rows <- intersect(rownames(file_part), rownames(feature_scHiCNorm))
        file_part <- file_part[common_rows, ]
        feature_scHiCNorm <- feature_scHiCNorm[common_rows, ]

        Y <- as.matrix(file_part)
        null_model <- prepare_data(file = Y, bin_pair_feature=feature_scHiCNorm,cell_feature=cell_feature)
        group_model <- prepare_data(file = Y, bin_pair_feature = feature_scHiCNorm, cell_feature=cell_feature,group = group )

        # model initialize and optimize
        init_null_model <- zinbInitialize(Y, null_model, BPPARAM = BPPARAM)
        opt_null_model <- zinbOptimize(init_null_model, Y, maxiter = 15, stop.epsilon = 0.0001)
        init_group_model <- zinbInitialize(Y, group_model, BPPARAM = BPPARAM)
        opt_group_model <- zinbOptimize(init_group_model, Y, maxiter = 15, stop.epsilon = 0.0001)

        # calculate parameters
        mu_null <- exp(opt_null_model$A %*% opt_null_model$alpha_mu + t(opt_null_model$B %*% opt_null_model$beta_mu))
        theta_null <- exp(opt_null_model$zeta)
        logitPi_null <- opt_null_model$A %*% opt_null_model$alpha_pi + t(opt_null_model$B %*% opt_null_model$beta_pi)

        mu_group <- exp(opt_group_model$A %*% opt_group_model$alpha_mu + t(opt_group_model$B %*% opt_group_model$beta_mu))
        theta_group <- exp(opt_group_model$zeta)
        logitPi_group <- opt_group_model$A %*% opt_group_model$alpha_pi + t(opt_group_model$B %*% opt_group_model$beta_pi)

        # loglike ratio test
        prob_null <- zinb.loglik.single(Y, mu_null, theta_null, logitPi_null)
        prob_group <- zinb.loglik.single(Y, mu_group, theta_group, logitPi_group)
        lrt_statistic <- -2 * (prob_null - prob_group)
        p_values <- pchisq(lrt_statistic, df = 1, lower.tail = FALSE)
        results <- data.frame(
            mean_group1 = rowMeans(Y[, group == 0]),
            mean_group2 = rowMeans(Y[, group == 1]),
            mean_mu_null = rowMeans(mu_null),
            mean_mu_group = rowMeans(mu_group),
            prob_null,
            prob_group,
            lrt_statistic,
            p_value = p_values
        )
        rownames(results) <- rownames(file_part)
        return(results)
    }, seq_along(file_list), file_list, SIMPLIFY = FALSE)


    final_result <- do.call(rbind, results_list)
    write.table(final_result, paste0(file_path, "/result2/chr", chr, "/chr",chr,".txt"))
}

DISHIC(file_path, feature_path, code_path, chr, cores, bin_size, limit_size, group_size)
