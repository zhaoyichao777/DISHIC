# --------------------------------------------------------------------------------
# This code is based on modifications made to the original zinbwave code.
# Original zinbwave code can be found at: https://github.com/biocore/zinbwave
# Modifications made to adapt to specific data structure and analysis requirements.
# --------------------------------------------------------------------------------
zinbInitialize <- function(Y, model, BPPARAM) {
    A <- model$A
    B <- model$B
    alpha_mu <- model$alpha_mu
    alpha_pi <- model$alpha_pi
    beta_mu <- model$beta_mu
    beta_pi <- model$beta_pi
    zeta <- model$zeta
    m <- model$m
    n <- model$n
    epsilon <- model$epsilon
    epsilon_alpha_mu <- model$epsilon_alpha_mu
    epsilon_beta_mu <- model$epsilon_beta_mu
    epsilon_alpha_pi <- model$epsilon_alpha_pi
    epsilon_beta_pi <- model$epsilon_beta_pi
    epsilon_zeta <- model$epsilon_zeta
    P <- Y > 0
    L <- matrix(NA, nrow = m, ncol = n)
    L[P] <- log(Y[P])
    Z <- 1 - P 

    ## initialize μ
    iter <- 0
    while (iter < 2) {
        # alpha_mu
        tBbeta_mu <- t(B %*% beta_mu)
        alpha_mu <- matrix(unlist(bplapply(seq(n), function(i) {
            solveRidgeRegression( 
                x = A[P[, i], , drop = FALSE], 
                y = L[P[, i], i] - tBbeta_mu[P[, i], i], 
                beta = alpha_mu[, i],
                epsilon = epsilon_alpha_mu,
                family = "gaussian"
            )
        }, BPPARAM = BPPARAM)), nrow = NCOL(A))
        # beta_mu
        Aalpha_mu <- A %*% alpha_mu # beta的一列一列去求（一个bin一个bin）
        beta_mu <- matrix(unlist(bplapply(seq(m), function(i) {
            solveRidgeRegression( # P[i, ]binpair i中，不为0的样本
                x = B[P[i, ], , drop = FALSE],
                y = L[i, P[i, ]] - Aalpha_mu[i, P[i, ]],
                beta = beta_mu[, i],
                epsilon = epsilon_beta_mu,
                family = "gaussian"
            )
        }, BPPARAM = BPPARAM)), nrow = NCOL(B))
        iter <- iter + 1
    }

    ## initialize pi
    iter <- 0

    while (iter < 2) {
        # alpha_pi
        tBbeta_pi <- t(B %*% beta_pi)
        alpha_pi <- matrix(unlist(bplapply(seq(n), function(i) {
            solveRidgeRegression(
                x = A,
                y = Z[, i],
                offset = tBbeta_pi[, i],
                beta = alpha_pi[, i],
                epsilon = epsilon_alpha_pi,
                family = "binomial"
            )
        }, BPPARAM = BPPARAM)), nrow = NCOL(A))
        # beta_pi
        Aalpha_pi <- A %*% alpha_pi
        beta_pi <- matrix(unlist(bplapply(seq(m), function(i) {
            solveRidgeRegression(
                x = B,
                y = Z[i, ],
                offset = Aalpha_pi[i, ],
                beta = beta_pi[, i],
                epsilon = epsilon_beta_pi,
                family = "binomial"
            )
        }, BPPARAM = BPPARAM)), nrow = NCOL(B))
        iter <- iter + 1
    }
    zeta <- rep(0, m)

    model <- create_model(
        A = A, B = B,
        alpha_mu = alpha_mu, alpha_pi = alpha_pi,
        beta_mu = beta_mu, beta_pi = beta_pi,
        zeta = zeta,
        m = m, n = n,
        epsilon = epsilon,
        epsilon_alpha_mu = epsilon_alpha_mu,
        epsilon_beta_mu = epsilon_beta_mu,
        epsilon_alpha_pi = epsilon_alpha_pi,
        epsilon_beta_pi = epsilon_beta_pi,
        epsilon_zeta = epsilon_zeta
    )
}
