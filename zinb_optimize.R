# -------------------------------------------------------------------------------- This code is based on modifications made to the original zinbwave code.  Original
# zinbwave code can be found at: https://github.com/biocore/zinbwave Modifications made to adapt to specific data structure and analysis requirements.
# --------------------------------------------------------------------------------
zinbOptimize <- function(model, file, maxiter = 25, stop.epsilon = 0.001) {
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

    epsilonAlpha <- c(epsilon_alpha_mu, epsilon_alpha_pi)
    epsilonBeta <- c(epsilon_beta_mu, epsilon_beta_pi)
    Y <- as.matrix(file)

    total.lik <- rep(NA, maxiter)

    iter <- 0

    while (iter < maxiter) {
        gc()
        print(paste("iter:", iter))
        mu <- exp(A %*% alpha_mu + t(B %*% beta_mu))
        logitPi <- A %*% alpha_pi + t(B %*% beta_pi)
        theta <- exp(zeta)
        loglik <- zinb.loglik(Y, mu, rep(theta, rep(n, m)), logitPi)
        penalty <- sum(epsilon_alpha_mu * (alpha_mu)^2)/2 + sum(epsilon_alpha_pi * (alpha_pi)^2)/2 + sum(epsilon_beta_mu * (beta_mu)^2)/2 + sum(epsilon_beta_pi * (beta_pi)^2)/2 +
            epsilon_zeta * var(zeta)/2
        loglik <- loglik - penalty
        total.lik[iter] <- loglik

        if (iter > 1) {
            if (abs((total.lik[iter] - total.lik[iter - 1])/total.lik[iter - 1]) < stop.epsilon) {
                print("break")
                break
            }
        }
        ########## zeta ############
        print("zeta")
        zeta <- optimizeDispersion(m, mu, logitPi, epsilon_zeta, Y)
        penalty <- sum(epsilon_alpha_mu * (alpha_mu)^2)/2 + sum(epsilon_alpha_pi * (alpha_pi)^2)/2 + sum(epsilon_beta_mu * (beta_mu)^2)/2 + sum(epsilon_beta_pi * (beta_pi)^2)/2 +
            epsilon_zeta * var(zeta)/2
        print(zinb.loglik(Y, mu, rep(exp(zeta), rep(n, m)), logitPi) - penalty)
        gc()
        ########## alpha ############
        print("alpha")
        estimate_alpha <- matrix(unlist(bplapply(seq(n), function(i) {
            optimright_fun(alpha_mu[, i], alpha_pi[, i], Y[, i], A = A, B = B[i, ], beta_mu, beta_pi, zeta, epsilonAlpha)
        }, BPPARAM = BPPARAM)), nrow = NROW(alpha_mu) + NROW(alpha_pi))

        ind <- 1
        alpha_mu <- estimate_alpha[ind:(ind + NROW(alpha_mu) - 1), , drop = FALSE]
        ind <- ind + NROW(alpha_mu)
        alpha_pi <- estimate_alpha[ind:(ind + NROW(alpha_pi) - 1), , drop = FALSE]


        itermu <- exp(A %*% alpha_mu + t(B %*% beta_mu))
        iterlP <- A %*% alpha_pi + t(B %*% beta_pi)
        penalty <- sum(epsilon_alpha_mu * (alpha_mu)^2)/2 + sum(epsilon_alpha_pi * (alpha_pi)^2)/2 + sum(epsilon_beta_mu * (beta_mu)^2)/2 + sum(epsilon_beta_pi * (beta_pi)^2)/2 +
            epsilon_zeta * var(zeta)/2
        print(zinb.loglik(Y, itermu, rep(exp(zeta), rep(n, m)), iterlP) - penalty)
        gc()
        ########## beta ############
        print("beta")
        estimate_beta <- matrix(unlist(bplapply(seq(m), function(i) {
            optimleft_fun(beta_mu[, i], beta_pi[, i], Y[i, ], A = A[i, ], B = B, alpha_mu, alpha_pi, zeta[i], n, epsilonBeta)
        }, BPPARAM = BPPARAM)), nrow = NROW(beta_mu) + NROW(beta_pi))
        ind <- 1
        beta_mu <- estimate_beta[ind:(ind + NROW(beta_mu) - 1), , drop = FALSE]
        ind <- ind + NROW(beta_mu)
        beta_pi <- estimate_beta[ind:(ind + NROW(beta_pi) - 1), , drop = FALSE]

        itermu <- exp(A %*% alpha_mu + t(B %*% beta_mu))
        iterlP <- A %*% alpha_pi + t(B %*% beta_pi)
        penalty <- sum(epsilon_alpha_mu * (alpha_mu)^2)/2 + sum(epsilon_alpha_pi * (alpha_pi)^2)/2 + sum(epsilon_beta_mu * (beta_mu)^2)/2 + sum(epsilon_beta_pi * (beta_pi)^2)/2 +
            epsilon_zeta * var(zeta)/2
        print(zinb.loglik(Y, itermu, rep(exp(zeta), rep(n, m)), iterlP) - penalty)
        iter <- iter + 1
    }
    opt_model <- create_model(A = A, B = B, alpha_mu = alpha_mu, alpha_pi = alpha_pi, beta_mu = beta_mu, beta_pi = beta_pi, zeta = zeta, m = m, n = n, epsilon = epsilon,
        epsilon_alpha_mu = epsilon_alpha_mu, epsilon_beta_mu = epsilon_beta_mu, epsilon_alpha_pi = epsilon_alpha_pi, epsilon_beta_pi = epsilon_beta_pi, epsilon_zeta = epsilon_zeta)
    return(opt_model)
}



optimright_fun <- function(alpha_mu, alpha_pi, Y, A, B, beta_mu, beta_pi, zeta, epsilonAlpha) {
    par <- c(alpha_mu, alpha_pi)
    par_length <- length(par)
    lower_bounds <- rep(-10, par_length)
    upper_bounds <- rep(10, par_length)
    optim(fn = zinb.loglik.regression, gr = zinb.loglik.regression.gradient, par = par, Y = Y, feature.mu = A, offset_mu = t(B %*% beta_mu), feature.pi = A, offset_pi = t(B %*%
        beta_pi), C.theta = zeta, epsilon = epsilonAlpha, control = list(fnscale = -1, trace = 0), method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)$par
}


optimleft_fun <- function(beta_mu, beta_pi, Y, A, B, alpha_mu, alpha_pi, zeta, n, epsilonBeta) {
    par <- c(beta_mu, beta_pi)
    par_length <- length(par)
    lower_bounds <- rep(-10, par_length)
    upper_bounds <- rep(10, par_length)
    optim(fn = zinb.loglik.regression, gr = zinb.loglik.regression.gradient, par = par, Y = t(Y), feature.mu = B, offset_mu = t(A %*% alpha_mu), feature.pi = B, offset_pi = t(A %*%
        alpha_pi), C.theta = matrix(zeta, nrow = n, ncol = 1), epsilon = epsilonBeta, control = list(fnscale = -1, trace = 0), method = "L-BFGS-B", lower = lower_bounds,
        upper = upper_bounds)$par
}


zinb.regression.parseModel <- function(para, feature.mu, offset_mu, feature.pi, offset_pi) {
    logMu <- offset_mu
    logitPi <- offset_pi
    dim.para <- rep(0, 2)
    start.para <- rep(NA, 2)
    i <- 0

    j <- ncol(feature.mu)
    if (j > 0) {
        logMu <- logMu + feature.mu %*% para[(i + 1):(i + j)]
        dim.para[1] <- j
        start.para[1] <- i + 1
        i <- i + j
    }

    j <- ncol(feature.pi)
    if (j > 0) {
        logitPi <- logitPi + feature.pi %*% para[(i + 1):(i + j)]
        dim.para[2] <- j
        start.para[2] <- i + 1
        i <- i + j
    }
    return(list(logMu = logMu, logitPi = logitPi, dim.para = dim.para, start.para = start.para))
}


zinb.loglik <- function(Y, mu, theta, logitPi) {
    # log-probabilities of counts under the NB model
    logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

    # contribution of zero inflation
    lognorm <- -log1pexp(logitPi)

    # log-likelihood
    sum(logPnb[Y > 0]) + sum(logPnb[Y == 0] + log1pexp(logitPi[Y == 0] - logPnb[Y == 0])) + sum(lognorm)
}

zinb.loglik.dispersion <- function(zeta, Y, mu, logitPi) {
    zinb.loglik(Y, mu, exp(zeta), logitPi)
}

zinb.loglik.dispersion.gradient <- function(zeta, Y, mu, logitPi) {
    theta <- exp(zeta)

    Y0 <- Y <= 0
    Y1 <- Y > 0
    has0 <- !is.na(match(TRUE, Y0))
    has1 <- !is.na(match(TRUE, Y1))

    grad <- 0
    if (has1) {
        grad <- grad + sum(theta * (digamma(Y[Y1] + theta) - digamma(theta) + zeta - log(mu[Y1] + theta) + 1 - (Y[Y1] + theta)/(mu[Y1] + theta)))
    }

    if (has0) {
        logPnb <- suppressWarnings(dnbinom(0, size = theta, mu = mu[Y0], log = TRUE))
        grad <- grad + sum(theta * (zeta - log(mu[Y0] + theta) + 1 - theta/(mu[Y0] + theta))/(1 + exp(logitPi[Y0] - logPnb)))
    }
    return(grad)
}
zinb.loglik.regression <- function(para, Y, feature.mu = matrix(nrow = length(Y), ncol = 0), offset_mu = matrix(0, nrow = length(Y), ncol = 1), feature.pi = matrix(nrow = length(Y),
    ncol = 0), offset_pi = matrix(0, nrow = length(Y), ncol = 1), C.theta = matrix(0, nrow = length(Y), ncol = 1), epsilon = 0) {

    r <- zinb.regression.parseModel(para = para, feature.mu = feature.mu, offset_mu = offset_mu, feature.pi = feature.pi, offset_pi = offset_pi)

    z <- zinb.loglik(Y, exp(r$logMu), exp(C.theta), r$logitPi)

    # Penalty
    z <- z - sum(epsilon * para^2)/2
    return(z)
}

zinb.loglik.regression.gradient <- function(para, Y, feature.mu = matrix(nrow = length(Y), ncol = 0), offset_mu = matrix(0, nrow = length(Y), ncol = 1), feature.pi = matrix(nrow = length(Y),
    ncol = 0), offset_pi = matrix(0, nrow = length(Y), ncol = 1), C.theta = matrix(0, nrow = length(Y), ncol = 1), epsilon = 0) {
    # Parse the model
    r <- zinb.regression.parseModel(para = para, feature.mu = feature.mu, offset_mu = offset_mu, feature.pi = feature.pi, offset_pi = offset_pi)

    theta <- exp(C.theta)
    mu <- exp(r$logMu)
    n <- length(Y)

    Y0 <- Y <= 0
    Y1 <- Y > 0
    has0 <- !is.na(match(TRUE, Y0))
    has1 <- !is.na(match(TRUE, Y1))

    need.wres.mu <- r$dim.para[1] > 0
    need.wres.pi <- r$dim.para[2] > 0

    muz <- 1/(1 + exp(-r$logitPi))
    clogdens0 <- dnbinom(0, size = theta[Y0], mu = mu[Y0], log = TRUE)
    lognorm <- -r$logitPi - log1pexp(-r$logitPi)

    dens0 <- muz[Y0] + exp(lognorm[Y0] + clogdens0)

    # mu
    if (need.wres.mu) {
        wres_mu <- numeric(length = n)
        if (has1) {
            wres_mu[Y1] <- Y[Y1] - mu[Y1] * (Y[Y1] + theta[Y1])/(mu[Y1] + theta[Y1])
        }
        if (has0) {
            wres_mu[Y0] <- -exp(-log(dens0) + lognorm[Y0] + clogdens0 + C.theta[Y0] - log(mu[Y0] + theta[Y0]) + log(mu[Y0]))
        }
    }

    # pi
    if (need.wres.pi) {
        wres_pi <- numeric(length = n)
        if (has1) {
            wres_pi[Y1] <- -muz[Y1]
        }
        if (has0) {
            wres_pi[Y0] <- (1 - exp(clogdens0)) * muz[Y0] * (1 - muz[Y0])/dens0
        }
    }


    grad <- numeric(0)

    # a_mu
    if (r$dim.para[1] > 0) {
        istart <- r$start.para[1]
        iend <- r$start.para[1] + r$dim.para[1] - 1
        grad <- c(grad, colSums(wres_mu * feature.mu) - epsilon[istart:iend] * para[istart:iend])
    }

    # a_pi
    if (r$dim.para[2] > 0) {
        istart <- r$start.para[2]
        iend <- r$start.para[2] + r$dim.para[2] - 1
        grad <- c(grad, colSums(wres_pi * feature.pi) - epsilon[istart:iend] * para[istart:iend])
    }

    return(grad)
}



optimizeDispersion <- function(m, mu, logitPi, epsilon, Y) {
    g <- suppressWarnings(optimize(f = zinb.loglik.dispersion, Y = Y, mu = mu, logitPi = logitPi, maximum = TRUE, interval = c(-50, 50)))
    zeta <- rep(g$maximum, m)

    locfun <- function(logt) {
        s <- sum(unlist(bplapply(seq(m), function(i) {
            zinb.loglik.dispersion(logt[i], Y[i, ], mu[i, ], logitPi[i, ])
        }, BPPARAM = BPPARAM)))
        if (m > 1) {
            s <- s - epsilon * var(logt)/2
        }
        return(s)
    }

    locgrad <- function(logt) {
        s <- unlist(bplapply(seq(m), function(i) {
            zinb.loglik.dispersion.gradient(logt[i], Y[i, ], mu[i, ], logitPi[i, ])
        }, BPPARAM = BPPARAM))
        if (m > 1) {
            s <- s - epsilon * (logt - mean(logt))/(m - 1)
        }
        return(s)
    }
    lower_bounds <- rep(-50, m)
    upper_bounds <- rep(50, m)
    zeta <- optim(par = zeta, fn = locfun, gr = locgrad, control = list(fnscale = -1, trace = 0, maxit = 25), method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)$par
    return(zeta)
}
