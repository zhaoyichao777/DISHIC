create_model <- function(A, B, alpha_mu, alpha_pi, beta_mu, beta_pi, zeta, m, n,epsilon,
                      epsilon_alpha_mu, epsilon_beta_mu, epsilon_alpha_pi,
                      epsilon_beta_pi,epsilon_zeta) {
    if (missing(m)) {
        if (!missing(A)) {
            m <- NROW(A)
        } else if (!missing(beta_mu)) {
            m <- NCOL(beta_mu)
        } else if (!missing(beta_pi)) {
            m <- NCOL(beta_pi)
        }
    }
    if (missing(n)) {
        if (!missing(B)) {
            n <- NROW(B)
        } else if (!missing(alpha_mu)) {
            n <- NCOL(alpha_mu)
        } else if (!missing(alpha_pi)) {
            n <- NCOL(alpha_pi)
        }
    }

    if (missing(A)) {
        A <- matrix(1, nrow = m, ncol = 1) # intercept
    }
    if (missing(B)) {
        B <- matrix(1, nrow = n, ncol = 1)
    }
    if (missing(alpha_mu)) {
        alpha_mu <- matrix(0, nrow = NCOL(A), ncol = n)
    }
    if (missing(alpha_pi)) {
        alpha_pi <- matrix(0, nrow = NCOL(A), ncol = n)
    }
    if (missing(beta_mu)) {
        beta_mu <- matrix(0, nrow = NCOL(B), ncol = m)
    }
    if (missing(beta_pi)) {
        beta_pi <- matrix(0, nrow = NCOL(B), ncol = m)
    }
    if (missing(zeta)) {
        zeta <- matrix(0, nrow = m, ncol = 1)
    }

    if (missing(epsilon)) {
        epsilon <- m
    }
    if (missing(epsilon_alpha_mu)) {
        epsilon_alpha_mu <- epsilon / n
        e <- rep(epsilon_alpha_mu, NCOL(A))
        e[1] <- 0
        epsilon_alpha_mu <- e
    }
    if (missing(epsilon_beta_mu)) {
        epsilon_beta_mu <- epsilon / m
        e <- rep(epsilon_beta_mu, NCOL(B))
        e[1] <- 0
        epsilon_beta_mu <- e
    }
    if (missing(epsilon_alpha_pi)) {
        epsilon_alpha_pi <- epsilon / n
        e <- rep(epsilon_alpha_pi, NCOL(A))
        e[1] <- 1e-3
        epsilon_alpha_pi <- e
    }
    if (missing(epsilon_beta_pi)) {
        epsilon_beta_pi <- epsilon / m
        e <- rep(epsilon_beta_pi, NCOL(B))
        e[1] <- 1e-3
        epsilon_beta_pi <- e
    }
    if (missing(epsilon_zeta)) {
        epsilon_zeta <- epsilon
    }

    # if (missing(O_mu)) {
    #     O_mu <- matrix(0, nrow = n, ncol = J)
    # }


    obj <- list(
        A = A,
        B = B,
        alpha_mu = alpha_mu,
        alpha_pi = alpha_pi,
        beta_mu = beta_mu,
        beta_pi = beta_pi,
        zeta = zeta,
        m = m,
        n = n,
        epsilon=epsilon,
        epsilon_alpha_mu=epsilon_alpha_mu,
        epsilon_beta_mu=epsilon_beta_mu,
        epsilon_alpha_pi=epsilon_alpha_pi,
        epsilon_beta_pi=epsilon_beta_pi,
        epsilon_zeta=epsilon_zeta
    )
    class(obj) <- "model"
    return(obj)
}

# print model
print.model <- function(model) {
    cat("Model Details:\n")
    cat("A:", model$A, "\n")
    cat("B:", model$B, "\n")
    cat("alpha_mu:", model$alpha_mu, "\n")
    cat("alpha_pi:", model$alpha_pi, "\n")
    cat("beta_mu:", model$beta_mu, "\n")
    cat("beta_pi:", model$beta_pi, "\n")
    cat("theta:", model$theta, "\n")
    cat("m:", model$m, "\n")
}
