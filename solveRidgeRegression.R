# --------------------------------------------------------------------------------
# This code is based on modifications made to the original zinbwave code.
# Original zinbwave code can be found at: https://github.com/biocore/zinbwave
# Modifications made to adapt to specific data structure and analysis requirements.
# --------------------------------------------------------------------------------
solveRidgeRegression <- function(x, y, beta = rep(0, NCOL(x)), epsilon = 1e-6, family = c("gaussian", "binomial"), offset) {
    # x:col-feature_number,row-cell/binpair
    family <- match.arg(family)

    # loglik
    f <- if (family == "gaussian") {
        function(b) { 
            eta <- x %*% b 
            l <- sum((y - eta)^2) / 2 
            l + sum(epsilon * b^2) / 2 
        }
    } else if (family == "binomial") {
        function(b) { 
            eta <- x %*% b + offset 
            l <- sum(-y * eta + log1pexp(eta))
            l + sum(epsilon * b^2) / 2 
        }
    }

    # gradient of loglik
    g <- if (family == "gaussian") {
        function(b) {
            eta <- x %*% b 
            l <- t(x) %*% (-y + eta)
            l + epsilon * b
        }
    } else if (family == "binomial") {
        function(b) {
            eta <- x %*% b + offset
            l <- t(x) %*% (-y + 1 / (1 + exp(-eta)))
            l + epsilon * b
        }
    }

    # optimize
    m <- optim(fn = f, gr = g, par = beta, control = list(trace = 0), method = "BFGS")
    return(m$par)
}
