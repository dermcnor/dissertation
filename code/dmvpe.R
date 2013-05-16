dmvpe <- function (x, mean, sigma, beta, log = FALSE){
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
	if (missing(beta)) {
        beta <- 1
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
	if (beta <= 0) {
        stop("beta must be positive")
    }
    distval <- mahalanobis(x, center = mean, cov = sigma)^beta
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
	k <- ncol(x)
    logretval <- log(k) + lgamma(k/2) - lgamma(1 + k/(2 * beta)) - (1 + k/(2 * beta)) * log(2) -(k * log(pi) + logdet + distval)/2
    if (log){
        return(logretval)
	}
    exp(logretval)
}
