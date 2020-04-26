#' Run Benjamini-Hochberg procedure with added noise
#'
#' Runs the altered Benjamini-Hochberg procedure, which adds a mean zero normal RV to each
#' beta hat to make the z-scores independent of each other.
#'
#' @param X A nxp design matrix.
#' @param y A vector of n responses.
#' @param fdr The desired false discovery rate.
#'
#' @return A list containing:
#'  \item{X}{n-by-p design matrix.}
#'  \item{y}{vector of observed responses. }
#'  \item{Z}{vector of test statistics.}
#'  \item{thres}{data-dependent threshold.}
#'  \item{selected}{list of selected variables.}
#'
#' @export
run_BHq_white_noise <- function(X, y, fdr) {
  n <- nrow(X)
  p <- ncol(X)
  Sigma_inv <- solve(t(X) %*% X )
  Sig_inv_diag <- diag(Sigma_inv)
  beta_ols <- Sigma_inv %*% t(X) %*% y
  sig_hat <- sqrt(sum((y - X %*% beta_ols)^2)/(n-p)) #estimate of sigma

  lambda_min <- min(eigen(Sigma)$values)
  cov <- sig_hat*((1/lambda_min)*diag(nrow=p) - Sigma_inv)
  Z_prime <- mvtnorm::rmvnorm(1, rep(0,p), cov)

  Z <- (beta_ols[,1]+Z_prime[1,])/(sig_hat*sqrt(1/lambda_min))
  t_max <- max(abs(Z))
  t_diff <- t_max/1000
  t_seq <- seq(t_diff,t_max,t_diff)

  get_ratio <- function(Z, t, p) {
    num <- p*(1-stats::pchisq(t^2, 1))
    denom <- sum(abs(Z) >= t)
    return(num/denom)
  }
  t_ratio <- sapply(t_seq, function(t) get_ratio(Z, t, p))
  thresh <- t_seq[min(which(t_ratio <= fdr))]

  selected <- which(abs(Z) >= thresh)
  res <- list(X = X,
              y = y,
              statistic = Z,
              threshold = thresh,
              selected = selected)

  return(res)
}
