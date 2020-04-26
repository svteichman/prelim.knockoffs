#' Run Benjamini-Hochberg procedure
#'
#' Runs the original Benjamini-Hochberg procedure.
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
run_BHq_base <- function(X, y, fdr) {
  n <- nrow(X)
  p <- ncol(X)
  Sigma_inv <- solve(t(X) %*% X )
  Sig_inv_diag <- diag(Sigma_inv)
  beta_ols <- Sigma_inv %*% t(X) %*% y
  sig_hat <- sqrt(sum((y - X %*% beta_ols)^2)/(n-p)) #estimate of sigma

  Z <- beta_ols/(sig_hat*sqrt(Sig_inv_diag))
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
