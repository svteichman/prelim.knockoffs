#' Create permuted variables
#'
#' Creates a permuted version of the fixed design matrix.
#'
#' @param X n-by-p design matrix, with \eqn{n > p}.
#' @param y response vector of length n.
#' This is used in the setting in which \eqn{n < 2p}, to estimate sigma hat and generate additional rows in X and y.
#' @return A list containing:
#'  \item{X}{n-by-p design matrix, rescaled so that \eqn{||X_j||^2_2 = 1}, and augmented if \eqn{n < 2p}.}
#'  \item{X_perm}{n-by-p matrix of permuted variables.}
#'  \item{y}{vector of observed responses (augmented if \eqn{n < 2p}). }
#'
#' @source
#' Barber and Candes,
#' Controlling the false discovery rate via knockoffs.
#' Ann. Statist. 43 (2015), no. 5, 2055--2085.
#' https://projecteuclid.org/euclid.aos/1438606853
#'
#' @details
#' This creates p permuted variables for a fixed design matrix of size n by p.
#'
#'
#' @export
create_permutations <- function(X, y) {
  #method = match.arg(method)
  n = nrow(X)
  p = ncol(X)

  if (n <= p) stop("number of observations (n) must be greater than number of variables (p)")
  if (n < 2*p) {
    warning('n < 2p, X and y will be augmented with (2p-n) rows',immediate.=T)

    extra <- 2*p - n #number of rows to augment by
    beta_ols <- solve(t(X) %*% X) %*% t(X) %*% y
    sig_hat <- sqrt(sum((y - X %*% beta_ols)^2)/(n-p)) #estimate of sigma

    X <- rbind(X, matrix(0,nrow = extra, ncol = p)) #augment X with 0's
    set.seed(0)
    y_extra <- stats::rnorm(extra, mean = 0, sd = sig_hat)
    y <- c(y, y_extra)
  }

  X.scaled <- scale(X, center=F, scale=sqrt(colSums(X^2)))
  X_scaled <- X.scaled[,]

  perm <- sample(1:n ,n, replace = F)
  X_perm <- X_scaled[perm,]
  res <- list(X_scaled, X_perm, y)
  names(res) <- c("X","X_perm","y")
  return(res)
}

