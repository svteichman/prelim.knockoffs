#' Create knockoffs
#'
#' Creates knockoff variables for fixed design matrix.
#'
#' @param X n-by-p design matrix, with \eqn{n > p}.
#' @param method either "equi" or "sdp" (default: "sdp").
#' These are two different methods to generate the p-dimensional s vector, which is used to constuct the knockoff variables.
#' @param y response vector of length n.
#' This is used in the setting in which \eqn{n < 2p}, to estimate sigma hat and generate additional rows in X and y.
#' @return A list containing:
#'  \item{X}{n-by-p design matrix, rescaled so that \eqn{||X_j||^2_2 = 1}, and augmented if \eqn{n < 2p}.}
#'  \item{X_knock}{n-by-p matrix of knockoff variables.}
#'  \item{y}{vector of observed responses (augmented if \eqn{n < 2p}). }
#'
#' @references
#'   Barber and Candes,
#'   Controlling the false discovery rate via knockoffs.
#'   Ann. Statist. 43 (2015), no. 5, 2055--2085.
#'   \href{https://projecteuclid.org/euclid.aos/1438606853}
#'
#' @details
#' This creates p knockoff variables for a fixed design matrix of size n by p. These knockoffs can be
#' used to provably control the false discovery rate (FDR) for the linear model with iid Gaussian errors.
#'
#' @examples
#' p <- 50; n <- 100; k <- 5
#' X <- matrix(rnorm(n*p), nrow = n)
#' true_covar <- sample(p, k)
#' beta <- 5 * (1:p %in% true_covar)
#' y <- X %*% beta + rnorm(n, mean = 0, sd = 1)
#' knock <- create_knockoffs(X, y, method = 'sdp')
#'
#' @export
create_knockoffs <- function(X, y, method = c('sdp','equi')) {
  method = match.arg(method)
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
    y_extra <- rnorm(extra, mean = 0, sd = sig_hat)
    y <- c(y, y_extra)
  }

  scale_factor <- sqrt(colSums(X^2))
  X_scaled <- t(t(X)/scale_factor) #normalize X such that X^TX = Sigma, Sigma_jj = 1 for all j

  X_knock = switch(match.arg(method),
              "equi" = create_equi(X_scaled),
              "sdp"  = create_sdp(X_scaled)
  )

  res <- list(X_scaled, X_knock, y)
  names(res) <- c("X","Xk","y")
  return(res)
}





