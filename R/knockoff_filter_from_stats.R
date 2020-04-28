#' Runs knockoff procedure abbreviated
#'
#' Runs knockoff procedure for Fixed-X from Barber & Candes (2015) but takes in X, knockoffs, y,
#' and statistics, and computes threshold and returns selected variables. This can be used when
#' using both knockoff and knockoff+ procedures, the output from one procedure can be put into this
#' function, which will compute the threshold and return selected variables from the other procedure.
#'
#' @param X n-by-p design matrix, with \eqn{n > p}.
#' @param y response vector of length n.
#' @param fdr desired false discovery rate.
#' @param plus True for knockoff+ procedure, false for knockoff procedure.
#' @param Xk n-by-p knockoff matrix.
#' @param stat W test statistics.
#'
#' @return A list containing:
#'  \item{X}{n-by-p design matrix, rescaled so that \eqn{||X_j||^2_2 = 1}, and augmented if \eqn{n < 2p}.}
#'  \item{Xk}{n-by-p matrix of knockoff variables.}
#'  \item{y}{vector of observed responses (augmented if \eqn{n < 2p}). }
#'  \item{W}{vector of test statistics.}
#'  \item{thresh}{threshold for variable selection.}
#'  \item{Selected}{list of selected variables.}
#'
#' @references
#' Barber and Candes,
#' Controlling the false discovery rate via knockoffs.
#' Ann. Statist. 43 (2015), no. 5, 2055--2085.
#' https://projecteuclid.org/euclid.aos/1438606853
#'
#' @details
#' This runs the knockoff filter for a given X matrix and y vector. This controls the false discovery
#' rate at a given level. This takes in the knockoffs and test statistics and returns the threshold
#' and selected variables. It should be used to reduce computation when running both knockoff and
#' knockoff+ procedures on the same data.
#'
#' @examples
#' p <- 50; n <- 100; k <- 5
#' X <- matrix(stats::rnorm(n*p), nrow = n)
#' true_covar <- sample(p, k)
#' beta <- 5 * (1:p %in% true_covar)
#' y <- X %*% beta + stats::rnorm(n, mean = 0, sd = 1)
#' results <- knockoff_filter(X, y, 0.20, TRUE, Xk, W)
#'
#' @export
knockoff_filter_from_stats <- function(X, y, fdr = .20, plus = TRUE, Xk, W) {
  thresh <- compute_threshold(W, fdr, plus)

  selected = select_variables(W, thresh)
  if (!is.null(names(X))) {
    names(selected) = names(X)[selected] }

  # Package up the results.
  res <- list(call = match.call(),
              X = X,
              Xk = Xk,
              y = y,
              statistic = W,
              threshold = thresh,
              selected = selected)
  return(res)
}

