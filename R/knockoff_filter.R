#' Runs knockoff procedure
#'
#' Runs knockoff procedure for Fixed-X from Barber & Candes (2015).
#'
#' @param X n-by-p design matrix, with \eqn{n > p}.
#' @param y response vector of length n.
#' @param method either "equi" or "sdp" (default: "sdp").
#' These are two different methods to generate the p-dimensional s vector, which is used to constuct the knockoff variables.
#' @param fdr desired false discovery rate.
#' @param plus True for knockoff+ procedure, false for knockoff procedure.
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
#' rate at a given level.
#'
#' @examples
#' p <- 50; n <- 100; k <- 5
#' X <- matrix(stats::rnorm(n*p), nrow = n)
#' true_covar <- sample(p, k)
#' beta <- 5 * (1:p %in% true_covar)
#' y <- X %*% beta + stats::rnorm(n, mean = 0, sd = 1)
#' results <- knockoff_filter(X, y, 'equi', 0.20, TRUE)
#'
#' @export
knockoff_filter <- function(X, y, method = c('sdp','equi'), fdr = .20, plus = T) {
  n <- nrow(X)
  p <- ncol(X)
  n_y <- length(y)
  if (n != n_y) stop("X and  y must have the same number of observations")

  knock <- create_knockoffs(X, y, method = method)
  X <- knock$X
  Xk <- knock$Xk
  y <- knock$y

  W <- get_stat_lambda_max(X, Xk, y)
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

