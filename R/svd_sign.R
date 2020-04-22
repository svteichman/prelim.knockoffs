#' SVD with fixed sign
#' Performs SVD and then adjusts the sign of each column of the U matrix so that the value
#' with the greatest magnitude is positive. It fixes the signs of the V matrix accordingly.
#'
#' @param X Design matrix to create knockoffs for.
#'
#' @return The SVD of X: U, a vector of singular values d, V, such that \eqn{X = Udiag\{d\}V^T}.
#'
#' @export
svd_sign <- function(X) {
  X_svd <- svd(X)

  p <- ncol(X)
  absmax <- function(x) { x[which.max( abs(x) )[1]]}
  cols_max_neg <- apply(X_svd$u, 2, absmax) < 0
  X_svd$u[,cols_max_neg] <- -X_svd$u[,cols_max_neg]
  X_svd$v[,cols_max_neg] <- -X_svd$v[,cols_max_neg]
  return(X_svd)
}
