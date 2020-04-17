#' Construct \eqn{\widetilde{U}}
#'
#' This function constructs \eqn{\widetilde{U}}, a nxp orthonormal matrix that is orthogonal to the column
#' space of X, such that \eqn{X^T\widetilde{U} = 0}.
#'
#' @param X Design matrix to create knockoffs for.
#'
#' @return The orthonormal nxp matrix \eqn{\widetilde{U}}.
get_u_tild <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  X_svd <- svd_sign(X)
  u <- X_svd$u
  u_aug <- cbind(u, matrix(0,n,p))
  Q <- qr.Q(qr(u_aug))
  return(Q[,(p+1):(2*p)])
}

