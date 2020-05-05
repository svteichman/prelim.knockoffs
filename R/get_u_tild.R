#' Construct \eqn{\widetilde{U}}
#'
#' This function constructs \eqn{\widetilde{U}}, a nxp orthonormal matrix that is orthogonal to the column
#' space of X, such that \eqn{X^T\widetilde{U} = 0}.
#'
#' @param X Design matrix to create knockoffs for.
#' @param randomize If true, U tilde is a random matrix. If false, U tilde is the second half of the
#' matrix Q, where Q is part of the QR decomposition of an nx2p matrix \eqn{[U\ 0]}, where \eqn{X = UDV^T}.
#' default is false.
#'
#' @return The orthonormal nxp matrix \eqn{\widetilde{U}}.
#'
#' @export
get_u_tild <- function(X, randomize = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  X_svd <- svd_sign(X)
  u <- X_svd$u
  u_aug <- cbind(u, matrix(0,n,p))
  Q <- qr.Q(qr(u_aug))
  u_tild <- Q[,(p+1):(2*p)]
  if (randomize == TRUE) {
    rand_mat <- matrix(stats::rnorm(p*p),nrow=p)
    rand_Q = qr.Q(qr(rand_mat))
    u_tild <- u_tild %*% rand_Q
  }
  return(u_tild)
}
