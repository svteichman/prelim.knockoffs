#' SDP Knockoffs
#' Creates fixed-X knockoffs using the SDP method
#'
#' @param X Design matrix to create knockoffs for.
#' @param randomize If true, U tilde is a random matrix. If false, U tilde is the second half of the
#' matrix Q, where Q is part of the QR decomposition of an nx2p matrix \eqn{[U\ 0]}, where \eqn{X = UDV^T}.
#' default is false.
#'
#' @return The knockoff nxp design matrix.
#'
#' @export
create_sdp <- function(X, randomize = FALSE) {
  rank <- qr(X)$rank
  p <- ncol(X)
  if (rank < p) warning('X is not full rank. Proceeding with knockoffs',immediate.=T)

  X_svd <- svd_sign(X)
  u <- X_svd$u
  d <- X_svd$d
  d_inv = 1 / d
  v <- X_svd$v
  u_tild <- get_u_tild(X, randomize = randomize)

  Sigma <- v %*% diag(d^2) %*% t(v)
  Sigma_inv <- v %*% diag(d_inv^2) %*% t(v)

  s <- solve_sdp(Sigma)
  s_small <- s < 1e-5
  s[s_small] <- 0

  C_svd <- svd_sign(2*diag(s) - (diag(s) %*% Sigma_inv %*% diag(s)))
  X_ko <- X - (X %*% Sigma_inv %*% diag(s)) +
    (u_tild %*% diag(sqrt(pmax(0, C_svd$d)))) %*% t(C_svd$v)

  return(X_ko)
}
