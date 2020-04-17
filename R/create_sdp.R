#' SDP Knockoffs
#' Creates fixed-X knockoffs using the SDP method
#'
#' @param X Design matrix to create knockoffs for.
#'
#' @return The knockoff nxp design matrix.
create_sdp <- function(X) {
  rank <- qr(X)$rank
  p <- ncol(X)
  if (rank < p) warning('X is not full rank. Proceeding with knockoffs',immediate.=T)

  X_svd <- svd_sign(X)
  u <- X_svd$u
  d <- X_svd$d
  v <- X_svd$v
  u_tild <- get_u_tild(X)

  Sigma <- v %*% diag(d^2) %*% t(v)
  Sigma_inv <- v %*% diag(1/(d^2)) %*% t(v)

  s <- create.solve_sdp(Sigma)
  s_small <- s < 1e-5
  s[s_small] <- 0

  CtC <- 2*diag(s)- diag(s) %*% Sigma_inv %*% diag(s)
  CtC_svd <- svd_sign(CtC)
  CtC_d <- CtC_svd$d
  CtC_v <- CtC_svd$v
  X_ko <- X %*% (diag(p) - Sigma_inv %*% diag(s)) + u_tild %*% diag(sqrt(CtC_d)) %*% t(CtC_v)

  return(X_ko)
}
