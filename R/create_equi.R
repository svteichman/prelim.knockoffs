#' Equi-correlated Knockoffs
#' Creates fixed-X knockoffs using the equi-correlated method
#'
#' @param X Design matrix to create knockoffs for.
#'
#' @return The knockoff nxp design matrix.
#'
#' @export
create_equi <- function(X) {
  rank <- qr(X)$rank
  p <- ncol(X)
  if (rank < p) stop('X is not full rank. Equicorrelated knockoffs will have low power because s will be near 0.')

  X_svd <- svd_sign(X)
  u <- X_svd$u
  d <- X_svd$d
  v <- X_svd$v
  u_tild <- get_u_tild(X)

  eigen <- d^2
  lambda_min <- min(eigen)
  s <- min(2*lambda_min, 1)

  C_part <- pmax(0, 2*s - (s/d)^2)
  X_ko <- (u %*% diag(d - s/d) + u_tild %*% diag(sqrt(C_part))) %*% t(v)
  return(X_ko)
}


