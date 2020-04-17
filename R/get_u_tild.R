get_u_tild <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  X_svd <- svd_sign(X)
  u <- X_svd$u
  u_aug <- cbind(u, matrix(0,n,p))
  Q <- qr.Q(qr(u_aug))
  return(Q[,(p+1):(2*p)])
}

