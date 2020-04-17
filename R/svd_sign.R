svd_sign <- function(X) {
  X_svd <- svd(X)

  p <- ncol(X)
  absmax <- function(x) { x[which.max( abs(x) )[1]]}
  cols_max_neg <- apply(X_svd$u, 2, absmax) < 0
  X_svd$u[,cols_max_neg] <- -X_svd$u[,cols_max_neg]
  X_svd$v[,cols_max_neg] <- -X_svd$v[,cols_max_neg]
  return(X_svd)
}
