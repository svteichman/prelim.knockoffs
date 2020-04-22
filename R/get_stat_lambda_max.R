#' Lambda max test statistic
#' Creates W_j test statistics. Performs penalized regression for a range of lambda values on the augmented
#' X matrix \eqn{[X\ \widetilde{X}]}. Gets \eqn{Z_j = \sup\{\lambda:\hat{\beta}_j \neq 0\}} for each
#' original and knockoff variables. Computes the test statistic \eqn{W_j = max(Z_j, Z_{j+p)*-1\cdot I(Z_j < Z_{j+p})})}.
#'
#' @param X n x p X matrix.
#' @param Xk n x p knockoff matrix.
#' @param y n dimensional y vector.
#'
#' @return A vector of p test statistics \eqn{W_j}.
#'
#' @export
get_stat_lambda_max <- function(X, Xk, y) {
  # construct X matrix
  X_aug <- cbind(X, Xk)
  n <- nrow(X_aug)
  p <- ncol(X_aug)
  X_aug <- scale(X_aug)[,] #standardize

  # construct list of potential lambdas
  num_lambda <- 500
  lambda_max <- max(abs(t(X_aug) %*% y))/n
  lambda_min <- lambda_max/(2e3)
  k <- (0:(num_lambda - 1))/num_lambda
  lambda <- lambda_max * (lambda_min/lambda_max)^k

  # fit penalized regression via glmnet
  mod <- glmnet::glmnet(X_aug, y, family = "gaussian", lambda = lambda)
  first_entry <- function(vect) {
    index <- ifelse(sum(abs(vect) > 0) == 0, 0, min(which(abs(vect) > 0)))
    return(index)
  }
  lam_indices <- apply(mod$beta, 1, function(x) first_entry(x))
  lambda <- c(0, lambda)
  lam_vals <- lambda[lam_indices+1]*n

  # calculate W_j's
  Z <- lam_vals[1:(p/2)]
  Z_tild <- lam_vals[(p/2 + 1):p]
  sign <- Z > Z_tild
  W <- pmax(Z, Z_tild)
  W <- -W*(1-sign) + W*1*sign
  set_zero <- which(W == 0)
  W[set_zero] <- 0

  return(W)
}
