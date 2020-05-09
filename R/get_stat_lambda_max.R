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

  # randomly permute original variables and knockoff variables
  # because glmnet fits LASSO with coordinate descent, order of variables affects coefficients
  swap <- rbinom(p/2,1,0.5)
  to_swap <- which(swap == 1)
  X_aug_s <- X_aug
  X_aug_s[,to_swap] <- X_aug[,to_swap + p/2]
  X_aug_s[,to_swap + p/2] <- X_aug[,to_swap]

  # construct list of potential lambdas
  num_lambda <- 500
  lambda_max <- max(abs(t(X_aug_s) %*% y))/n
  lambda_min <- lambda_max/(2e3)
  k <- (0:(num_lambda - 1))/num_lambda
  lambda <- lambda_max * (lambda_min/lambda_max)^k

  # fit penalized regression via glmnet
  mod <- glmnet::glmnet(X_aug_s, y, family = "gaussian", lambda = lambda, intercept=T,
                        standardize=F, standardize.response=F)
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
  sign <- sign(Z - Z_tild)
  W <- pmax(Z, Z_tild)
  W <- W*sign

  W_unswap <- W*(1-2*swap)
  return(W_unswap)
}



swap = rbinom(ncol(X),1,0.5)
swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
X.swap  = X * (1-swap.M) + Xk * swap.M
Xk.swap = X * swap.M + Xk * (1-swap.M)

W = W * (1-2*swap)
