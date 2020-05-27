#' Check Positive definite
#' Checks if the given matrix is positive definite.
#'
#' @param M Matrix to check whether it is positive definite.
#' @param tol Tolerance when comparing smallest eigenvalue to 0.
#'
#' @return TRUE if the matrix is PD, FALSE otherwise.
#'
#' @export
check_PD <- function(M, tol = 1e-9) {
  lambda_min <- min(eigen(M)$values)
  PD <- lambda_min > (tol*10)
  return(PD)
}
