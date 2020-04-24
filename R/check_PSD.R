#' Check Positive semi-definite
#' Checks if the given matrix is positive semi-definite.
#'
#' @param M Matrix to check whether it is positive semi-definite.
#'
#' @return TRUE if the matrix is PSD, FALSE otherwise.
#'
#' @export
check_PSD <- function(M, tol = 1e-9) {
  lambda_min <- min(eigen(M)$values)
  PSD <- lambda_min > (tol*10)
  return(PSD)
}
