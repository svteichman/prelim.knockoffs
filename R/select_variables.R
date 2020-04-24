#' Variable selection
#'
#' Selects variables using the W_j statistics and threshold provided.
#'
#' @param W the set of p test statistics.
#' @param thresh threshold for variable selection.
#'
#' @return A list of selected variables.
#'
#' @export
select_variables <- function(W, thresh) {
  selected <- sort(which(W >= thresh))
  return(selected)
}
