#' Knockoff threshold
#'
#' Calculates the threshold T for variable selection for the Knockoff and Knockoff+ procedures.
#'
#' @param W the set of p test statistics.
#' @param q the desired FDR.
#' @param plus a boolean T/F, true for the Knockoff+ procedure and false for the Knockoff procedure.
#'
#' @return A threshold to use for variable selection.
#'
#' @export
compute_threshold <- function(W, q, plus) {
  # get set of possible threshold values
  set_W <- sort(unique(abs(W)))
  set_W <- set_W[set_W > 0]

  # account for the additional 1 in the numerator for Knockoff+ procedure
  if (plus == T) {offset = 1}
  if (plus == F) {offset = 0}

  ratio <- function(W, t, offset) {
    num <- offset + sum(W <= -t)
    denom <- max(sum(W >= t),1)
    return(num/denom)
  }

  t_ratios <- lapply(set_W, function(t) ratio(W, t, offset))
  ind <- ifelse(sum(t_ratios < q) == 0, 0, min(which(t_ratios < q)))
  return(set_W[ind])
}

