#' Run Benjamini-Hochberg procedure
#'
#' Runs the either the original Benjamini-Hochberg procedure, the BH procedure with a log correction,
#' or the BH procedure with whitened noise.
#'
#' @param X A nxp design matrix.
#' @param y A vector of n responses.
#' @param fdr The desired false discovery rate.
#' @param method either "base", "log", or "white_noise" (default: "base").
#'
#' @return A list containing:
#'  \item{X}{n-by-p design matrix.}
#'  \item{y}{vector of observed responses. }
#'  \item{Z}{vector of test statistics.}
#'  \item{thres}{data-dependent threshold.}
#'  \item{selected}{list of selected variables.}
#'
#' @export
run_BHq <- function(X, y, fdr, method = c('base','log','white_noise')) {
  res = switch(match.arg(method),
                   "base" = run_BHq_base(X,y,fdr),
                   "log"  = run_BHq_log(X,y,fdr),
                   "white_noise"  = run_BHq_white_noise(X,y,fdr))
  return(res)
}
