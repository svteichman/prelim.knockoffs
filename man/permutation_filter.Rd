% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutation_filter.R
\name{permutation_filter}
\alias{permutation_filter}
\title{Runs permutation procedure}
\usage{
permutation_filter(X, y, fdr = 0.2, plus = TRUE)
}
\arguments{
\item{X}{n-by-p design matrix, with \eqn{n > p}.}

\item{y}{response vector of length n.}

\item{fdr}{desired false discovery rate.}

\item{plus}{True for knockoff+ procedure, false for knockoff procedure.}
}
\value{
A list containing:
 \item{X}{n-by-p design matrix, rescaled so that \eqn{||X_j||^2_2 = 1}, and augmented if \eqn{n < 2p}.}
 \item{X_perm}{n-by-p matrix of knockoff variables.}
 \item{y}{vector of observed responses (augmented if \eqn{n < 2p}). }
 \item{W}{vector of test statistics.}
 \item{thresh}{threshold for variable selection.}
 \item{Selected}{list of selected variables.}
}
\description{
Runs knockoff procedure with permuted variables instead of knockoff variables
for Fixed-X from Barber & Candes (2015).
}
\details{
This runs the knockoff procedure with permuted variables instead of knockoff variables
for a given X matrix and y vector. This controls the false discovery rate at a given level.
}
\examples{
p <- 50; n <- 100; k <- 5
X <- matrix(stats::rnorm(n*p), nrow = n)
true_covar <- sample(p, k)
beta <- 5 * (1:p \%in\% true_covar)
y <- X \%*\% beta + stats::rnorm(n, mean = 0, sd = 1)
results <- permutation_filter(X, y, 0.20, TRUE)

}
\references{
Barber and Candes,
Controlling the false discovery rate via knockoffs.
Ann. Statist. 43 (2015), no. 5, 2055--2085.
https://projecteuclid.org/euclid.aos/1438606853
}
