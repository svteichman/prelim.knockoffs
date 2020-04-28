#' SDP using Rdsdp
#' Solves the optimation problem to choose s for the SDP knockoff method, using a SDP solver
#' implemented in Rdsdp. Specifically, this solves the problem:
#' \eqn{\text{Maximize}\sum_{j=1}^p s_j} subject to the linear constraint \eqn{0 \leq s_j \leq 1\ \forall\ s_j}
#' and the PSD constraint \eqn{2\Sigma - \text{diag}\{s\} \succeq 0}.
#'
#' @param Sigma Normalized empirical covariance matrix to optimize s for.
#' @param gaptol The tolerance for the duality gap as a fraction of the objective functions (a parameter in Rdsdp).
#' @param maxit The maximum number of iterations allowed (a paramter in Rdsdp).
#'
#' @return The vector of s values that solve the above optimization problem.
#'
#' @export
solve_sdp <- function(Sigma, gaptol = 1e-6, maxit = 1000) {
  p <- nrow(Sigma)
  corr <- stats::cov2cor(Sigma)
  # Make linear cone blocks
  # Constraint that each s_j >= 0
  A1 <- diag(-1,p)
  C1 <- rep(0,p)
  # Constraint that each s_j <= 1
  A2 <- diag(1,p)
  C2 <- rep(1,p)

  # Make PSD cone block
  A_PSD <- matrix(0, nrow = p, ncol = p*p)
  for (i in 1:p) {
    A_PSD[i,((i-1)*p+i)] <- 1
  }
  C_PSD <- c(2*corr)

  # Combine linear and PSD cone blocks
  A <- cbind(A1, A2, A_PSD)
  C <- matrix(c(C1, C2, C_PSD),1)

  # Describe the size of each block (linear, then PSD)
  K <- list(p, 2*p)
  names(K) <- c('s','l')

  # Make right hand side of constraint (all 1's)
  b <- rep(1, p)

  # Set options and run Rdsdp
  options <- list(gaptol, maxit, 0, 0, 0)
  names(options) <- c("gaptol","maxit","logsummary","outputstats","print")
  sol = Rdsdp::dsdp(A,b,C,K,options)
  s <- sol$y

  # Fix numerical issues in domain
  s[s > 1] <- 1
  s[s < 0] <- 0

  # Check if solution is feasible
  if(sol$STATS$stype == "Infeasible") {
    warning('Dual problem is infeasible.')
  }
  if(sol$STATS$stype == "Unbounded") {
    warning('Dual problem is unbounded.')
  }

  # Check if resulting matrix is PSD. If not, lower the s values gradually
  s_factor <- 1e-8
  if (!check_PSD(2*corr - diag(s))) {
    PSD <- FALSE
    lim <- 0.1
    while (PSD == FALSE & s_factor <= lim) {
      if (check_PSD(2*corr - diag(s*(1-s_factor)), tol = 1e-9)) {
        PSD <- TRUE
        s <- s*(1-s_factor)
      }
      s_factor <- s_factor*10
    }
  }
  if (s_factor == 1) {
    s <- s*0
  }
  s[s<0] <- 0

  # Check for power
  if (sum(s) == 0) {
    warning("Knockoffs have no power.")
  }

  # return optimal s
  return(s*diag(Sigma))
}
