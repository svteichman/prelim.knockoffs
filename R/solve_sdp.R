solve_sdp <- function(Sigma, gaptol = 1e-6, maxit = 1000, verbose = FALSE) {
  p <- nrow(Sigma)

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
  C_PSD <- c(2*Sigma)

  # Combine linear and PSD cone blocks
  A <- cbind(A1, A2, A_PSD)
  C <- c(C1, C2, C_PSD)

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

  # COMPENSATE FOR NUMERICAL ISSUES - FEASIBILITY
  # COME BACK HERE IF PROCEDURE IS AT ALL DIFFERENT FROM THEIR METHOD

  # Check for power
  if (sum(s) == 0) {
    warning("Knockoffs have no power.")
  }

  # return optimal s
  return(s)
}
