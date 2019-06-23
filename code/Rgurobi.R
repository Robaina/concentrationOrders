solveGurobiProblem <- function(c, Gp, hp, Ap, bp,
  binaryVariables=NULL, modelsense="min", N_mets=NULL) {

  # Solves a MILP through Gurobi using a GLPK input

  library(gurobi)

  N_ineqs <- length(hp)
  N_eqs <- length(bp)
  N_vars <- ncol(Gp)
  vartype <- rep("C", N_vars)

  if (!is.null(binaryVariables)) {
    vartype[binaryVariables] <- "B"
  }

  # Prepare gurobi input
  model <- list()
  params <- list()
  model$A <- rbind(Gp, Ap)
  model$sense <- c(rep("<", N_ineqs), rep("=", N_eqs))
  model$rhs <- c(hp, bp)
  model$obj <- c
  model$modelsense <- modelsense
  model$vtype <- vartype
  #model$lb <- v_lb
  #model$ub <- v_ub
  params$OutputFlag <- 1

  gursol <- gurobi(model, params)
  res <- list()
  res$logx <- res$x[1:N_mets]
  res$status <- gursol$status
  return(res)
}
