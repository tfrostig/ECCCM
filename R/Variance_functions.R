
findOmegaBeta <- function(beta.mc, omega, sigma, n.o) {
  omega.beta       <-  cbind('omega_beta'     = omega %*% beta.mc,
                             'omega_beta_var' = (sigma^2 / n.o) * diag(omega %*% omega %*% omega))
}


#' estimateVarAdd - estimation of additional variance
#' @param beta.mc The estimate of coefficient of joint regression (transformed marginal beta)
#' @param beta.omega Transformed `beta.mc` basically, just `beta.omega %*% omega`, where omega is the estimated inverse covariance matrix.
#' @param x.r The reference panel.
#' @param ind.beta.mc Indices of which `beta.mc` to consider, in the bayesian scenario determined by vector gamma.
#' @param ind.beta.omega Indices of which `beta.omega` to consider
#' @param omega - Inverse of estimated covaraince matrix
#' @param n.o - Number of observation in original study
#' @return A list containing the varius variance terms, including the naive. `total` is the combined added variance from
#' using the reference panel
#' @export

estimateVarAdd <- function(beta.mc, beta.omega, x.r, ind.beta.mc, ind.beta.omega, omega, n.o) {
  ### Term 1
  term.1.var.r <- varFirstTerm(beta.mc, omega, x = x.r, ind = ind.beta.mc)  / nrow(x.r)
  ### Term 2
  term.2.var.r <- varSecondTerm(beta.omega, x = x.r, ind = ind.beta.omega) / n.o
  ### Naive var
  return(list('var.term.1' = term.1.var.r,
              'var.term.2' = term.2.var.r,
              'total'      = (term.1.var.r + term.2.var.r)))

}
