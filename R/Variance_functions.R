
findOmegaBeta <- function(beta.mc, omega, sigma, n.o) {
  omega.beta       <-  cbind('omega_beta'     = omega %*% beta.mc,
                             'omega_beta_var' = (sigma^2 / n.o) * diag(omega %*% omega %*% omega))
}

### Old version of matrix estimation
estimateVariance <- function(g.r, sig.r, inv.r, beta.joint.est, n.o, method = 'moment') {
  n.r   <- nrow(g.r)
  p     <- ncol(g.r)
  if (method == 'moment') {
    if (ncol(g.r) == 1) {
      sig.est.vec     <- apply(g.r, 1, function(x) as.vector(x %*% t(x))) ### Strange behaviour of apply with 1 dimension
      A.r             <- drop((t(sig.est.vec) %*% sig.est.vec / n.r) - sig.r * sig.r)
      kron.beta.inv.r <- beta.joint.est * drop(inv.r)
      term.2.var.r    <- A.r * t(kron.beta.inv.r) %*% kron.beta.inv.r / n.r
      ## Term 3
      beta.inv.r      <- drop(t(beta.joint.est) * drop(inv.r))
      term.3.var.r    <- A.r * t(beta.inv.r) %*% beta.inv.r / n.o
    }
    if (ncol(g.r) > 1) {
      ### Variance of covariance matrix
      outer.prod    <- outerProdRow(g.r)
      A.r           <- findCovVar(outer.prod, cov(g.r))
      ## Term 2
      kron.beta.inv.r <- kronecker(t(beta.joint.est), inv.r)
      term.2.var.r    <- kron.beta.inv.r %*% A.r %*% t(kron.beta.inv.r) / n.r
      ## Term 3
      beta.inv.r           <- t(beta.joint.est) %*% inv.r
      kron.beta.inv.diag.r <- kronecker(beta.inv.r, diag(p))
      term.3.var.r         <- kron.beta.inv.diag.r %*% A.r %*% t(kron.beta.inv.diag.r) / n.o
    }
  }
  if (method == 'normal') {
    c1.r          <- ((n.r - p) * (n.r - p - 1) * (n.r - p - 3))^(-1)
    d1.r          <- ((n.r - p) * (n.r - p - 1)^2 * (n.r - p - 3))^(-1)
    ## Term 2
    term.2.var.r  <- n.r^2 *  (c1.ref * drop(t(beta.joint.est) %*% sig.r %*% beta.joint.est) * inv.r +
                                 (c1.r + 2 * d1.r) * beta.joint.est %*% t(beta.joint.est))
    ## Term 3
    term.3.var.r  <- (drop(t(beta.joint.est) %*% inv.r %*% beta.joint.est) *
                        sig.r + beta.joint.est %*% t(beta.joint.est)) / n.o
  }
  return(list('var.term.2' = term.2.var.r,
              'var.term.3' = term.3.var.r,
              'total'      = (term.2.var.r) + (term.3.var.r)))
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
  term.1.var.r <- varFirstTermNew(beta.mc, omega, x = x.r, ind = ind.beta.mc)  / nrow(x.r)
  ### Term 2
  term.2.var.r <- varSecondTermNew(beta.omega, x = x.r, ind = ind.beta.omega) / n.o
  ### Naive var
  return(list('var.term.1' = term.1.var.r,
              'var.term.2' = term.2.var.r,
              'total'      = (term.1.var.r + term.2.var.r)))

}
