
findOmegaBeta <- function(beta.mc, omega, sigma, n.o) {
  omega.beta       <-  cbind('omega_beta'     = omega %*% beta.mc,
                             'omega_beta_var' = (sigma^2 / n.o) * diag(omega %*% omega %*% omega))
}


#' estimateVarAdd - estimation of additional variance
#' @param beta.mc The estimate of coefficient of joint regression (transformed marginal beta)
#' @param beta.omega Transformed `beta.mc` basically, just `beta.omega %*% omega`, where omega is the estimated inverse covariance matrix.
#' @param cov.list.boot Bootstrap of scaled covariance matrices.
#' @param ind.beta.mc Indices of which `beta.mc` to consider, in the bayesian scenario determined by vector gamma.
#' @param ind.beta.omega Indices of which `beta.omega` to consider
#' @param omega - Inverse of estimated covaraince matrix
#' @param n.o - Number of observation in original study
#' @return A list containing the varius variance terms, including the naive. `total` is the combined added variance from
#' using the reference panel
#' @export

estimateVarAdd <- function(beta.mc, cov.list.boot, omega, n.o, n.r, ind.vec) {
  ### Finding the variance of (beta' \prod \Sigma^{-1}) A (beta' prod \Sigma^{-1})
  temp.var     <- varFirstTerm(beta.mc, omega, cov_list_scaled = cov.list.boot, ind = ind.vec)
  ### Term 1
  term.1.var.r <- temp.var / n.r
  ### Term 2
  term.2.var.r <- temp.var / n.o
  ### Naive var
  return(list('var.term.1' = term.1.var.r,
              'var.term.2' = term.2.var.r,
              'total'      = (term.1.var.r + term.2.var.r)))

}


#### Bootstrap functions
#' Create bootstrap for matrix
#' @param x matrix of n rows
#' @return returns a matrix with n rows sampled from x
createBootStrap <- function(x) {
  return(x[sample(nrow(x), nrow(x), replace = TRUE), ])
}

#### Create bootstrap covariance matrices
#' Create list of boostraped scaled covariance matrices
#' @param x matrix of n rows
#' @param B integer, specifying the number of bootstrap samples to create
#' @return returns a list of bootstrapped scaled covariance matrices
createListCov <- function(x, B = nrow(x)) {
  res.list <-  vector("list", length = B)
  for (i in 1:B) {
    res.list[[i]] <- cor(createBootStrap(x))
  }
  return(res.list)
}



