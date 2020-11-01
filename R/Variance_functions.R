#' estimateVarAdd - estimation of additional variance
#' @param beta.mc The estimate of coefficient of joint regression (transformed marginal beta)
#' @param cov.list List of covariance matrices
#' @param cov.mat Estimated covarinace matrix
#' @param n.o Number of observation in original study
#' @param n.r Number of observation in reference study
#' @param ind.vec Indices of which `beta.mc` to consider
#' @param diag.flag - Flag, to estimate the variance of the covarinace matrix based on the main diagonal, or full estimation.
#' @return A list containing the varius variance terms, including the naive. `total` is the combined added variance from
#' using the reference panel
#' @export

estimateVarAdd <- function(beta.mc, cov.list, cov.mat, n.o, n.r, ind.vec, diag.flag = F) {
  ### Finding the variance of (beta' \prod \Sigma^{-1}) A (beta' prod \Sigma^{-1})
  if (!diag.flag) {
    temp.var     <- varBeta(beta     = beta.mc,
                            cov_mat  = cov.mat,
                            cov_list = cov.list,
                            ind = ind.vec,
                            nr  = n.r)
  }
  if (diag.flag) {
    temp.var     <- varBetaDiag(beta     = beta.mc,
                                cov_mat  = cov.mat,
                                cov_list = cov.list,
                                ind = ind.vec,
                                nr  = n.r)
  }
  ### Term 1
  term.1.var.r <- temp.var / n.r
  ### Term 2
  term.2.var.r <- temp.var / n.o
  ### Naive var
  return(list('var.term.1' = term.1.var.r,
              'var.term.2' = term.2.var.r,
              'total'      = (term.1.var.r + term.2.var.r)))

}






#' estimateVarAddGauss - estimation of additional variance
#' @param beta.mc The estimate of coefficient of joint regression (transformed marginal beta)
#' @param cov.list.boot Bootstrap of scaled covariance matrices.
#' @param ind.beta.mc Indices of which `beta.mc` to consider, in the bayesian scenario determined by vector gamma.
#' @param ind.beta.omega Indices of which `beta.omega` to consider
#' @param omega - Inverse of estimated covaraince matrix
#' @param n.o - Number of observation in original study
# #' @param diag.flag - Flag, to estimate the variance of the covarinace matrix based on the main diagonal, or full estimation.
#' @return A list containing the varius variance terms, including the naive. `total` is the combined added variance from
#' using the reference panel
#' @export

estimateVarAddGauss <- function(beta.mc, cov.mat, n.o, n.r, ind.vec) {
  ### Finding the variance of (beta' \prod \Sigma^{-1}) A (beta' prod \Sigma^{-1})
  temp.var     <- varBetaGauss(beta     = beta.mc,
                               cov_mat  = cov.mat,
                               ind = ind.vec,
                               nr  = n.r)
  ### Term 1
  term.1.var.r <- temp.var / n.r
  ### Term 2
  term.2.var.r <- temp.var / n.o
  ### Naive var
  return(list('var.term.1' = term.1.var.r,
              'var.term.2' = term.2.var.r,
              'total'      = (term.1.var.r + term.2.var.r)))

}






