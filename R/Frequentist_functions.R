#' analyzeRef - Conduct full analysis based on reference panel
#' @param marg.beta.hat Marginal regression coeffcients, inner product of G_o and y divided by the number of observations in original study
#' @param x.r Reference panel
#' @param n.o The number of observations in original study
#' @param sigma.method String indicating how to estimate sigma, currently supports three methods, 'conservative', sigma is estimated as 1.
#' 'estimate' the usual estimator of the variance (assuming y is scaled sigma = 1 - var(x.r beta_mc)), and 'semi.conservative' which is similar to
#' the 'estimate' method, but the beta_mc is thresholded.
#' @param method.filter multiple adjustment method to use in thresholding of the coefficient (see p.adjust)
#' @param method.test multiple adjustment method  to use when testing the adjusted coefficients (see p.adjust)
#' @param qu P-value threshold.
#' @param explained.omega.beta Minimal weight for thresholding omega.beta
#' @param max.omega.beta Maximum number of cofficients taken from omega.beta
#' @param is.scaled Boolean flag, is original study X and y are scaled. For now only supports TRUE.
#' @return list containing `test.correct` data.frame with the adjusted coefficeint and adjusted variance, with two sided testing p-value
#' `test.naive` data.frame with the adjusted coefficeint but not adjusted variance, with two sided testing p-value,
#'  `add.var` the additional variance, `sigma` the estimated sigma, and
#'  `var.omega.beta` a vector with the number of coefficient used and the proportion of weight used.
#' @export


analyzeRef <- function(marg.beta.hat,
                       x.r,
                       n.o,
                       sigma.method         = 'conservative',
                       method.filter        = 'none',
                       method.test          = 'BH',
                       method.threshold     = 'soft',
                       qu                   = 0.05,
                       explained.omega.beta = 0.9,
                       max.omega.beta       = 20,
                       is.scaled            = TRUE) {
  p   <- ncol(x.r)
  n.r <- nrow(x.r)
  cov.list.r    <- covMatMaker(x.r,is.scaled)
  marg.to.joint <- marginalToJoint(marg.beta.hat = marg.beta.hat,
                                   x.r = x.r,
                                   n.o = n.o,
                                   cov.r = cov.list.r$cov,
                                   inv.r = cov.list.r$omega,
                                   sigma = 1)
  threshold.beta.est <- rep(0, p)
  ### Transform Beta to multivarite
  if (sigma.method == 'conservative') {
    test.df <- testCoef(est.beta = marg.to.joint$est.beta.hat,
                        var.beta = marg.to.joint$naive.var.beta.hat,
                        method   = method.filter)
    ind.beta.pass                     <- which(test.df[ ,3] < qu)
    sigma.est                         <- 1
    threshold.beta.est[ind.beta.pass] <-  marg.to.joint$est.beta.hat[ind.beta.pass]
  }
  if (sigma.method == 'estimate') {
    sigma.est         <- drop(sqrt(1 - var(x.r %*% marg.to.joint$est.beta.hat)))
    test.df           <- testCoef(est.beta = marg.to.joint$est.beta.hat,
                                  var.beta = marg.to.joint$naive.var.beta.hat * sigma.est,
                                  method   = method.filter)
    ind.beta.pass <- which(test.df[ ,3] < qu)
    threshold.beta.est[ind.beta.pass] <-  marg.to.joint$est.beta.hat[ind.beta.pass]
  }
  if (sigma.method == 'semi.conservative') {
    test.df <- testCoef(est.beta = marg.to.joint$est.beta.hat,
                        var.beta = marg.to.joint$naive.var.beta.hat,
                        method   = method.filter)
    ind.beta.pass                     <- which(test.df[ ,3] < qu)
    threshold.beta.est[ind.beta.pass] <-  marg.to.joint$est.beta.hat[ind.beta.pass]
    sigma.est         <- drop(sqrt(1 - var(x.r %*% threshold.beta.est)))
  }
  if (length(ind.beta.pass) == 0) {
    return(list(test.correct    = test.df,
                add.var         = NULL,
                test.naive      = test.df,
                var.omega.beta  = c('number.of.coef' = NULL,
                                    'prop.weight'    = NULL),
                var.beta        = c('number.of.coef' = NULL,
                                    'prop.weight'    = NULL),
                sigma.est           = sigma))
  }
  if (length(ind.beta.pass) > 0) {
    beta.omega                        <- cov.list.r$omega %*% threshold.beta.est
    weight.beta.omega                 <- findWeights(beta.omega,
                                                     var.explained = explained.omega.beta,
                                                     max.weights   = max.omega.beta)
    ind.betaomega.pass                <- weight.beta.omega$ind
    # weight.beta                       <- findWeights(test.df[ ,3],
    #                                                  var.explained = explained.beta,
    #                                                  max.weights   = max.beta)
    # ind.beta.pass                    <- weight.beta$ind ## removed due to being to conservative.
    est.var <- estimateVarAdd(beta.mc        = threshold.beta.est,
                              beta.omega     = beta.omega,
                              x.r            = x.r,
                              ind.beta.mc    = ind.beta.pass,
                              ind.beta.omega = ind.betaomega.pass,
                              omega          = cov.list.r$omega,
                              n.o            = n.o)

    test.correct.df <- testCoef(est.beta = threshold.beta.est,
                                var.beta = marg.to.joint$naive.var.beta.hat + diag(est.var$total),
                                method   = method.test)

    return(list(test.correct    = test.correct.df,
                add.var         = est.var$total,
                test.naive      = test.df,
                var.omega.beta  = c('number.of.coef' = weight.beta.omega$cut.off,
                                    'prop.weight'    = sum(beta.omega[ind.betaomega.pass]^2) / sum(beta.omega^2)),
                var.beta        = c('number.of.coef' = length(ind.beta.pass),
                                    'prop.weight'    = sum(test.df[ind.beta.pass,3]^2) / sum(test.df[ ,3]^2)),
                sigma           = sigma.est))
  }
}


#' testCoef - Two sided testing of coefficients
#' @param est.beta Estimated coefficients to test
#' @param var.beta The variance of coefficients
#' @param method  Method of multiplicity correction, adjusted according to
#' @return A list containing the varius variance terms, including the naive. `total` is the combined added variance from
#' using the reference panel
#' @export
#'
testCoef <- function(est.beta, var.beta, method = 'BH') {
  test.statistic <- est.beta / sqrt(var.beta)
  pval.no.adjust <- (2 * (1 - pnorm(abs(test.statistic))))
  pval.adjust    <- p.adjust(pval.no.adjust, method = method)
  res.mat        <- cbind('test.statistic'  = test.statistic,
                          'pval.unadjusted' = pval.no.adjust,
                          'temp'            = pval.adjust)
  colnames(res.mat)[3] <- paste('pval.adjusted', method, sep = '.')
  return(data.frame(res.mat))
}


#' covMatMakter - calculates covariance matrix
#' @param x.r Reference panel
#' @param is.scaled Flag, was the data scaled beforehand, currently package only support TRUE
#' @return list of covariance matrix and inverse covariance matrix
#' @export

covMatMaker <- function(x.r, is.scaled = TRUE) {
  n.r <- nrow(x.r)
  if (is.scaled) {
    cov.r    <- t(x.r) %*% x.r / (n.r - 1) ## Scaled, mean = 0.
    inv.r    <- (n.r - 1) * solve(t(x.r) %*% x.r)
  }
  if (!is.scaled) {
    cov.r    <- cov(x.r)
    inv.r    <- solve(cov(x.r))
  }
  return(list('cov'   = cov.r,
              'omega' = inv.r))
}

#' marginalToJoint - Transforms maringal coefficients estimates to joint coefficients estimates
#' @param marg.beta.hat Estiamted coefficeint resulting from marginal regression
#' @param x.r Reference panel
#' @param cov.r Estimated covariance of x.r
#' @param inv.r Inverse of cov.r
#' @param n.o Number of observation in original study
#' @param sigma Estimated variance of epsilon, since y is assumed to be normalized then defaults to 1
#' @param is.scaled Boolean flag, is original study X and y are scaled. For now only supports TRUE.
#' @return data.frame containing column `est.beta.hat` estimate of beta, and `naive.var.beta.hat`
#' the variance of the estimate not taking into account the use of the reference panel
#' @export
marginalToJoint <- function(marg.beta.hat, x.r, cov.r, inv.r, n.o,
                            sigma = 1, is.scaled = TRUE) {
  if (is.scaled) {
    est.beta <- drop(inv.r %*% marg.beta.hat)
  }
  if (!is.scaled) {
    norm.x.r <- diag(apply(x.r, 2, var))
    est.beta <- drop(norm.x.r %*% inv.r %*% marg.beta.hat)
  }
  var.beta    <- sigma^2 * diag(inv.r) / n.o
  return(data.frame('est.beta.hat'       = est.beta,
                    'naive.var.beta.hat' = var.beta))
}



