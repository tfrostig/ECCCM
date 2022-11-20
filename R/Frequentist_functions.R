#' analyzeRef - Conduct full analysis based on reference panel
#' @param marg.beta.hat Marginal regression coeffcients, inner product of G_o and y divided by the number of observations in original study
#' @param x.r Reference panel must be centered
#' @param n.o The number of observations in original study
#' @param sigma.method String indicating how to estimate sigma, currently supports three methods, 'conservative', sigma is estimated as 1.
#' 'estimate' the usual estimator of the variance (assuming y is scaled sigma = 1 - var(x.r beta_mc)), and 'semi.conservative' which is similar to
#' the 'estimate' method, but the beta_mc is thresholded.
#' @param method.filter multiple adjustment method to use in thresholding of the coefficient (see p.adjust)
#' @param method.test multiple adjustment method  to use when testing the adjusted coefficients (see p.adjust)
#' @param threshold.cov.seq user-defined threshold value. If it is a vector of regularization values, it automatically selects one that minimizes cross validation risk.
#' @param qu P-value threshold.
# @param is.scaled Boolean flag, is original study X and y are scaled. For now only supports TRUE.
#' @return list containing `test.correct` data.frame with the adjusted coefficeint and adjusted variance, with two sided testing p-value.
#' `test.naive` data.frame with the adjusted coefficeint but not adjusted variance, with two sided testing p-value.
#'  `add.var` the additional variance, `sigma` the estimated sigma.
#'  `var.omega.beta` a vector with the number of coefficient used and the proportion of weight used.
#' @export

analyzeRef <- function(marg.beta.hat,
                       x.r,
                       n.o,
                       sigma.method         = 'conservative',
                       method.filter        = 'none',
                       method.test          = 'BH',
                       qu                   = 0.05,
                       to.diag              = FALSE,
                       is.marg              = TRUE) {
  p   <- ncol(x.r)
  n.r <- nrow(x.r)
  cov.list.r  <- list('cov' = cov(x.r), 'cor' = cor(x.r))
  if (is.marg) {
    print('Transforming marginal coefficient into conditional')
    marg.to.joint <- marginalToJoint(marg.beta.hat = marg.beta.hat,
                                     n.o = n.o,
                                     cor.r = cov.list.r$cor,
                                     sigma = 1)
  }
  if (!is.marg) {
    marg.to.joint <- data.frame('est.beta.hat'       = marg.beta.hat,
                                'naive.var.beta.hat' = diag((1 / n.o) * solve(cov.list.r$cor)))
  }
  #### estimating variance and threshold beta
  threshold.result <- findSigmaThreshold(multi.beta.hat = marg.to.joint$est.beta.hat,
                                         multi.beta.hat.var = marg.to.joint$naive.var.beta.hat,
                                         x.r = x.r,
                                         qu = qu, method.filter = method.filter, sigma.method = sigma.method)


  #### dealing with case where non of the coefficients passed
  if (length(threshold.result$ind.beta.est) == 0) {
    return(list(test.correct    = threshold.result$test.df,
                add.var         = NULL,
                test.naive      = threshold.result$test.df,
                var.beta        = c('number.of.coef' = NULL,
                                    'prop.weight'    = NULL),
                sigma.est       = threshold.result$sigma.est))
  }
  if (length(threshold.result$ind.beta.est) > 0) {
    ### estimating additional variance
    cov.list           <- outerProdRow(x.r)
    est.var            <- estimateVarAdd(beta.mc        = threshold.result$threshold.beta,
                                         cov.list       = cov.list,
                                         cov.mat        = cov.list.r$cov,
                                         n.o            = n.o,
                                         n.r            = n.r,
                                         ind.vec        = threshold.result$ind.beta.est,
                                         diag.flag      = to.diag)
    ### corrected testing of coeffcients
    added.var       <- threshold.result$sigma.est^2 * marg.to.joint$naive.var.beta.hat + diag(est.var$total)
    test.correct.df <- testCoef(est.beta = threshold.result$threshold.beta,
                                var.beta = added.var,
                                method   = method.test)
    ### weight of non-thresholded coefficient out of all original
    prop <- sum(threshold.result$test.df[threshold.result$ind.beta.est, 3]^2) / sum(threshold.result$test.df[ ,3]^2)
    return(list(test.correct                = test.correct.df,
                add.var                     = est.var$total,
                test.naive                  = threshold.result$test.df,
                var.beta                    = c('number.of.coef' = length(threshold.result$ind.beta.pass),
                                                'prop.weight'    = prop),
                sigma                       = threshold.result$sigma.est))
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
  res.mat        <- cbind('beta'            = est.beta,
                          'test.statistic'  = test.statistic,
                          'pval.unadjusted' = pval.no.adjust,
                          'pval.adjusted'   = pval.adjust)
  colnames(res.mat)[4] <- paste('pval.adjusted', method, sep = '.')
  return(data.frame(res.mat))
}


#' covMatMakter - calculates covariance matrix
#' @param x.r Reference panel
#' @param is.scaled Flag, was the data scaled beforehand, currently package only support TRUE
#' @param method.threshold which type of threshold to apply on the covariance matrix, currently supports `soft`, `hard` and `none`. Default to `soft`.
#' @param cv.threshold the number of repetitions for 2-fold random cross validations for each threshold value.
#' @return list of covariance matrix and inverse covariance matrix
#' @export

covMatMaker <- function(x.r, method.threshold = 'none', cv.threshold, threshold.seq) {
  n.r <- nrow(x.r)
  if (method.threshold == 'soft') {
    threshold.para <- CovTools::CovEst.soft(x.r, thr = threshold.seq, nCV = cv.threshold)
    cov.r          <- threshold.para$S
    inv.r          <- solve(cov.r)
  }
  if (method.threshold == 'hard') {
    threshold.para <- CovTools::CovEst.hard(x.r, thr = threshold.seq, nCV = cv.threshold)
    cov.r          <- threshold.para$S
    inv.r          <- solve(cov.r)
  }
  if (method.threshold == 'none') {
    threshold.para <- NULL
    cov.r          <- cov(x.r)
    inv.r          <- solve(cov.r)
  }
  if (!(method.threshold %in% c('soft', 'hard', 'none'))) {
    stop('method.threshold must be one of soft/hard/none')
  }
  return(list('cov'       = cov.r,
              'omega'     = inv.r,
              'threshold' = threshold.para))
}

#' marginalToJoint - Transforms marginal coefficients estimates to joint coefficients estimates
#' @param marg.beta.hat Estimated coefficients resulting from marginal regression
#' @param x.r Reference panel
#' @param cov.r Estimated covariance of x.r
#' @param inv.r Inverse of cov.r
#' @param n.o Number of observation in original study
#' @param sigma Estimated variance of epsilon, since y is assumed to be normalized then defaults to 1
#' @return data.frame containing column `est.beta.hat` estimate of beta, and `naive.var.beta.hat`
#' the variance of the estimate not taking into account the use of the reference panel
#' @export
marginalToJoint <- function(marg.beta.hat,
                            cor.r,
                            n.o,
                            sigma = 1) {
  inv.r    <- solve(cor.r)
  est.beta <- drop(inv.r %*% marg.beta.hat)
  var.beta <- sigma^2 * diag(inv.r) / n.o
  return(data.frame('est.beta.hat'       = est.beta,
                    'naive.var.beta.hat' = var.beta))
}


#' analyzeRefGauss - Conduct full analysis based on reference panel under the assumption of Gaussian reference
#' @param marg.beta.hat Marginal regression coefficients, inner product of G_o and y divided by the number of observations in original study
#' @param ld.mat Covariance matrix based on the reference panel
#' @param n.r The number of observation in reference panel
#' @param n.o The number of observations in original study
#' @param sigma.method String indicating how to estimate sigma, currently supports three methods, 'conservative', sigma is estimated as 1.
#' 'estimate' the usual estimator of the variance (assuming y is scaled sigma = 1 - var(x.r beta_mc)), and 'semi.conservative' which is similar to
#' the 'estimate' method, but the beta_mc is thresholded.
#' @param method.filter multiple adjustment method to use in thresholding of the coefficient (see p.adjust)
#' @param method.test multiple adjustment method  to use when testing the adjusted coefficients (see p.adjust)
#' @param qu Coefficients with p-value lower than 'qu' are thresholded to 0.
#' @return list containing `test.correct` data.frame with the adjusted coefficient and adjusted variance, with two sided testing p-value.
#' `test.naive` data.frame with the adjusted coefficient but not adjusted variance, with two sided testing p-value.
#' `add.var` the additional variance, `sigma` the estimated sigma.
#' `var.omega.beta` a vector with the number of coefficient used and the proportion of weight used.
analyzeRefGauss<- function(marg.beta.hat,
                           ld.mat,
                           n.o,
                           n.r,
                           sigma.method         = 'conservative',
                           method.filter        = 'none',
                           method.test          = 'BH',
                           qu                   = 1) {
  p             <- ncol(ld.mat)
  cov.list.r    <- list('cov' = ld.mat, 'omega' = solve(ld.mat))
  marg.to.joint <- marginalToJoint(marg.beta.hat = marg.beta.hat,
                                   n.o = n.o,
                                   cor.r = cov2cor(cov.list.r$cov),
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
    sigma.est         <- drop(sqrt(1 - marg.to.joint$est.beta.hat %*% cov.list.r$cov %*% marg.to.joint$est.beta.hat))
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
    sigma.est         <- drop(sqrt(1 - marg.to.joint$est.beta.hat %*% cov.list.r$cov %*% marg.to.joint$est.beta.hat))
  }
  if (length(ind.beta.pass) == 0) {
    return(list(test.correct    = test.df,
                add.var         = NULL,
                test.naive      = test.df,
                sigma.est       = sigma.est))
  }
  if (length(ind.beta.pass) > 0) {
    ### Applying the shrinkage of beta
    est.var <- estimateVarAddGauss(beta.mc     = ((n.r) / (n.r - p - 1)) * threshold.beta.est,
                                   cov.mat     = cov.list.r$cov,
                                   n.o         = n.o,
                                   n.r         = n.r,
                                   ind.vec     = ind.beta.pass)

    test.correct.df <- testCoef(est.beta = threshold.beta.est,
                                var.beta = marg.to.joint$naive.var.beta.hat + diag(est.var$total),
                                method   = method.test)

    return(list(test.correct                = test.correct.df,
                add.var                     = est.var,
                test.naive                  = test.df,
                sigma                       = sigma.est))
  }
}



### find variance according to specification
findSigmaThreshold <- function(multi.beta.hat, multi.beta.hat.var, x.r, qu, method.filter, sigma.method) {
  threshold.beta.est <- rep(0, length(multi.beta.hat))
  ### Transform Beta to multivariate
  if (sigma.method == 'conservative') {
    test.df <- testCoef(est.beta = multi.beta.hat,
                        var.beta = multi.beta.hat.var,
                        method   = method.filter)
    ind.beta.pass                     <- which(test.df[ ,3] < qu)
    sigma.est                         <- 1
    threshold.beta.est[ind.beta.pass] <-  multi.beta.hat[ind.beta.pass]
  }
  if (sigma.method == 'estimate') {
    sigma.est         <- drop(sqrt(1 - var(x.r %*% multi.beta.hat)))
    test.df           <- testCoef(est.beta = multi.beta.hat,
                                  var.beta = multi.beta.hat.var * sigma.est,
                                  method   = method.filter)
    ind.beta.pass <- which(test.df[ ,3] < qu)
    threshold.beta.est[ind.beta.pass] <-  multi.beta.hat[ind.beta.pass]
  }
  if (sigma.method == 'semi.conservative') {
    test.df <- testCoef(est.beta = multi.beta.hat,
                        var.beta = multi.beta.hat.var,
                        method   = method.filter)
    ind.beta.pass                     <- which(test.df[ ,3] < qu)
    threshold.beta.est[ind.beta.pass] <-  multi.beta.hat[ind.beta.pass]
    sigma.est         <- drop(sqrt(1 - var(x.r %*% threshold.beta.est)))
  }
  return(list('test.df'        = test.df,
              'threshold.beta' = threshold.beta.est,
              'ind.beta.est'   = ind.beta.pass,
              'sigma.est'      = sigma.est))
}

