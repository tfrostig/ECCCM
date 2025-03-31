#### Sanity check everything works
library(ECCCM)
library(MASS)

calcFDR <- function(hyp.rej, true.rej, hypo.num){
  if (length(hyp.rej) == 0) {
    return(0)
  } else {
    return(length(intersect(hyp.rej, setdiff(1:hypo.num, true.rej))) / length(hyp.rej))
  }
}


### WARNING: uses hard coded methods for filter and estimaton of sigma
estimateAdditionalVariance <- function(coeff, x.ref, n.org, qu, cor.list, is.marg) {
  ### Dealing with case where all coefficients are thesholded to 0 (no added)
  add.var <- 0
  ### Correcting covariance matrix
  test.emp.all.est   <-  ECCCM:::analyzeRef(marg.beta.hat = coeff,
                                            x.r           = x.ref,
                                            n.o           = n.org,
                                            method.filter = 'BH',
                                            sigma.method  = 'estimate',
                                            qu            = qu,
                                            to.diag       = FALSE,
                                            is.marg       = is.marg)

  if (!is.null(test.emp.all.est$add.var)) {
    add.var <- test.emp.all.est$add.var
  }
  correct.var <- test.emp.all.est$sigma^2 * (cor.list$cor_omega / n.org) + add.var
  return(list('var'   = correct.var,
              'sigma' = test.emp.all.est$sigma))
}



coef.vec <- c(0, 0, 1 ,0, 1, 0, 0, 0, 0, 0)
n.o <- 80000
n.r <- 1000
qu  <- 0.05

expected.qu <- qu * sum(coef.vec == 0) / length(coef.vec)

cov.mat <- matrix(0.7, ncol = length(coef.vec), nrow =  length(coef.vec))
diag(cov.mat) <- 1

temp.vec <- NULL
temp.orc <- NULL
temp.nai <- NULL
temp.var.orc <- NULL
cas.snp <- which(coef.vec != 0)

for (i in 1:100) {
  #x   <- matrix(rnorm((n.o + n.r) * 5), ncol = 5)
  x   <-  mvrnorm(n.o + n.r, rep(0, ncol(cov.mat)), cov.mat)
  x.o <- x[1:n.o, ]
  x.r <- x[(n.o + 1):(n.o + n.r), ]

  y           <- x %*% coef.vec + rnorm(n.o + n.r, 0, 11)
  sd.y.oracle <- as.numeric(sqrt(t(coef.vec) %*% cov.mat %*% coef.vec + 11^2))
  sd.y        <- sd(y)

  x.o <- scale(x.o)
  y.o <- scale(y[1:n.o])
  y.r <- scale(y[(n.o + 1):(n.o + n.r)])

  beta.org    <- solve(t(x.o) %*% x.o) %*% t(x.o) %*% y.o
  beta.report <- drop(solve(diag(apply(x.o, 2, var))) %*% t(x.o) %*% y.o / n.o)

  #### Oracle
  oracle.lm   <- summary(lm(y.o ~ x.o))
  pval.oracle <- p.adjust(oracle.lm$coefficients[2:6, 4], 'BH')

  #### Naive
  marg.to.joint <- ECCCM::marginalToJoint(marg.beta.hat = beta.report,
                                          n.o = n.o,
                                          cor.r = cor(x.r),
                                          sigma = 1)

  cor.list.oracle  <- list('cor'       = cov2cor(cov.mat),
                           'cor_omega' = solve(cov2cor(cov.mat)))

  ### oracle estimator variance
  correct.var.oracle <- ECCCM:::estimateVarAddGauss(coef.vec / sd.y.oracle,
                                                    cov.mat = cov.mat,
                                                    n.o = n.o,
                                                    n.r = n.r,
                                                    which(coef.vec != 0))$total ## which(coef.vec != 0)

  test.naive.df      <- ECCCM::testCoef(est.beta = marg.to.joint$est.beta.hat,
                                   var.beta = marg.to.joint$naive.var.beta.hat,
                                   method   = 'BH')

  test.var.oracle.df <- ECCCM::testCoef(est.beta     = marg.to.joint$est.beta.hat,
                                        var.beta     = marg.to.joint$naive.var.beta.hat + diag(correct.var.oracle),
                                        method       = 'BH')

  pval.method        <- ECCCM:::analyzeRefGauss(marg.beta.hat = beta.report,
                                                ld.mat        = cov(x.r),
                                                n.o           = n.o,
                                                n.r           = n.r,
                                                method.filter = 'BH',
                                                sigma.method  = 'estimate',
                                                qu            = qu)

  temp.naive       <- test.naive.df$pval.adjusted.BH
  temp.correct     <- pval.method$test.correct[ ,4]
  temp.var.oracle  <- test.var.oracle.df$pval.adjusted.BH


  temp.vec[i]     <- calcFDR(which(temp.correct < qu), cas.snp, length(coef.vec))
  temp.orc[i]     <- calcFDR(which(pval.oracle < qu), cas.snp, length(coef.vec))
  temp.nai[i]     <- calcFDR(which(temp.naive < qu), cas.snp, length(coef.vec))
  temp.var.orc[i] <- calcFDR(which(temp.var.oracle < qu), cas.snp, length(coef.vec))
}


mean(temp.vec)
mean(temp.orc)
mean(temp.nai)
mean(temp.var.orc)

