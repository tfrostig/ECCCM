# Functions  --------------------------------------------------------------
### Input:  m   - scalar, number of SNPs
###         k   - scalar, maximum possible number of causal SNPs
### Output: p.k - scalar, probability of each v
priorK <- function(m, k) {
  p.k <- NA
  p.k <- choose(m, k) * ((1 / m)^k) * (((m - 1) / m)^(m - k)) ### Prior on number of causal SNPs
  return(p.k)
}
### Input:  z.hat   - vector of numeric marginal regressions
###         s.gamma - scalar, parameter of explained variance
###         R       - LD matrix
###         k       - number of casual SNPs
###         n       - scalar, number of observations
### Output: posteror.gamma - vector of posterior probability of SNPs to be causal
FineMap <- function(z.hat, s.gamma, R, k = 1, n) {
  m             <- length(z.hat)
  null.mat      <- matrix(0, nrow = m, ncol = m)
  comb.list     <- list()
  full.post.gamma <- NULL
  all.p.k       <- priorK(m, 1:k) ### K should be chosen carefuly.
  all.p.k       <- all.p.k / sum(all.p.k)
  for (ik in 1:k) {
    all.variation <- choose(m, ik)
    if (all.variation > 10000) {
      warning('number of variation exceeds 10,000, consider reduce m and/or k')
    }
    posterior.gamma <- rep(NA, m)
    p.k  <-  all.p.k[ik] ### K should be chosen carefuly.
    all.combinations <- utils::combn(m, ik) ## Enumerates all possible combinations
    comb.list[[ik]]  <- t(all.combinations)
    for (i in 1:all.variation) {
      gamma.vec         <- all.combinations[ ,i]
      delta.gamma       <- null.mat
      diag(delta.gamma)[gamma.vec] <- 1
      sigma.gamma       <- n * s.gamma^2 * delta.gamma
      posterior.gamma[i] <- (1 / all.variation) * p.k * dmvnorm(x      = z.hat,
                                                                mean   = rep(0, m),
                                                                sigma  = R + R %*% sigma.gamma %*% R)
    }
    full.post.gamma <- c(full.post.gamma, posterior.gamma)
  }
  return(cbind(Matrix::bdiag(comb.list), 'posterior' = full.post.gamma))
}

### Input:  z.hat   - vector of numeric marginal regressions
###         s.gamma - scalar, parameter of explained variance
###         x.ref   - genome reference
###         k       - maximum number of casual SNPs to consider
###         n.org   - scalar, number of observations in original research
### Output: posteror.gamma - vector of posterior probability of SNPs to be causal

FineMapCorrect <- function(z.hat, s.gamma, x.ref, k = 1, sigma = 1, n.org) {
  m             <- length(z.hat)
  n.ref         <- nrow(x.ref)
  null.mat      <- matrix(0, nrow = m, ncol = m)
  comb.list     <- list()
  full.post.gamma <- NULL
  all.p.k       <- priorK(m, 1:k) ### K should be chosen carefuly.
  all.p.k       <- all.p.k / sum(all.p.k)
  cov.mat.list  <- ECCM::outerProdRow(x.ref)
  cov.x.ref     <- cov(x.ref)
  omega.ref     <- solve(cov.x.ref)
  lambda.hat    <- n.org^(-1/2) * omega.ref %*% z.hat
  lambda.hat[sort(abs(lambda.hat), index.return = TRUE)$ix[1:(m - k - 1)]] <- 0 ### allow only k + 1 non zero coefficients (consider up to 3).
  A             <- estimateVariance(g.r = x.ref,
                                    sig.r = cov.x.ref,
                                    inv.r = omega.ref,
                                    beta.joint.est = lambda.hat,
                                    n.o = n.org,
                                    method = 'moment')
  added.var     <- A$var.term.2 + A$var.term.3
  for (ik in 1:k) {
    all.variation <- choose(m, ik)
    if (all.variation > 10000) {
      warning('number of variation exceeds 10,000, consider reduce m and/or k')
    }
    posterior.gamma <- rep(NA, m)
    p.k  <-  all.p.k[ik] ### K should be chosen carefuly.
    all.combinations <- utils::combn(m, ik) ## Enumerates all possible combinations
    comb.list[[ik]]  <- t(all.combinations)
    for (i in 1:all.variation) {
      add.var.z.hat <- n.org * (cov.x.ref %*% added.var %*% cov.x.ref)
      gamma.vec         <- all.combinations[ ,i]
      delta.gamma       <- null.mat
      diag(delta.gamma)[gamma.vec] <- 1
      sigma.gamma <- n.org * s.gamma^2 * delta.gamma
      posterior.gamma[i] <- (1 / all.variation) * p.k * dmvnorm(x      = drop(z.hat),
                                                                mean   = rep(0, m),
                                                                sigma  = cov.x.ref + add.var.z.hat + cov.x.ref %*% sigma.gamma %*% cov.x.ref)
    }
    full.post.gamma <- c(full.post.gamma, posterior.gamma)
  }
  return(cbind(Matrix::bdiag(comb.list), 'posterior' = full.post.gamma))
}
