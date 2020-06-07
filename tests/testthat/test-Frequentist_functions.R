test_that("no rejection do not lead to error", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  testthat::expect_null(ECCCM::analyzeRef(rep(0, 10), x.r, 100)$var.omega.beta)
})


test_that("Function runs with various parameters", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  testthat::expect_true(class(ECCCM::analyzeRef(marg.beta.hat = rep(0.1, 10),
                                                x.r = x.r,
                                                n.o = 100,
                                                sigma.method = 'conservative',
                                                method.filter = 'none',
                                                method.test = 'BH',
                                                qu = 0.05)) == 'list')
})


test_that("Function runs with various parameters", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  testthat::expect_true(class(ECCCM::analyzeRef(marg.beta.hat = rep(0.1, 10),
                                                x.r = x.r,
                                                n.o = 100,
                                                sigma.method = 'estimate',
                                                method.filter = 'none',
                                                method.test = 'BH',
                                                qu = 0.05)) == 'list')
})

test_that("Function runs with various parameters", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  testthat::expect_true(class(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                                x.r = x.r,
                                                n.o = 100,
                                                sigma.method = 'conservative',
                                                method.filter = 'none',
                                                method.test = 'none',
                                                qu = 0.05)) == 'list')
})

test_that("Function runs with various parameters", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  testthat::expect_true(class(ECCCM::analyzeRef(marg.beta.hat = rep(0.1, 10),
                                                x.r = x.r,
                                                n.o = 100,
                                                sigma.method = 'semi.conservative',
                                                method.filter = 'BH',
                                                method.test = 'BH',
                                                qu = 0.05)) == 'list')
})




test_that("Function runs with various parameters", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  testthat::expect_true(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                                x.r = x.r,
                                                n.o = 100,
                                                sigma.method = 'conservative',
                                                method.filter = 'none',
                                                method.test = 'BH',
                                                qu = 0.05)$var.beta[1] > 1)
})



test_that("Functions works with all various options of regularization", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  testthat::expect_true(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                          x.r = x.r,
                                          n.o = 100,
                                          sigma.method = 'conservative',
                                          method.filter = 'none',
                                          method.test = 'BH',
                                          method.threshold = 'none',
                                          qu = 0.05)$var.beta[1] > 1)

  testthat::expect_true(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                          x.r = x.r,
                                          n.o = 100,
                                          sigma.method = 'conservative',
                                          method.filter = 'none',
                                          method.test = 'BH',
                                          method.threshold = 'soft',
                                          qu = 0.05)$var.beta[1] > 5)

  testthat::expect_true(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                          x.r = x.r,
                                          n.o = 100,
                                          sigma.method = 'conservative',
                                          method.filter = 'none',
                                          method.test = 'BH',
                                          method.threshold = 'hard',
                                          qu = 0.05)$var.beta[2] > 0.5)

  testthat::expect_error(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                           x.r = x.r,
                                           n.o = 100,
                                           sigma.method = 'conservative',
                                           method.filter = 'none',
                                           method.test = 'BH',
                                           method.threshold = 'nonsense',
                                           qu = 0.05))
})

test_that("analyzeRefGauss function runs with various parameters", {
  x.r <- matrix(rnorm(1000), ncol = 10, nrow = 100)
  cov.mat <- cov(x.r)
  testthat::expect_true(class(ECCCM::analyzeRefGauss(marg.beta.hat = rep(0.1, 10),
                                                     ld.mat = cov.mat,
                                                     n.o = 100,
                                                     n.r = 100,
                                                     sigma.method = 'semi.conservative',
                                                     method.filter = 'BH',
                                                     method.test = 'BH',
                                                     qu = 0.05)) == 'list')
})


test_that("analyzeRefGauss maintains the speficied FDR", {
  requireNamespace("MASS", quietly = TRUE)
  cov.mat <- (0.8^abs(outer(1:20, 1:20, '-')))
  emp.fdr <- rep(0, 1000)
  for (i in 1:1000) {
    x.r        <- scale(MASS::mvrnorm(300, rep(0, 20), cov.mat))
    x.o        <- scale(MASS::mvrnorm(900, rep(0, 20), cov.mat))
    y          <- x.o %*% c(rep(0, 19), 2) + rnorm(900)
    beta.marg  <- drop(t(x.o) %*% y) / nrow(x.o)
    beta.mult  <- nrow(x.r) * solve(t(x.r) %*% x.r) %*% beta.marg
    rej.vec    <- ECCCM::analyzeRefGauss(marg.beta.hat = beta.marg,
                                         ld.mat = cov(x.r),
                                         n.o = 900,
                                         n.r = 300,
                                         sigma.method = 'conservative',
                                         method.filter = 'none',
                                         method.test = 'BH',
                                         qu = 1)$test.correct[ ,3] < 0.05
    emp.fdr[i] <- sum(rej.vec[1:19]) / sum(rej.vec)
  }
  testthat::expect_true(mean(emp.fdr) < 0.05)
})


