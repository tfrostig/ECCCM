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
                                                qu = 0.05,
                                                explained.omega.beta = 0.5)$var.omega.beta[2] > 0.5)
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
                                          qu = 0.05,
                                          explained.omega.beta = 0.5)$var.omega.beta[2] > 0.5)

  testthat::expect_true(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                          x.r = x.r,
                                          n.o = 100,
                                          sigma.method = 'conservative',
                                          method.filter = 'none',
                                          method.test = 'BH',
                                          method.threshold = 'soft',
                                          qu = 0.05,
                                          explained.omega.beta = 0.5)$var.omega.beta[2] > 0.5)

  testthat::expect_true(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                          x.r = x.r,
                                          n.o = 100,
                                          sigma.method = 'conservative',
                                          method.filter = 'none',
                                          method.test = 'BH',
                                          method.threshold = 'hard',
                                          qu = 0.05,
                                          explained.omega.beta = 0.5)$var.omega.beta[2] > 0.5)

  testthat::expect_error(ECCCM::analyzeRef(marg.beta.hat = rep(1, 10),
                                           x.r = x.r,
                                           n.o = 100,
                                           sigma.method = 'conservative',
                                           method.filter = 'none',
                                           method.test = 'BH',
                                           method.threshold = 'nonsense',
                                           qu = 0.05,
                                           explained.omega.beta = 0.5))
})

