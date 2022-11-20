library(MASS)
#### Create E matrix
createE <- function(n, i, j) {
  E      <- matrix(0, ncol = n, nrow = n)
  E[i,j] <- 1
  return(E)
}

#### Create permutation matrix
createK <- function(n) {
  zero.vec <- rep(0, n)
  K        <- matrix(0, ncol = n^2, nrow = n^2)
  for (i in 1:n) {
    for (j in 1:n) {
      E <- createE(n, i, j)
      K <- K + kronecker(E, t(E))
    }
  }
  return(K)
}

#### Create elimination matrix Md Vec(A) = V(diag(A))
createMd <- function(n) {
  Md <- 0
  for (i in 1:n) {
    E  <- (createE(n, i, i))
    Md <- Md + kronecker(E, E)
  }
  return(Md)
}

#### Create Ms
createMs <- function(n) {
  return(0.5 * (diag(n^2) + createK(n)))
}

forbeniusNorm <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  return(sqrt(sum((x - y)^2)))
}

### Parameters
n.est   <- 100000
n.org   <- 10000
n.ref   <- 5000
p       <- 5
maf.vec <- runif(p)
cov.mat <- 0.7^outer(1:p, 1:p, function(x, y) abs(y - x))
diag(cov.mat) <- sqrt(sqrt(rchisq(p, n.ref)))
eff.dim <- 1
h       <- 0.05
omega   <- solve(cov.mat)

############ Testing
x.ref       <- mvrnorm(n.ref, rep(0, p), cov.mat)
x.ref       <- scale(x.ref, scale = FALSE)
cov.x.ref   <- cov(x.ref)
cor.x.ref   <- cor(x.ref)
p           <- ncol(cov.x.ref)
### Helping matrices
k       <- createK(p)
Md      <- createMd(p)
Ms      <- createMs(p)
### Expression
A <- diag(p^2) - 0.5 * kronecker(cor.x.ref, diag(p)) %*% Md - 0.5 * kronecker(diag(p), cor.x.ref) %*% Md
first.term <- (diag(p^2) - Ms %*% kronecker(diag(p), cor.x.ref) %*% Md)





test_that("Verify that matrix P created correctly", {
  testthat::expect_true(forbeniusNorm(ECCCM:::findPij(est_cor_mat = cor(x.ref), 1, 1),
                                      A[(p + 1):(2 * p), (p + 1):(2 * p)]) <= 10^(-14))
  testthat::expect_true(forbeniusNorm(ECCCM:::findPij(est_cor_mat = cor(x.ref), 2, 1),
                                      A[(2 * p + 1):(3 * p), (p + 1):(2 * p)]) <= 10^(-14))
  testthat::expect_true(forbeniusNorm(ECCCM:::findPij(est_cor_mat = cor(x.ref), 1, 2),
                                      A[(p + 1):(2 * p), (2 * p + 1):(3 * p)]) <= 10^(-14))
})


test_that("Verify that matrix PLambda created correctly", {
  ### Plambda
  PLa <- A %*% kronecker(diag(1 / sqrt(diag(cov.x.ref))), diag(1 / sqrt(diag(cov.x.ref))))
  inv_sig  <- Matrix::sparseMatrix(x = 1 / sqrt(diag(cov.x.ref)), i = 1:p, j = 1:p, dims = c(p, p))

  testthat::expect_true(forbeniusNorm(ECCCM:::findPLambdaij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 1),
                                      PLa[(p + 1):(2 * p), (p + 1):(2 * p)]) <= 10^(-14))
  testthat::expect_true(forbeniusNorm(ECCCM:::findPLambdaij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 2, 1),
                                      PLa[(2 * p + 1):(3 * p), (p + 1):(2 * p)]) <= 10^(-14))
  testthat::expect_true(forbeniusNorm(ECCCM:::findPLambdaij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 2),
                                      PLa[(p + 1):(2 * p), (2 * p + 1):(3 * p)]) <= 10^(-14))
})


test_that("Verify that matrix LambdaP created correctly", {
  ### Plambda
  LaP <- kronecker(diag(1 / sqrt(diag(cov.x.ref))), diag(1 / sqrt(diag(cov.x.ref)))) %*% t(A)
  inv_sig  <- Matrix::sparseMatrix(x = 1 / sqrt(diag(cov.x.ref)), i = 1:p, j = 1:p, dims = c(p, p))
  testthat::expect_true(forbeniusNorm(ECCCM:::findLambdaPij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 1),
                                      LaP[(p + 1):(2 * p), (p + 1):(2 * p)]) <= 10^(-14))
  testthat::expect_true(forbeniusNorm(ECCCM:::findLambdaPij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 2, 1),
                                      LaP[(2 * p + 1):(3 * p), (p + 1):(2 * p)]) <= 10^(-14))
  testthat::expect_true(forbeniusNorm(ECCCM:::findLambdaPij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 2),
                                      LaP[(p + 1):(2 * p), (2 * p + 1):(3 * p)]) <= 10^(-14))
})

test_that("Verify covariance to correlation", {
  testthat::expect_true(forbeniusNorm(ECCCM:::covToCor(cov.mat), cov2cor(cov.mat)) <= 10^-14)
})

test_that("Verify extracting sparse diagonal (square root)", {
  temp.mat <- matrix(0, ncol(cov.mat), ncol(cov.mat))
  diag(temp.mat) <- diag(cov.mat)^(-0.5)
  testthat::expect_true(forbeniusNorm(temp.mat, temp.mat) < 10^-14)
})


### Variance of correlation matrix
test_that("Verify calculation of correlation matrix covarinace matrix", {

  cov.mat.list <- ECCCM:::outerProdRow(x.ref)
  var.cor.vec.mat <- var(t(replicate(1000,
                                     as.vector(cor(mvrnorm(n.ref, rep(0, p), cov.mat))), simplify = TRUE)))
  inv_sig  <- Matrix::sparseMatrix(x = 1 / sqrt(diag(cov.x.ref)), i = 1:p, j = 1:p, dims = c(p, p))

  testthat::expect_true(forbeniusNorm(
    ECCCM:::findRijV2(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 0, j = 0) / n.ref,
    var.cor.vec.mat[1:p, 1:p]
  ) < 10^-3)

  testthat::expect_true(forbeniusNorm(
    ECCCM:::findRijV2(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 1, j = 2) / n.ref,
    var.cor.vec.mat[(p + 1):(2 * p), (2 * p + 1):(3 * p)]
  ) < 10^-3)

  testthat::expect_true(forbeniusNorm(
    ECCCM:::findRijV2(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 2, j = 1) / n.ref,
    var.cor.vec.mat[(2 * p + 1):(3 * p), (p + 1):(2 * p)]
  ) < 10^-3)

})
