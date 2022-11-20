### Checking computation time of analyzeRef
library(MASS)
library(tictoc)
### Parameters
n.est   <- 10000
n.org   <- 10000
n.ref   <- 2700
p       <- 100
maf.vec <- runif(p)
cov.mat <- 0.7^outer(1:p, 1:p, function(x, y) abs(y - x))
diag(cov.mat) <- 1 + sample(rpois(p, 5))
eff.dim <- 1
h       <- 0.05
omega   <- solve(cov.mat)


############ Testing
x.ref       <- mvrnorm(n.ref, rep(0, p), cov.mat)
x.ref       <- scale(x.ref, scale = FALSE)
cov.x.ref   <- (t(x.ref) %*% x.ref) / n.ref
cor.x.ref   <- cor(x.ref)
p           <- ncol(cov.x.ref)
beta.vec    <- c(1, 1, rep(0, p - 2))
cor.mat     <- cov2cor(cov.mat)
omega.cor    <- solve(cor.mat)


microbenchmark::microbenchmark(varBetaGauss(c(c(1, 1), rep(0, p - 2)), cov.mat, c(1, 2), 100), times = 10)


### Computation of estimated variance of covariance matrix
microbenchmark::microbenchmark(cov.mat.list <- outerProdRow(x.ref))
microbenchmark::microbenchmark(cov.var <- findCovVar(cov.mat.list, cov.x.ref))

### Computation
microbenchmark::microbenchmark((varBeta(c(c(1, 1), rep(0, p - 2)), cov.mat, cov.mat.list, c(1, 2), 3)) / n.ref)
microbenchmark::microbenchmark(varBetaGauss(c(c(1, 1), rep(0, p - 2)), cov.mat, c(1, 2), 100))


microbenchmark::microbenchmark(findRij(est_cor_mat = cov2cor(cov.mat), inv_d_sig = inv_sig,
                                       cov_list = cov.mat.list, i = 0, j = 0))

microbenchmark::microbenchmark(ECCCM:::findCovTwoInd(cov.mat.list, 0, 0))
microbenchmark::microbenchmark(findPLambdaij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 1))
