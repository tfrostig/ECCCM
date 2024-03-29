// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

/// Finding the covariance of two rows of covariance matrix. /arma::mat
//estCov(arma::mat g_r, arma::mat cov_est, int i, int j) {

// [[Rcpp::export]]
List outerProdRow(arma::mat X) {
  int n = X.n_rows;
  List res_list = List(n);
  for (int i = 0; i < n; i++) {
    res_list[i] = X.row(i).t() * X.row(i);
  }
  return res_list;
}

// [[Rcpp::export]]
arma::mat findCovVar(List list_cov_mat, arma::mat cov_mat) {
  int nl = list_cov_mat.size();
  int i;
  arma::mat examp_mat = list_cov_mat[0];
  int p = examp_mat.n_cols;
  arma::mat sol = arma::mat(pow(p, 2), pow(p, 2), fill::zeros);
  arma::mat mean_cov_mat = vectorise(cov_mat) * vectorise(cov_mat).t();
  for (i = 0; i < nl; i++) {
    arma::mat temp_mat = list_cov_mat[i];
    arma::vec temp_vec = vectorise(temp_mat);
    sol += temp_vec * temp_vec.t();
  }
  sol = (sol / nl) - mean_cov_mat;
  return sol;
}

// [[Rcpp::export]]
arma::mat findCovTwoInd(List cov_list, int ind_a, int ind_b) {
  int n = cov_list.size();
  arma::mat temp_mat = cov_list[0];
  int p = temp_mat.n_cols;
  arma::mat sol = arma::mat(p, p, fill::zeros);
  arma::vec sol_a(p, fill::zeros);
  arma::vec sol_b(p, fill::zeros);
  for (int i = 0; i < n; i++) {
    arma::mat temp_mat = cov_list[i];
    sol   += temp_mat.col(ind_a) * temp_mat.col(ind_b).t();
    sol_a += temp_mat.col(ind_a);
    sol_b += temp_mat.col(ind_b);
  }
  arma::mat temp_sol = (sol / n) - (sol_a / n) * (sol_b / n).t();
  return temp_sol;
}


// [[Rcpp::export]]
arma::mat findCovTwoIndV2(arma::cube outer_cube, arma::mat cov_mat, int ind_a, int ind_b) {
  int n = outer_cube.n_elem_slice;
  arma::mat col_ind_a = outer_cube.col(ind_a);
  arma::mat col_ind_b = outer_cube.col(ind_b);
  return (col_ind_a * col_ind_b.t() / n) - cov_mat.col(ind_a) * cov_mat.col(ind_b).t();
}


// [[Rcpp::export]]
double quadForm(arma::vec x, arma::mat S) {
  double sol = as_scalar(x.t() * S * x);
  return sol;
}

// [[Rcpp::export]]
arma::vec powVec(arma::vec x, double power_of) {
  int n = x.size();
  for (int i = 0; i < n; i++) {
    x[i] = std::pow(x[i], power_of);
  }
  return x;
}

// [[Rcpp::export]]
arma::mat addVarGauss(arma::mat cov_mat, arma::mat omega, arma::vec beta, int n_r, int n_o) {
  int p = cov_mat.n_cols;
  double c_const = pow((n_r - p) * (n_r - p - 1) * (n_r - p - 3), -1);
  double d_const = pow((n_r - p) * pow(n_r - p - 1, 2) * (n_r - p - 3), -1);
  arma::mat beta_outer = beta * beta.t();
  arma::mat sol_a = (quadForm(beta, cov_mat) * omega + beta_outer) / n_o;
  arma::mat sol_b = pow(n_r, 2)* (c_const * quadForm(beta, cov_mat) * omega + (c_const + 2 * d_const) * beta_outer);
  arma::mat sol_c = ((n_r - 1) * pow(n_r, 2) * c_const / n_o) * omega;
  return(sol_a + sol_b + sol_c);
}


// [[Rcpp::export]]
arma::sp_mat createE(int dim_size, int i, int j) {
  arma::sp_mat E(dim_size, dim_size);
  E(i,j) = 1;
  return E;
}

// [[Rcpp::export]]
arma::sp_mat findPij(arma::mat est_cor_mat, int first_ind, int second_ind) {
  int p              = est_cor_mat.n_cols;
  arma::mat res      = arma::zeros(p, p);
  if (first_ind == second_ind) {
    res.col(first_ind) = -0.5 * est_cor_mat.col(first_ind);
    res(first_ind, first_ind) += -0.5 * est_cor_mat(first_ind, first_ind);
    res.diag() += 1;
  }
  if (first_ind != second_ind) {
    res(second_ind, second_ind) = -0.5 * est_cor_mat(first_ind, second_ind);
  }
  arma::sp_mat res_sparse(res);
  return res_sparse;
}

// [[Rcpp::export]]
arma::sp_mat findPLambdaij(arma::mat est_cor_mat, arma::sp_mat inv_d_sig, int first_ind, int second_ind) {
  int p              = est_cor_mat.n_cols;
  arma::sp_mat res(p, p);
  if (first_ind == second_ind) {
    res.col(first_ind) = -0.5 * est_cor_mat.col(first_ind);
    res(first_ind, first_ind) += -0.5 * est_cor_mat(first_ind, first_ind);
    res.diag() += 1;
  }
  if (first_ind != second_ind) {
    res(second_ind, second_ind) = -0.5 * est_cor_mat(first_ind, second_ind);
  }
  res = res * inv_d_sig * inv_d_sig(second_ind, second_ind);
  return res;
}

// [[Rcpp::export]]
arma::sp_mat findLambdaPij(arma::mat est_cor_mat, arma::sp_mat inv_d_sig, int first_ind, int second_ind) {
  return findPLambdaij(est_cor_mat, inv_d_sig, second_ind, first_ind).t();
}

// [[Rcpp::export]]
arma::mat findRij(arma::mat est_cor_mat, arma::sp_mat inv_d_sig, List cov_list, int i, int j) {
  int p = est_cor_mat.n_cols;
  arma::mat res = arma::zeros(p, p);
  for (int k = 0; k < p; k++) {
    for (int l = 0; l < p; l++) {
      res += findPLambdaij(est_cor_mat, inv_d_sig, i, l) * findCovTwoInd(cov_list, l, k) * findLambdaPij(est_cor_mat, inv_d_sig, k, j);
    }
  }
  return res;
}


// [[Rcpp::export]]
arma::mat findVarGaussCovInd(arma::mat cov_mat, int i, int j) {
  arma::mat outer_mat = cov_mat.col(j) * cov_mat.col(i).t();
  arma::mat sol = arma::as_scalar(cov_mat(j, i)) * cov_mat + outer_mat;
  return sol;
}


// [[Rcpp::export]]
arma::mat findRijGauss(arma::mat est_cor_mat, arma::sp_mat inv_d_sig, arma::mat est_cov_mat, int i, int j) {
  int p = est_cor_mat.n_cols;
  arma::mat res = arma::zeros(p, p);
  for (int k = 0; k < p; k++) {
    for (int l = 0; l < p; l++) {
      res += findPLambdaij(est_cor_mat, inv_d_sig, i, l) * findVarGaussCovInd(est_cov_mat, l, k) * findLambdaPij(est_cor_mat, inv_d_sig, k, j);
    }
  }
  return res;
}

// [[Rcpp::export]]
arma::mat findRijDiagonal(arma::mat est_cor_mat, arma::sp_mat inv_d_sig, List cov_list, int i, int j) {
  int p = est_cor_mat.n_cols;
  arma::mat res = arma::zeros(p, p);
  for (int k = 0; k < p; k++) {
      res += findPLambdaij(est_cor_mat, inv_d_sig, i, k) * findCovTwoInd(cov_list, k, k) * findLambdaPij(est_cor_mat, inv_d_sig, k, j);
    }
  return res;
}


// [[Rcpp::export]]
arma::mat covToCor(arma::mat cov_mat) {
  arma::vec sd_vec = powVec(cov_mat.diag(), -0.5);
  arma::mat cor_mat = cov_mat % (sd_vec * sd_vec.t());
  return cor_mat;
}

// [[Rcpp::export]]
arma::sp_mat diagSqrtSparse(arma::mat temp_mat) {
  int p = temp_mat.n_cols;
  arma::sp_mat sol(p, p);
  arma::vec sd_vec = powVec(temp_mat.diag(), -0.5);
  sol.diag() = sd_vec;
  return sol;
}

// [[Rcpp::export]]
arma::mat varBeta(arma::vec beta, arma::mat cov_mat, List cov_list, arma::vec ind) {
  int p = beta.size();
  int k = ind.size();
  int i_ind = 0;
  int j_ind = 0;
  ind = ind - 1; // indices are from R
  arma::mat sol          = arma::mat(p, p, fill::zeros);
  arma::mat cor_mat      = covToCor(cov_mat);
  arma::sp_mat inv_d_sig = diagSqrtSparse(cov_mat);
  arma::mat omega        = arma::inv_sympd(cor_mat);
  // omega = omega.clean(10^-8); // thresholded omega
  for (int i = 0; i < k; i++) {
    i_ind = ind.at(i);
    for (int j = 0; j < k; j++) {
      j_ind = ind.at(j);
      sol  +=  beta.at(i_ind) * beta.at(j_ind) * findRij(cor_mat, inv_d_sig, cov_list, i_ind, j_ind);
    }
  }
  return omega * sol * omega;
}



// [[Rcpp::export]]
arma::mat varBetaGauss(arma::vec beta, arma::mat cov_mat, arma::vec ind) {
  int p = beta.size();
  int k = ind.size();
  int i_ind = 0;
  int j_ind = 0;
  ind = ind - 1; // indices are from R
  arma::mat sol          = arma::mat(p, p, fill::zeros);
  arma::mat cor_mat      = covToCor(cov_mat);
  arma::sp_mat inv_d_sig = diagSqrtSparse(cov_mat);
  arma::mat omega        = arma::inv_sympd(cor_mat);
  //omega = omega.clean(10^-8); // thresholded omega
  for (int i = 0; i < k; i++) {
    i_ind = ind.at(i);
    for (int j = 0; j < k; j++) {
      j_ind = ind.at(j);
      sol  +=  beta.at(i_ind) * beta.at(j_ind) * findRijGauss(cor_mat, inv_d_sig, cov_mat, i_ind, j_ind);
    }
  }
  return omega * sol * omega;
}

// [[Rcpp::export]]
arma::mat varBetaDiag(arma::vec beta, arma::mat cov_mat, List cov_list, arma::vec ind, int nr) {
  int p = beta.size();
  int k = ind.size();
  int i_ind = 0;
  int j_ind = 0;
  ind = ind - 1; // indices are from R
  arma::mat sol          = arma::mat(p, p, fill::zeros);
  arma::mat cor_mat      = covToCor(cov_mat);
  arma::sp_mat inv_d_sig = diagSqrtSparse(cov_mat);
  arma::mat omega        = arma::inv_sympd(cor_mat);
  //omega = omega.clean(10^-8); // thresholded omega
  for (int i = 0; i < k; i++) {
    i_ind = ind.at(i);
    for (int j = 0; j < k; j++) {
      j_ind = ind.at(j);
      sol  +=  beta.at(i_ind) * beta.at(j_ind) * findRijDiagonal(cor_mat, inv_d_sig, cov_list, i_ind, j_ind);
    }
  }
  return omega * sol * omega;
}

/*** R
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

library(MASS)
library(tictoc)
### Parameters
n.est   <- 10000
n.org   <- 10000
n.ref   <- 2700
p       <- 5
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

### Helping matrices
k       <- createK(p)
Md      <- createMd(p)
Ms      <- createMs(p)
### Expression
A <- (diag(p^2) - Ms %*% kronecker(diag(p), cor.x.ref) %*% Md)
first.term   <- (diag(p^2) - Ms %*% kronecker(diag(p), cor.x.ref) %*% Md)
cov.mat.list <- outerProdRow(x.ref)
inv_sig      <- Matrix::sparseMatrix(x = 1 / sqrt(diag(cov.x.ref)), i = 1:p, j = 1:p, dims = c(p, p))


var.beta <-  var(t(
  replicate(1000, as.vector(omega.cor %*% cor(MASS::mvrnorm(n.ref, rep(0, p), cov.mat)) %*% beta.vec))))

var.cov.vec.mat <- var(t(
  replicate(1000, as.vector(cov(mvrnorm(n.ref, rep(0, p), cov.mat))), simplify = TRUE)))

var.cor.vec.mat <- var(t(replicate(1000,
                                   as.vector(cor(mvrnorm(n.ref, rep(0, p), cov.mat))), simplify = TRUE)))

cov.var <- findCovVar(cov.mat.list, cov.x.ref)
var.cor.formula  <- A %*% kronecker(inv_sig, inv_sig) %*% cov.var %*% kronecker(inv_sig, inv_sig) %*% t(A)


### Checking
### Variacne of covariance matrix
# 1
findCovTwoInd(cov.mat.list, 1, 1)
cov.var[(p + 1):(2 * p), (p + 1):(2 * p)]
var.cov.vec.mat[(p + 1):(2 * p), (p + 1):(2 * p)] * n.ref
(kronecker(cov.mat, cov.mat) %*% (diag(p^2) + k))[(p + 1):(2 * p), (p + 1):(2 * p)]

# 2
findCovTwoInd(cov.mat.list, 2, 1) /
cov.var[(2 * p + 1):(3 * p), (p + 1):(2 * p)]
# 3
findCovTwoInd(cov.mat.list, 1, 2) /
cov.var[(p + 1):(2 * p), (2 * p + 1):(3 * p)]



#### A
findPij(est_cor_mat = cor(x.ref), 1, 1)
A[(p + 1):(2 * p), (p + 1):(2 * p)]
# 2
findPij(est_cor_mat = cor(x.ref), 2, 1)
A[(2 * p + 1):(3 * p), (p + 1):(2 * p)]
# 3
findPij(est_cor_mat = cor(x.ref), 1, 2)
A[(p + 1):(2 * p), (2 * p + 1):(3 * p)]

#### A transpose
# 1
t(as.matrix(findPij(est_cor_mat = cor(x.ref), 1, 1)))
t(A)[(p + 1):(2 * p), (p + 1):(2 * p)]
# 2
t(as.matrix(findPij(est_cor_mat = cor(x.ref), 1, 2)))
t(A)[(2 * p + 1):(3 * p), (p + 1):(2 * p)]
# 3
t(as.matrix(findPij(est_cor_mat = cor(x.ref), 2, 1)))
t(A)[(p + 1):(2 * p), (2 * p + 1):(3 * p)]

### Plambda
PLa <- A %*% kronecker(diag(1 / sqrt(diag(cov.x.ref))), diag(1 / sqrt(diag(cov.x.ref))))
### Checking
# 1
findPLambdaij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 1)
PLa[(p + 1):(2 * p), (p + 1):(2 * p)]
# 2
findPLambdaij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 2, 1)
PLa[(2 * p + 1):(3 * p), (p + 1):(2 * p)]
# 3
findPLambdaij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 2)
PLa[(p + 1):(2 * p), (2 * p + 1):(3 * p)]

### LambdaP
LaP <- kronecker(diag(1 / sqrt(diag(cov.x.ref))), diag(1 / sqrt(diag(cov.x.ref)))) %*% t(A)
### Checking
# 1
findLambdaPij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 1)
LaP[(p + 1):(2 * p), (p + 1):(2 * p)]

# 2
findLambdaPij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 2, 1)
LaP[(2 * p + 1):(3 * p), (p + 1):(2 * p)]

# 3
findLambdaPij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, 1, 2)
LaP[(p + 1):(2 * p), (2 * p + 1):(3 * p)]



### Variance of correlation matrix
### Checking
(findRij(est_cor_mat = cov2cor(cov.mat), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 0, j = 0) / n.ref)
var.cor.vec.mat[1:p, 1:p]
var.cor.formula[1:p, 1:p] / n.ref

(findRij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 1, j = 2) / n.ref)
var.cor.vec.mat[(p + 1):(2 * p), (2 * p + 1):(3 * p)]
var.cor.formula[(p + 1):(2 * p), (2 * p + 1):(3 * p)] / n.ref

#### changing covariance matrix to correlation matrix
covToCor(cov.mat)
cov2cor(cov.mat)

### extracting sparse diagonal matrix from covariance matrix
diagSqrtSparse(cov.mat)
diag(cov.mat)^(-0.5)

### Checking final version
beta.vec <- c(1, 1, rep(0, p - 2))
cor.mat  <- cov2cor(cov.mat)
omega.cor <- solve(cor.mat)

kronecker(t(beta.vec), omega.cor) %*% var.cor.formula %*% kronecker(beta.vec, omega.cor) / n.ref
var.beta
(varBeta(beta.vec, cov.mat, cov.mat.list, c(1, 2), 3)) / n.ref






### Speed considerations
# microbenchmark::microbenchmark(
#   findRij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 1, j = 1) / n.ref,
#   findRijV2(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 1, j = 1) / n.ref
#
# )


#### Checking Gaussian Rij assumption
cov.mat <- 0.7^outer(1:p, 1:p, function(x, y) abs(y - x))
full.var <- kronecker(cov.mat, cov.mat) + kronecker(cov.mat, cov.mat) %*% createK(ncol(cov.mat))

# 1
findVarGaussCovInd(cov_mat = cov.mat, 1, 1)
full.var[(p + 1):(2 * p), (p + 1):(2 * p)]

# 2
findVarGaussCovInd(cov_mat = cov.mat, 2, 1)
full.var[(2 * p + 1):(3 * p), (p + 1):(2 * p)]

# 3
findVarGaussCovInd(cov_mat = cov.mat, 1, 2)
full.var[(p + 1):(2 * p), (2 * p + 1):(3 * p)]




findRijGauss

##### Comparing
(findRij(est_cor_mat = cor(x.ref), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 1, j = 2) / n.ref)

microbenchmark::microbenchmark(findCovVar(cov.mat.list, cov(x.ref))[1:5, 1:5])

cov_mat = cov(x.ref)
microbenchmark::microbenchmark(ECCCM:::findCovTwoInd(cov.mat.list, 0, 0), ECCCM:::findVarGaussCovInd(cov_mat, 0, 0))

microbenchmark::microbenchmark((findRij(est_cor_mat = cov2cor(cov.mat), inv_d_sig = inv_sig, cov_list = cov.mat.list, i = 0, j = 0) / n.ref),
findRijGauss(cov_mat, inv_d_sig = inv_sig, est_cor_mat = cov2cor(cov_mat), i = 0, j = 0))

*/

