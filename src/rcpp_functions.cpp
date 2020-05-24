// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


/// Finding the covariance of two rows of covariance matrix. /arma::mat
//estCov(arma::mat g_r, arma::mat cov_est, int i, int j) {

// [[Rcpp::export]]
List outerProdRow(arma::mat X) {
  int n = X.n_rows;
  int i;
  List res_list = List(n);
  for (i = 0; i < n; i++) {
    res_list[i] = X.row(i).t() * X.row(i);
  }
  return res_list;
}



// find the covariance matrix between the rows of m and k of a covariance matrix
// [[Rcpp::export]]
arma::mat covRows(List cov_mats, int m, int k) {
  int nl = cov_mats.size();
  int i;
  arma::mat examp_mat = cov_mats[0];
  int p = examp_mat.n_cols;
  arma::mat sol = arma::mat(p, p, fill::zeros);
  arma::mat mean_cov_mat = arma::mat(p, p, fill::zeros);
  for (i = 0; i < nl; i++) {
    arma::mat temp_mat = cov_mats[i];
    sol += temp_mat.col(m) * temp_mat.col(k).t();
    // Rcout << "The value of sol : " << temp_mat.col(m) << "\n";
    mean_cov_mat += temp_mat;
  }
  mean_cov_mat = mean_cov_mat / nl;
  sol = ((sol / nl) - mean_cov_mat.col(m) * mean_cov_mat.col(k).t());
  return sol;
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
arma::mat findCovBayes(List list_cov_mat, arma::mat cov_mat, arma::mat delta_gamma) {
  int i = 0;
  int p = cov_mat.n_cols;
  arma::mat sol = arma::mat(p ,p, fill::zeros);
  for (i = 0; i < p; i++) {
      sol += delta_gamma.at(i, i) * 1 * covRows(list_cov_mat, i, i);
  }
  return sol;
}


// [[Rcpp::export]]
arma::mat findCovByInd(arma::mat x, int ind) {
  arma::vec x_ind = x.col(ind);
  int nrow = x_ind.size();
  int ncol = x.n_cols;
  arma::rowvec mean_cov = zeros<rowvec>(ncol);
  int j;
  for (j = 0; j < nrow; j++) {
    x.row(j) *= x_ind.at(j);
    mean_cov += x.row(j);
  }
  mean_cov = mean_cov / nrow;
  x.each_row() -= mean_cov;
  return x;
}

// [[Rcpp::export]]
arma::mat findCovTwoInd(arma::mat x, int ind_a, int ind_b) {
  int n = x.n_rows;
  return findCovByInd(x, ind_a).t() * findCovByInd(x, ind_b) / n;
}

// [[Rcpp::export]]
arma::mat varFirstTerm(arma::vec beta, arma::mat omega,  arma::mat x, arma::vec ind) {
  int i = 0;
  int j = 0;
  int p = beta.size(); // thresholded omega_beta
  int k = ind.size();
  int i_ind = 0;
  int j_ind = 0;
  ind = ind - 1; // indices are from R
  arma::mat sol = arma::mat(p, p, fill::zeros);
  for (i = 0; i < k; i++) {
    i_ind = ind.at(i);
    for (j = 0; j < k; j++) {
      j_ind = ind.at(j);
      sol  +=  beta.at(i_ind) * beta.at(j_ind) * omega * findCovTwoInd(x, i_ind, j_ind) * omega;
    }
  }
  return sol;
}

// Input: omega_beta: linear transformed Sigma^{-1} * beta
//        list_cov_mat: list of rows outerproduct (only)
// [[Rcpp::export]]
arma::mat varSecondTerm(arma::vec omega_beta, arma::mat x, arma::vec ind) {
  int i = 0;
  int j = 0;
  int p = omega_beta.size(); // thresholded omega_beta
  int k = ind.size();
  int i_ind = 0;
  int j_ind = 0;
  ind = ind - 1; // indices are from R
  arma::mat sol = arma::mat(p, p, fill::zeros);
  for (i = 0; i < k; i++) {
    i_ind = ind.at(i);
    for (j = 0; j < k; j++) {
      j_ind = ind.at(j);
      sol +=  omega_beta.at(i_ind) * omega_beta.at(j_ind) * findCovTwoInd(x, i_ind, j_ind);
    }
  }
  return sol;
}

// [[Rcpp::export]]
double quadForm(arma::vec x, arma::mat S) {
  double sol = as_scalar(x.t() * S * x);
  return sol;
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





/*** R
x.ref <- MASS::mvrnorm(100, rep(0, 3), diag(3))
x.ref <- scale(x.ref)
cov.x.ref <- cov(x.ref)
m <- 3
gamma.vec <- c(1, rep(0, 9))
cov.mat.list <- outerProdRow(x.ref)
delta.gamma       <- matrix(0, ncol = 10, nrow = 10)
diag(delta.gamma)[gamma.vec] <- 1
a.var <- findCovVar(cov.mat.list, cov.x.ref)
covRows(cov.mat.list, 1, 1)

x <- c(1,2)
S <- diag(2)

quadForm(x, S)
addVarGauss(S, S, x, 5, 5)

est.beta.hat <- x
cov.list.r <- list('cov' = S, 'omega' = S)
n.r <- 5
n.o <- 5

term.2.var.theoertical <- (drop(t(est.beta.hat) %*% cov.list.r$cov %*% est.beta.hat) * cov.list.r$omega + est.beta.hat %*% t(est.beta.hat)) / n.r
term.3.var.theoertical <- (drop(t(est.beta.hat) %*% cov.list.r$omega %*% est.beta.hat) * cov.list.r$cov + est.beta.hat %*% t(est.beta.hat)) / n.o

term.2.var.theoertical + term.3.var.theoertical


  */

