// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// outerProdRow
List outerProdRow(arma::mat X);
RcppExport SEXP _ECCCM_outerProdRow(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(outerProdRow(X));
    return rcpp_result_gen;
END_RCPP
}
// covRows
arma::mat covRows(List cov_mats, int m, int k);
RcppExport SEXP _ECCCM_covRows(SEXP cov_matsSEXP, SEXP mSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cov_mats(cov_matsSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(covRows(cov_mats, m, k));
    return rcpp_result_gen;
END_RCPP
}
// findCovVar
arma::mat findCovVar(List list_cov_mat, arma::mat cov_mat);
RcppExport SEXP _ECCCM_findCovVar(SEXP list_cov_matSEXP, SEXP cov_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type list_cov_mat(list_cov_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_mat(cov_matSEXP);
    rcpp_result_gen = Rcpp::wrap(findCovVar(list_cov_mat, cov_mat));
    return rcpp_result_gen;
END_RCPP
}
// findCovBayes_old
arma::mat findCovBayes_old(List list_cov_mat, arma::mat cov_mat, arma::mat delta_gamma, bool only_diag);
RcppExport SEXP _ECCCM_findCovBayes_old(SEXP list_cov_matSEXP, SEXP cov_matSEXP, SEXP delta_gammaSEXP, SEXP only_diagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type list_cov_mat(list_cov_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_mat(cov_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_gamma(delta_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type only_diag(only_diagSEXP);
    rcpp_result_gen = Rcpp::wrap(findCovBayes_old(list_cov_mat, cov_mat, delta_gamma, only_diag));
    return rcpp_result_gen;
END_RCPP
}
// findCovBayes
arma::mat findCovBayes(List list_cov_mat, arma::mat cov_mat, arma::mat delta_gamma);
RcppExport SEXP _ECCCM_findCovBayes(SEXP list_cov_matSEXP, SEXP cov_matSEXP, SEXP delta_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type list_cov_mat(list_cov_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_mat(cov_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_gamma(delta_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(findCovBayes(list_cov_mat, cov_mat, delta_gamma));
    return rcpp_result_gen;
END_RCPP
}
// varFirstTerm
arma::mat varFirstTerm(arma::vec beta, arma::mat omega, List list_cov_mat, arma::vec ind);
RcppExport SEXP _ECCCM_varFirstTerm(SEXP betaSEXP, SEXP omegaSEXP, SEXP list_cov_matSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< List >::type list_cov_mat(list_cov_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(varFirstTerm(beta, omega, list_cov_mat, ind));
    return rcpp_result_gen;
END_RCPP
}
// varSecondTerm
arma::mat varSecondTerm(arma::vec omega_beta, List list_cov_mat, arma::vec ind);
RcppExport SEXP _ECCCM_varSecondTerm(SEXP omega_betaSEXP, SEXP list_cov_matSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type omega_beta(omega_betaSEXP);
    Rcpp::traits::input_parameter< List >::type list_cov_mat(list_cov_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(varSecondTerm(omega_beta, list_cov_mat, ind));
    return rcpp_result_gen;
END_RCPP
}
// findCovByInd
arma::mat findCovByInd(arma::mat x, int ind);
RcppExport SEXP _ECCCM_findCovByInd(SEXP xSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(findCovByInd(x, ind));
    return rcpp_result_gen;
END_RCPP
}
// findCovTwoInd
arma::mat findCovTwoInd(arma::mat x, int ind_a, int ind_b);
RcppExport SEXP _ECCCM_findCovTwoInd(SEXP xSEXP, SEXP ind_aSEXP, SEXP ind_bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ind_a(ind_aSEXP);
    Rcpp::traits::input_parameter< int >::type ind_b(ind_bSEXP);
    rcpp_result_gen = Rcpp::wrap(findCovTwoInd(x, ind_a, ind_b));
    return rcpp_result_gen;
END_RCPP
}
// varFirstTermNew
arma::mat varFirstTermNew(arma::vec beta, arma::mat omega, arma::mat x, arma::vec ind);
RcppExport SEXP _ECCCM_varFirstTermNew(SEXP betaSEXP, SEXP omegaSEXP, SEXP xSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(varFirstTermNew(beta, omega, x, ind));
    return rcpp_result_gen;
END_RCPP
}
// varSecondTermNew
arma::mat varSecondTermNew(arma::vec omega_beta, arma::mat x, arma::vec ind);
RcppExport SEXP _ECCCM_varSecondTermNew(SEXP omega_betaSEXP, SEXP xSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type omega_beta(omega_betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(varSecondTermNew(omega_beta, x, ind));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ECCCM_outerProdRow", (DL_FUNC) &_ECCCM_outerProdRow, 1},
    {"_ECCCM_covRows", (DL_FUNC) &_ECCCM_covRows, 3},
    {"_ECCCM_findCovVar", (DL_FUNC) &_ECCCM_findCovVar, 2},
    {"_ECCCM_findCovBayes_old", (DL_FUNC) &_ECCCM_findCovBayes_old, 4},
    {"_ECCCM_findCovBayes", (DL_FUNC) &_ECCCM_findCovBayes, 3},
    {"_ECCCM_varFirstTerm", (DL_FUNC) &_ECCCM_varFirstTerm, 4},
    {"_ECCCM_varSecondTerm", (DL_FUNC) &_ECCCM_varSecondTerm, 3},
    {"_ECCCM_findCovByInd", (DL_FUNC) &_ECCCM_findCovByInd, 2},
    {"_ECCCM_findCovTwoInd", (DL_FUNC) &_ECCCM_findCovTwoInd, 3},
    {"_ECCCM_varFirstTermNew", (DL_FUNC) &_ECCCM_varFirstTermNew, 4},
    {"_ECCCM_varSecondTermNew", (DL_FUNC) &_ECCCM_varSecondTermNew, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ECCCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
