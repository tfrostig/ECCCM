## usethis namespace: start
#' @useDynLib ECCCM, .registration = TRUE
## usethis namespace: end
NULL


## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL


usethis::use_package('CovTools')
usethis::use_package('Rcpp')
usethis::use_package('RcppArmadillo')
usethis::use_package('MASS', 'Suggests')
