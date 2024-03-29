% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Frequentist_functions.R
\name{analyzeRef}
\alias{analyzeRef}
\title{analyzeRef - Conduct full analysis based on reference panel}
\usage{
analyzeRef(
  marg.beta.hat,
  x.r,
  n.o,
  sigma.method = "conservative",
  method.filter = "none",
  method.test = "BH",
  qu = 0.05,
  to.diag = FALSE,
  is.marg = TRUE
)
}
\arguments{
\item{marg.beta.hat}{Marginal regression coeffcients, inner product of G_o and y divided by the number of observations in original study}

\item{x.r}{Reference panel must be centered}

\item{n.o}{The number of observations in original study}

\item{sigma.method}{String indicating how to estimate sigma, currently supports three methods, 'conservative', sigma is estimated as 1.
'estimate' the usual estimator of the variance (assuming y is scaled sigma = 1 - var(x.r beta_mc)), and 'semi.conservative' which is similar to
the 'estimate' method, but the beta_mc is thresholded.}

\item{method.filter}{multiple adjustment method to use in thresholding of the coefficient (see p.adjust)}

\item{method.test}{multiple adjustment method  to use when testing the adjusted coefficients (see p.adjust)}

\item{qu}{P-value threshold.}

\item{threshold.cov.seq}{user-defined threshold value. If it is a vector of regularization values, it automatically selects one that minimizes cross validation risk.}
}
\value{
list containing `test.correct` data.frame with the adjusted coefficeint and adjusted variance, with two sided testing p-value.
`test.naive` data.frame with the adjusted coefficeint but not adjusted variance, with two sided testing p-value.
 `add.var` the additional variance, `sigma` the estimated sigma.
 `var.omega.beta` a vector with the number of coefficient used and the proportion of weight used.
}
\description{
analyzeRef - Conduct full analysis based on reference panel
}
