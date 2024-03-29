% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Frequentist_functions.R
\name{analyzeRefGauss}
\alias{analyzeRefGauss}
\title{analyzeRefGauss - Conduct full analysis based on reference panel under the assumption of Gaussian reference}
\usage{
analyzeRefGauss(
  marg.beta.hat,
  ld.mat,
  n.o,
  n.r,
  sigma.method = "conservative",
  method.filter = "none",
  method.test = "BH",
  qu = 1
)
}
\arguments{
\item{marg.beta.hat}{Marginal regression coefficients, inner product of G_o and y divided by the number of observations in original study}

\item{ld.mat}{Covariance matrix based on the reference panel}

\item{n.o}{The number of observations in original study}

\item{n.r}{The number of observation in reference panel}

\item{sigma.method}{String indicating how to estimate sigma, currently supports three methods, 'conservative', sigma is estimated as 1.
'estimate' the usual estimator of the variance (assuming y is scaled sigma = 1 - var(x.r beta_mc)), and 'semi.conservative' which is similar to
the 'estimate' method, but the beta_mc is thresholded.}

\item{method.filter}{multiple adjustment method to use in thresholding of the coefficient (see p.adjust)}

\item{method.test}{multiple adjustment method  to use when testing the adjusted coefficients (see p.adjust)}

\item{qu}{Coefficients with p-value lower than 'qu' are thresholded to 0.}
}
\value{
list containing `test.correct` data.frame with the adjusted coefficient and adjusted variance, with two sided testing p-value.
`test.naive` data.frame with the adjusted coefficient but not adjusted variance, with two sided testing p-value.
`add.var` the additional variance, `sigma` the estimated sigma.
`var.omega.beta` a vector with the number of coefficient used and the proportion of weight used.
}
\description{
analyzeRefGauss - Conduct full analysis based on reference panel under the assumption of Gaussian reference
}
