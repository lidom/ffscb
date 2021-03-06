% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confidence_band.R
\name{confidence_band_fragm}
\alias{confidence_band_fragm}
\title{Makes confidence bands for fragmentary functional data}
\usage{
confidence_band_fragm(x, diag.cov.x, tau = NULL, t0 = NULL,
  df = NULL, type = c("FFSCB.z", "FFSCB.t", "KR.z", "KR.t", "naive.t"),
  conf.level = 0.95, n_int = 10, tol = NULL)
}
\arguments{
\item{x}{Functional parameter estimate (for instance, the empirical mean function). It can be either a vector or \link{fd} object from \link{fda}.}

\item{diag.cov.x}{diag(Cov(x)), in which x is the functional estimator (for instance, the covariance function of the empirical mean function). It can be either matrix or \link{bifd} object from \link{fda}. The eigen decomposition of Cov(X) can be used instead.}

\item{tau}{Pointwise standard deviation of the standardized and differentiated sample functions. Can be estimated by tau_fun().}

\item{t0}{Parameter t0 of the fast and fair simultaneous confidence bands.}

\item{df}{Degrees of freedom parameter for the t-distribution based bands 'FFSCB.t', 'KR.t', and 'naive.t'. If x is the empirical mean function, set df=n-1, where n denotes the sample size.}

\item{type}{The band(s) to be constructed.
\itemize{
  \item FFSCB.z : Fast'n'Fair (adaptive) simultaneous confidence band based for a Gaussian functional parameter estimate.
  \item FFSCB.t : Fast'n'Fair (adaptive) simultaneous confidence band based for a t-distributed functional parameter estimate.
  \item KR.z : The constant simultaneous confidence band based on the classical Kac-Rice (KR) formula for Gaussian random functions.
  \item KR.t : The constant simultaneous confidence band based on the classical Kac-Rice (KR) formula for t-distributed random functions.
}}

\item{conf.level}{A vector of confidence levels for the bands to achieve.}

\item{n_int}{Number of intervals for the piecewise linear confidence bounds.}

\item{tol}{tolerance 'tol' parameter used by stats::uniroot(). If tol=NULL, we use tol=.Machine$double.eps^0.32 which increases the default accuracy used by uniroot().}
}
\value{
confidence_band_fragm
}
\description{
Makes confidence bands for fragmentary functional data
}
\references{
\itemize{
   \item Liebl, D. and Reimherr, M. (2019). Fast and fair simultaneous confidence bands.
}
}
