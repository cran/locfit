\name{lscv}
\alias{lscv}
\title{
Least Squares Cross Validation Statistic.
}
\usage{
lscv(x, \dots,en=2)
}
\description{
  The calling sequence for \code{lscv} matches those for the
  \code{\link{locfit}} or \code{\link{locfit.raw}} functions.
  Note that this function is only designed for density estimation
  in one dimension. The returned object contains the
  least squares cross validation score for the fit.

  The computation of \(\int \hat f(x)^2 dx\) is performed numerically.
  For kernel density estimation, this is unlikely to agree exactly
  with other LSCV routines, which may perform the integrtion analytically.
}
\arguments{
  \item{x}{model formula}
  \item{...}{other arguments to locfit}
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{lscvplot}}
}

\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
