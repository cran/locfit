\name{gcv}
\alias{gcv}
\title{
Compute generalized cross-validation statistic.
}
\usage{
gcv(x, \dots)
}
\description{
  The calling sequence for \code{gcv} matches those for the
  \code{\link{locfit}} or \code{\link{locfit.raw}} functions.
  The fit is not returned; instead, the returned object contains
  Wahba's generalized cross-validation score for the fit.

  The GCV score is exact (up to numerical roundoff) if the
  \code{ev="data"} argument is provided. Otherwise, the residual
  sum-of-squares and degrees of freedom are computed using locfit's
  standard interpolation based approximations.

  For likelihood models, GCV is computed using -2 times the log-likelihood
  in place of the residual sum of squares; I know of no results validating
  this interpretation.
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gcvplot}}
}

\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
