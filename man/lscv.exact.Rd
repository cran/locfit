\name{lscv.exact}
\alias{lscv.exact}
\title{
Exact LSCV Calculation
}
\usage{
lscv.exact(x, alpha, ...)
}
\description{
  This function performs the exact computation of the least squares
  cross validation statistic for one-dimensional kernel density estimation
  and a constant bandwidth.

  At the time of writing, it is implemented only for the Gaussian
  kernel (with the standard deviation of 0.4; Locfit's standard).
}
\arguments{
  \item{x}{Numeric data vector}
  \item{alpha}{The bandwidth. Note that this function only works with
    a constant bandwidth, and in accordance with Locfit's standard, this
    can be specified as \code{alpha=c(0,h)}. It can also be given as
    a single number, which will be interpreted as a constant
    (not NN) bandwidth. }
  \item{...}{not used.}
}
\value{
  A vector of the LSCV statistic and the fitted degrees of freedom.
}
\examples{
data(geyser, package="locfit")
lscv.exact(geyser,alpha=0.25)
# equivalent form using lscv
lscv(geyser,alpha=0.25,exact=TRUE)
}
\seealso{
  \code{\link{lscv}},
  \code{\link{lscvplot}}
}

\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
