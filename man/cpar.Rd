\name{cpar}
\alias{cpar}
\title{
Conditionally parametric term for a Locfit model.
}
\usage{
locfit(y~x1+cpar(x2))
}
\description{
  A term entered in a \code{\link{locfit}} model formula using
  \code{cpar} will result in a fit that is conditionally parametric.
  If the formula is \code{y~x1+cpar(x2)}, the resulting surface,
  for fixed \code{x1}, will be quadratic (using the default \code{deg=2})
  in \code{x2}. Note that the coefficients of the quadratic may depend
  on \code{x1}; this distinguishes a conditionally parametric model
  from a semiparametric model.

}

\arguments{
  \item{x2}{numeric variable.}
}

\examples{
data(ethanol)
# fit a conditionally parametric model
fit <- locfit(NOx~E+cpar(C),data=ethanol)
plot(fit)
# one way to force a parametric fit with locfit
fit <- locfit(NOx~cpar(E),data=ethanol)
}

\seealso{
  \code{\link{locfit}}
}

\keyword{models}
