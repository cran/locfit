\name{lfbas}
\alias{lfbas}
\title{
User-specified basis functions for Locfit.
}
\usage{
locfit(...,basis)
}
\description{
  By default, Locfit uses a polynomial basis as the fitting functions.
  An alternative set of basis functions can be specified through
  the \code{basis} argument to \code{\link{locfit.raw}}.
  \code{lfbas} is a wrapper function used internally in Locfit's processing.

  To use the \code{basis} argument, you must write a function
  \code{f(x,t)} to evaluate the basis function at a fitting
  point \code{t} and data point(s) \code{x}. See below for an example.
  Note that the basis function will be called with multiple
  data points simultaneously, so should assume \code{x} is a matrix
  with \code{d} columns, where \code{d} is the number of independent
  variables. The fitting point \code{t} will always be a single point,
  in a vector of length \code{d}.

  The basis function should return a matrix, with each column being the
  evaluation of one fitting function.

  To ensure that the returned fit estimates the mean function, the
  first component of the basis should be 1; the remaining components
  should be 0 when \code{x=t}. To help ensure Locfit's interpolation
  routines are meaningful, the next \code{d} components should represent
  partial derivatives at \code{x=t}. Specifically, df(x,t)/dx[i],
  evaluated at \code{x=t}, should be the unit vector with 1 in the
  (i+1)st position.

  Violation of these rules can be useful, if functionals
  other than the mean are of interest. But this will require extreme
  care.

  Specifying a user basis results in extensive use of the
  \code{call_S} function, which may be slow. Speed may also
  be significantly affected by different ways of writing
  the basis.
}
\examples{
# Specify a bivariate linear with interaction basis.
data(ethanol)
my.basis <- function(x,t)
{
  u1 <- x[, 1] - t[1]
  u2 <- x[, 2] - t[2]
  cbind(1, u1, u2, u1 * u2)
}
fit <- locfit(NOx~E+C, data=ethanol, scale=0, basis=my.basis)
# With this basis, Locfit's standard interpolation and plot methods
# should be reasonable.
plot(fit,get.data=T)

# Estimation of change points. This provides an alternative to using
# left() and right(), and can easily be modified to detecting
# a change in slopes or other parameters. Note that the first
# component is the indicator of x>t, so the coefficient estimates
# the size of the change, assuming the change occurs at t.
data(penny)
my.basis <- function(x,t) cbind(x>t,1,x-t)
xev <- (1945:1988) + 0.5
fit <- locfit(thickness~year, data=penny, alpha=c(0,10), ev=xev, basis=my.basis)
# The plot will show peaks where change points are most likely.
# in S4, S-Plus 5 etc,
# plot(preplot(fit,where="fitp")^2, type="b") is an alternative.
plot(xev, predict(fit,where="fitp")^2, type="b")

# Estimate the mean function using local linear regression, with
# discontinuities at 1958.5 and 1974.5.
# The basis functions must consist of the constant 1, the linear term
# x-t, and indicator functions for two of the three sections.
# Note the care taken to ensure my.basis(t,t) = c(1,0,0,0) for all t.
my.basis <- function(x,t)
{ ret <- NULL
  if (t<1958.5) ret <- cbind(1, x>=1958.5, x>1974.5, x-t)
  if (t>1974.5) ret <- cbind(1, x<=1974.5, x<1958.5, x-t)
  if (is.null(ret))
    ret <- cbind(1, x<1958.5, x>1974.5, x-t)
  ret
}
fit <- locfit(thickness~year, data=penny, alpha=c(0,10), ev=xev, basis=my.basis)
plot(preplot(fit,where="fitp", get.data=T))
}
\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}}
}

\keyword{locfit}
