\name{left}
\alias{left}
\title{
One-sided left smooth for a Locfit model.
}
\usage{
left(x)
%locfit(y~left(x), ev="data", deg=1)
}
\description{
  The \code{left()} function is used in a locfit model formula
  to specify a one-sided smooth: when fitting at a point \eqn{x},
  only data points with \eqn{x_i \le x} should be used.
  This can be useful in estimating points of discontinuity,
  and in cross-validation for forecasting a time series.

  When using this function, it will usually be necessary to specify an
  evaluation structure, since the fit is not smooth and locfit's
  interpolation methods are unreliable. Also, it is usually best
  to use \code{deg=0} or \code{deg=1}, otherwise the fits may be too
  variable. If nearest neighbor bandwidth specification is used,
  it does not recognize \code{left()}.

}

\arguments{
  \item{x}{numeric variable.}
}

\examples{
# compute left and right smooths
data(penny)
xev <- (1945:1988)+0.5
fitl <- locfit(thickness~left(year), alpha=c(0,10), deg=1, ev=xev, data=penny)
fitr <- locfit(thickness~right(year),alpha=c(0,10), deg=1, ev=xev, data=penny)
# plot the squared difference, to show the change points.
plot( xev, (predict(fitr,where="ev") - predict(fitl,where="ev"))^2 )
}

\seealso{
  \code{\link{lfbas}},
  \code{\link{locfit}},
  \code{\link{right}}
}

\keyword{models}
