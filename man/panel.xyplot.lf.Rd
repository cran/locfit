\name{panel.xyplot.lf}
\alias{panel.xyplot.lf}
\title{
Locfit panel function
}
\usage{
panel.xyplot.lf(x, y, subscripts, clo, cup, wh, type, ...)
}
\description{
Panel function used by \code{\link{plot.locfit.3d}} for one dimensional
plots.
}
\arguments{
  \item{x}{x data}
  \item{y}{y data}
  \item{subscripts}{subscripts giving indices in the original data frame.}
  \item{clo}{lower limits of the confidence band.}
  \item{cup}{upper limits of the confidence band.}
  \item{wh}{\code{"f"} or \code{"d"}.}
  \item{type}{Plot type, to be passed to \code{panel.xyplot}.}
  \item{...}{Other arguments to be passed to \code{panel.xyplot}.}
}
  \seealso{
  \code{\link{plot.locfit.3d}}
}

\keyword{smooth}
