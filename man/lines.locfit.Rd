\name{lines.locfit}
\alias{lines.locfit}
\title{
Add `locfit' line to existing plot
}
\usage{
\method{lines}{locfit}(x, m, tr=x$trans, ...)
}
\description{
Adds a Locfit line to an existing plot.
}
\arguments{
\item{x}{
\code{"locfit"} object. Should be a model with one predictor.
}
\item{m}{
Number of points to evaluate the line at. Default: 100.
}
\item{tr}{
Transformation function to use for plotting. Default is the inverse
link function, or the identity function if derivatives are requested.
}
\item{...}{
Other arguments to the default \code{lines} function.
}}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}
}
\keyword{smooth}
