\name{lines.locfit}
\alias{lines.locfit}
\title{
Add `locfit' line to existing plot
}
\usage{
lines(x, m, \dots)
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
\item{...}{
Other arguments to the default \code{\link{lines}} function.
}}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}, \code{\link{lines}}
}
\keyword{locfit}
