\name{plot.locfit.2d}
\title{
Plot a two-dimensional "preplot.locfit" object.
}
\usage{
plot.locfit.2d(x, type, xlab, \dots)
}
\arguments{
\item{x}{
Two dimensional \code{"preplot.locfit"} object.
}
\item{type}{
Type of plot. Can be \code{"contour"} (the default); 
\code{"persp"} or \code{"image"}.
}
\item{xlab}{
Labels for the predictor variable axis; should be a character
vector of length 2.
}
\item{...}{
Additional arguments to the \code{\link{contour}},
\code{\link{persp}} or \code{\link{image}} functions.
}
}

\description{
This function is not usually called directly. It will be called automatically
when plotting one-dimensional \code{"locfit"} or \code{"preplot.locfit"}
objects.
}

\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}, \code{\link{preplot.locfit}}
}
\keyword{locfit}
