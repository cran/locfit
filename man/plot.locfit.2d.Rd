\name{plot.locfit.2d}
\alias{plot.locfit.2d}
\title{
Plot a two-dimensional "preplot.locfit" object.
}
\usage{
plot.locfit.2d(x, pv, tv, main, type = "contour", xlab, ylab,
               zlab = x$yname, ...)
}
\arguments{
\item{x}{
Two dimensional \code{"preplot.locfit"} object.
}
\item{pv}{Panel variables, to be varied within each panel of a plot.
          May be specified as a character vector, or variable numbers.
          There must be one or two panel variables; default is all
          variables in one or two dimensions; Variable 1 in three
          or more dimensions.  May by specified using either variable
          numbers or names.}
\item{tv}{Trellis variables, to be varied from panel to panel of the
          plot.}
\item{main}{main title for the plot.}
\item{type}{
Type of plot. Can be \code{"contour"} (the default);
\code{"persp"} or \code{"image"}.
}
\item{xlab}{
Labels for the predictor variable axis; should be a character
vector of length 2.
}
\item{ylab}{Label for the y-axis.}
\item{zlab}{If \code{type="persp"}, the label for the z-axis.}
\item{...}{
Additional arguments to the \code{contour},
\code{persp} or \code{image} functions.
}
}

\description{
This function is not usually called directly. It will be called automatically
when plotting one-dimensional \code{"locfit"} or \code{"preplot.locfit"}
objects.
}

\seealso{
  \code{\link{locfit}},
  \code{\link{plot.locfit}},
  \code{\link{preplot.locfit}}
}
\keyword{smooth}
