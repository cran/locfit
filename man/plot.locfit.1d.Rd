\name{plot.locfit.1d}
\alias{plot.locfit.1d}
\title{
Plot a one dimensional "preplot.locfit" object.
}
\usage{
plot.locfit.1d(x, pv, tv, main = "", xlab = "default", ylab = x$yname,
     type = "l", ylim, add = FALSE, lty = 1, col = 1, ...)
}
\arguments{
\item{x}{
One dimensional \code{"preplot.locfit"} object.
}
\item{pv}{Panel variables, to be varied within each panel of a plot.  May be
        specified as a character vector, or variable numbers.  There must be
        one or two panel variables; default is all variables in one or two
        dimensions; Variable 1 in three or more dimensions.
May by specified using either variable numbers or names.}
\item{tv}{Trellis variables, to be varied from panel to panel of the plot.}
\item{main}{Main title for the plot.}
\item{xlab}{label for the x-axis.}
\item{ylab}{label for the y-axis.}
\item{type}{type of plot (\code{"l"} for line, \code{"p"} for points,
            etc.)}
\item{ylim}{limit (min and max) for the y-axis.}
\item{add}{
If \code{TRUE}, the plot will be added to the existing plot.
}
\item{lty}{line type.}
\item{col}{color.}
\item{...}{
Additional arguments to the \code{plot} function.
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
