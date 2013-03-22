\name{plot.locfit}
\title{
Plot an object of class locfit.
}
\usage{
plot.locfit(x, xlim, pv, tv, mpv, mtv, deriv, se.fit, tr,
        what, get.data, f3d, \dots)
}
\arguments{
\item{x}{
locfit object.
}
\item{xlim}{
Plotting limits. Eg. \code{xlim=c(0,0,1,1)} plots over the unit square in
two dimensions. Default is bounding box of the data.
}
\item{pv}{
Panel variables, to be varied within each panel of a plot.  May be
        specified as a character vector, or variable numbers.  There must be
        one or two panel variables; default is all variables in one or two
        dimensions; Variable 1 in three or more dimensions.
May by specified using either variable numbers or names.
}
\item{tv}{
Trellis variables, to be varied from panel to panel of the plot.
}
\item{mpv}{
Controls the plot resolution (within panels, for trellis displays).
Default is 100 points in one dimension; 40 points (per dimension)
in two or more dimensions.
}
\item{mtv}{
Number of points for trellis variables; default 6.
}
\item{deriv}{
derivative argument, to plot local slopes and curvatures.
Default \code{numeric(0)}. See \code{predict.locfit} documentation for more details.
}
\item{se.fit}{
Should variability bands be plotted? Default FALSE. This is
        presently available for one predictor and one dimensional
        cross sections.
}
\item{tr}{
Transformation function to use for plotting. Default is the inverse
link function, or the identity function if derivatives are requested.
}
\item{what}{
What to plot. See \code{predict.locfit}.
}
\item{get.data}{
If \code{TRUE}, original data is added to the plot. Default: \code{FALSE}.
}
\item{f3d}{
Force the \code{locfit.3d} class on the prediction object, thereby generating
a trellis style plot. Default: \code{FALSE}, unless a \code{tv} argument is provided.
}
\item{...}{
Other arguments to \code{plot.locfit.1d}, \code{plot.locfit.2d} or
\code{plot.locfit.3d} as appropriate.
}}
\description{
The \code{plot.locfit} function generates grids of ploting points, followed
by a call to \code{predict.locfit}. The returned object is then passed to
\code{plot.locfit.1d}, \code{plot.locfit.2d} or \code{plot.locfit.3d} as appropriate.
}
\examples{
x <- rnorm(100)
y <- dnorm(x)+rnorm(100)/5
plot(locfit(y~x),vband=T)
x <- cbind(rnorm(100),rnorm(100))
plot(locfit(~x),type="persp")
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit.1d}}, \code{\link{plot.locfit.2d}},
\code{\link{plot.locfit.3d}}, \code{\link{lines.locfit}}, \code{\link{predict.locfit}}
}
\keyword{locfit}


