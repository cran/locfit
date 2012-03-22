\name{panel.locfit}
\alias{panel.locfit}
\title{
Locfit panel function
}
\usage{
panel.locfit(x, y, subscripts, z, xyz.labs, xyz.axes, xyz.mid, xyz.minmax, 
  xyz.range, col.regions, at, drape, contour, region, groups, ...)
}
\description{
This panel function can be used to add locfit fits to plots generated
by trellis.

Currently it works with \code{xyplot} for 1-d fits
and crudely with \code{wireframe} for 2-d fits.
}
\arguments{
  \item{...}{Most Locfit arguments can be passed through \code{xyplot}.}
}
\examples{
stopifnot(require(lattice))

# a simple multi-panel display
data(ethanol, package="locfit")
xyplot(NOx ~ E | C, data=ethanol, panel=panel.locfit)

# The second example uses some Locfit optional arguments.
# Note we can pass the alpha (bandwidth) and family arguments directly to
# xyplot. The cens argument must be given in full; not as a data frame variable.
# The resulting plot does not (yet) distinguish the censored points, but
# the fit will correctly apply censoring.
data(border, package="locfit")
xyplot(runs ~ day, data=border, panel=panel.locfit, family="poisson",
 alpha=0.3, cens=border$no)
}
\seealso{
  \code{\link{locfit}}, \code{\link{plot.locfit.3d}}, \code{xyplot}.
}

\keyword{internal}
