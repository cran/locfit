\name{panel.locfit}
\alias{panel.locfit}
\title{
Locfit panel function
}
\usage{
panel.locfit(x, y, subscripts, z, xyz.labs, xyz.axes, xyz.mid, xyz.minmax,
  xyz.range, col.regions, at, drape, contour, region, groups, ...)
%xyplot(formula, data, panel=panel.locfit, ...)
}
\description{
This panel function can be used to add \code{locfit} fits to plots generated
by trellis.

At the time of writing, it works with \code{xyplot} for 1-d fits
and crudely with \code{wireframe} for 2-d fits.
}
\arguments{
  \item{x, y}{x and y data.}
  \item{subscripts}{Subscripts of indices in the original data frame.}
  \item{z}{z data.}
  \item{xyz.labs, xyz.axes, xyz.mid, xyz.minmax, xyz.range, col.regions,
  at, drape, contour, region}{Arguments to be passed to 3D
  trellis panel functions.}
  \item{groups}{Currently not used.}
  \item{...}{Most Locfit arguments can be passed through \code{xyplot}.}
}
\examples{
library("lattice")

# a simple multi-panel display
data(ethanol)
xyplot(NOx~E|C,data=ethanol,panel=panel.locfit)

# The second example uses some Locfit optional arguments.
# Note we can pass the alpha (bandwidth) and family arguments directly to
# xyplot. The cens argument must be given in full; not as a data frame variable.
# The resulting plot does not (yet) distinguish the censored points, but
# the fit will correctly apply censoring.
data(border)
xyplot(runs~day,data=border,panel=panel.locfit,family="poisson",
 alpha=0.3,cens=border$no)
}
\seealso{
  \code{\link{locfit}}, \code{\link{plot.locfit.3d}}, \code{xyplot}.
}

\keyword{smooth}
