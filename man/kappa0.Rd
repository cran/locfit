\name{kappa0}
\alias{kappa0}
\title{
Critical Values for Simultaneous Confidence Bands.
}
\usage{
kappa0(formula, cov=0.95, \dots)
}
\arguments{
\item{formula}{
Local regression model formula.
}
\item{cov}{
Coverage Probability for critical values.
}
\item{ldots}{
Other arguments to \code{\link{locfit}}.
}
}
\description{
  The geometric constants for simultaneous confidence bands are computed,
  as described in Sun and Loader (1994) (bias adjustment is not implemented
  here). These are then passed to the \code{\link{crit}} function, which
  computes the critical value for the confidence bands.

  The method requires both the weight diagrams l(x), the
  derivative l'(x) and second derivatives l''(x). These are
  implemented exactly for a constant bandwidth; that
  is, \code{alpha=c(0,h)} for some \code{h}. For nearest
  neighbor bandwidths, the computations are approximate.

  The theoretical justification for the bands are computed using
  the spherical symmetry of the Normal distributions. For non-normal
  distributions, and likelihood models, one relies on central
  limit and related theorems...

  Computation uses the product Simpson's rule to evaluate the
  multidimensional integrals. Expect this to be slow in more
  than one dimension. The \code{mint} argument controls the
  precision.
}
\value{
  A list with components for the critical value, geometric constants,
  e.t.c. Can be passed directly to \code{\link{plot.locfit}} as the
  \code{crit} argument.
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}},
\code{\link{crit}}, \code{\link{crit<-}}.
}
\examples{
# compute and plot simultaneous confidence bands
data(ethanol)
fit <- locfit(NOx~E,data=ethanol)
crit(fit) <- kappa0(NOx~E,data=ethanol)
plot(fit,crit=crit,band="local")
}
\references{
Sun, J. and Loader, C. (1994). Simultaneous confidence bands for linear
  regression and smoothing. Annals of Statistics 22, 1328-1345.
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
