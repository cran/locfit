\name{plot.scb}
\alias{plot.scb}
\title{
Plot method for simultaneous confidence bands
}
\usage{
plot.scb(x, add=F, ...)
}
\description{
  Plot method for simultaneous confidence bands created by the
  \code{\link{scb}} function.
}
\arguments{
\item{x}{ \code{"scb"} object created by \code{\link{scb}()}.  }
\item{add}{If true, bands will be added to the existing plot.}
}
\examples{
# corrected confidence bands for a linear logistic model
data(insect)
fit <- scb(deaths ~ lconc, type=4, w=nins, data=insect,
   deg=1, family="binomial", kern="parm")
plot(fit)
}
\seealso{
\code{\link{scb}}
}
\keyword{methods}
