\name{plot.gcvplot}
\alias{plot.gcvplot}
\title{
Produce a cross-validation plot.
}
\usage{
plot.gcvplot(x,...)
}
\description{
Plots the value of the GCV (or other statistic) in a \code{gcvplot} object
against the degrees of freedom of the fit.
}
\arguments{
  \item{x}{ A "gcvplot" object, produced by \code{\link{gcvplot}},
     \code{\link{aicplot}} etc.}
  \item{...}{ Other arguments to \code{\link{plot}} .}
}

\examples{
data(ethanol)
plot(gcvplot(NOx~E,data=ethanol,alpha=seq(0.2,1.0,by=0.05)))
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gcv}},
  \code{\link{aicplot}},
  \code{\link{cpplot}},
  \code{\link{gcvplot}},
  \code{\link{lcvplot}}
}

\keyword{methods}
