\name{plot.gcvplot}
\alias{plot.gcvplot}
\title{
Produce a cross-validation plot.
}
\usage{
plot.gcvplot(x,...)
}
\arguments{
  \item{x}{ A "gcvplot" object, produced by \code{\link{gcvplot}},
     \code{\link{aicplot}} etc.
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

\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
