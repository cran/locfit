\name{summary.gcvplot}
\alias{summary.gcvplot}
\title{
Summary method for a gcvplot structure.
}
\usage{
summary(object)
}
\description{
  Computes a short summary for a generalized cross-validation plot structure
}
\arguments{
  \item{object}{A "gcvplot" structure produced by a call to
    \code{\link{gcvplot}}, \code{\link{cpplot}} e.t.c.
}
}

\value{
  A matrix with two columns; one row for each fit computed in the
  \code{\link{gcvplot}} call.
  The first column is the fitted degrees
  of freedom; the second is the GCV or other criterion computed.
}

\examples{
data(ethanol)
summary(gcvplot(NOx~E,data=ethanol,alpha=seq(0.2,1.0,by=0.05)))
}

\seealso{
  \code{\link{locfit}},
  \code{\link{gcv}},
  \code{\link{gcvplot}}
}

\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
