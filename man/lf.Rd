\name{lf}
\alias{lf}
\title{
Locfit term in Additive Model formula
}
\usage{
lf(..., alpha=0.7, deg=2, scale=1, kern="tcub", ev="tree", maxk=100)
}
\description{
  This function is used to specify a smooth term in a \code{gam()}
  model formula.

  At the time of writing, \code{gam} is not available in R.
}
\arguments{
  \item{...}{numeric predictor variable(s)}
  \item{alpha}{and other arguments are as in \code{\link{locfit.raw}()}}.
}
\examples{
# fit an additive semiparametric model to the ethanol data.
if (exists("gam")) {  # shut up R's checks
fit <- gam(NOx~lf(E)+C,data=ethanol)
plot(fit)
}
}
\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gam.lf}},
  \code{gam}
}

\keyword{models}
% Converted by Sd2Rd version 0.2-a5.
