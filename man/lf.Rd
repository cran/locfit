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
  \item{alpha}{
  Smoothing parameter. A single number (e.g. \code{alpha=0.7})
  is interpreted as a nearest neighbor fraction. With two
  componentes (e.g. \code{alpha=c(0.7,1.2)}), the first component
  is a nearest neighbor fraction, and the second component is
  a fixed component. A third component is the penalty term in locally
  adaptive smoothing.}
  \item{deg}{degree of polynomials to be fitted locally.}
  \item{scale}{A scale to apply to each variable. This is especially important for
  multivariate fitting, where variables may be measured in
  non-comparable units. It is also used to specify the frequency
  for \code{\link{ang}} terms. If \code{scale=FALSE} (the default) no scaling
  is performed. If \code{scale=TRUE}, marginal standard deviations are used.
  Alternatively, a numeric vector can provide scales for the
  individual variables.}
  \item{kern}{Weight function, default = \code{"tcub"}.
  Other choices are \code{"rect"}, \code{"trwt"}, \code{"tria"},
  \code{"epan"}, \code{"bisq"} and \code{"gauss"}. Choices may be restricted
  when derivatives are required; e.g. for confidence bands and some
  bandwidth selectors.}
  \item{ev}{Evaluation Structure, default = \code{"tree"}. Also available are
  \code{"phull"}, \code{"data"}, \code{"grid"}, \code{"kdtree"},
  \code{"kdcenter"} and \code{"crossval"}. \code{ev="none"} gives no
  evaluation points, effectively producing the global parametric fit.
  A vector or matrix of evaluation points can also be provided.}
  \item{maxk}{Controls space assignment for evaluation structures.
  For the adaptive evaluation structures, it is impossible to be sure
  in advance how many vertices will be generated. If you get
  warnings about `Insufficient vertex space', Locfit's default assigment
  can be increased by increasing \code{maxk}. The default is \code{maxk=100}.}

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
