\name{predict.locfit}
\alias{predict.locfit}
\title{
Prediction from a Locfit object.
}
\usage{
\method{predict}{locfit}(object, newdata, where="fitp", se.fit=FALSE,
               band="none", what="coef", ...)
}
\description{
  The \code{\link{locfit}} function computes a local fit at a selected set
of points (as defined by the \code{ev} argument). The \code{predict.locfit}
function is used to interpolate from these points to any other points.
The method is based on cubic hermite polynomial interpolation, using the
estimates and local slopes at each fit point.

  The motivation for this two-step procedure is computational speed.
Depending on the sample size, dimension and fitting procedure, the
local fitting method can be expensive, and it is desirable to keep the
number of points at which the direct fit is computed to a minimum.
The interpolation method used by \code{predict.locfit()} is usually
much faster, and can be computed at larger numbers of points.
}
\arguments{
\item{object}{
Fitted object from \code{\link{locfit}()}.
}
\item{newdata}{
Points to predict at. Can be given in several forms: vector/matrix;
list, data frame.
}
\item{where}{
An alternative to \texttt{newdata}. Choices include \code{"grid"}
for the grid \code{\link{lfmarg}(object)}; \code{"data"} for the original
data points and \code{"fitp"} for the direct fitting points (ie. no
interpolation).
}
\item{se.fit}{
If \code{TRUE}, standard errors are computed along with the fitted values.
}
\item{band}{
  Compute standard errors for the fit and include confidence bands on
  the returned object. Default is \code{"none"}. Other
  choices include \code{"global"} for bands using a global variance estimate;
  \code{"local"} for bands using a local variance estimate and \code{"pred"}
  for prediction bands (at present, using a global variance estimate).
  To obtain the global variance estimate for a fit, use \code{\link{rv}}.
  This can be changed with \code{\link{rv<-}}. Confidence bands, by default,
  are 95%, based on normal approximations and neglecting bias.
  To change the critical value or confidence level, or to obtain
  simultaneous instead of pointwise confidence,
  the critical value stored on the fit must be changed. See the
  \code{\link{kappa0}} and \code{\link{crit}} functions.
}
\item{what}{
What to compute predicted values of. The default, \code{what="coef"}, works
  with the fitted curve itself. Other choices include \code{"nlx"} for the
  length of the weight diagram; \code{"infl"} for the influence function;
  \code{"band"} for the bandwidth; \code{"degr"} for the local polynomial
  degree; \code{"lik"} for the maximized local likelihood; \code{"rdf"}
  for the local residual degrees of freedom and \code{"vari"} for the
  variance function. The interpolation algorithm for some of these quantities
  is questionable.
}
\item{...}{
Additional arguments to \code{\link{preplot.locfit}}.
}}
\value{
If \code{se.fit=F}, a numeric vector of predictors.
If \code{se.fit=T}, a list with components \code{fit}, \code{se.fit} and
\code{residual.scale}.
}
\examples{
data(ethanol)
fit <- locfit(NOx~E,data=ethanol)
predict(fit,c(0.6,0.8,1.0))
}
\keyword{smooth}
