\name{residuals.locfit}
\alias{residuals.locfit}
\title{
Fitted values for a `"locfit"' object.
}
\usage{
residuals(object, data, cv=F, studentize=F, type="dev")
}
\arguments{
\item{object}{
\code{"locfit"} object.
}
\item{data}{
The data frame for the original fit. Usually, shouldn't be needed.
}
\item{cv}{
If \code{TRUE}, leave-one-out cross validated residuals are approximated.
}
\item{studentize}{
If \code{TRUE}, residuals are studentized.
}
\item{type}{
  Type of fit or residuals to compute. The default is \code{"fit"} for
  \code{fitted.locfit}, and \code{"dev"} for \code{\link{residuals.locfit}}.
  Other choices include \code{"pear"} for Pearson residuals; \code{"raw"}
  for raw residuals, \code{"ldot"} for likelihood derivative;
  \code{"d2"} for the deviance residual squared; \code{lddot} for the
  likelihood second derivative. Generally, \code{type} should only be
  used when \code{what="coef"}.
}
}
\description{\code{residuals.locfit} is implemented as a front-end to
  \code{\link{fitted.locfit}}, with the \code{type} argument set.
\value{
A numeric vector of the residuals.
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
