\name{predict.locfit}
\alias{predict.locfit}
\title{
Prediction from a Locfit object.
}
\usage{
predict(object, newdata, se.fit=F, \dots)
}
\arguments{
\item{object}{
Fitted object from \code{\link{locfit}()}.
}
\item{newdata}{
Points to predict at. Can be given in several forms: vector/matrix;
list, data frame.
}
\item{se.fit}{
If \code{TRUE}, standard errors are computed along with the fitted values.
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
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
