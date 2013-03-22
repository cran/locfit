\name{preplot.locfit}
\alias{preplot.locfit}
\title{
Prediction from a Locfit object.
}
\usage{
preplot(object, newdata, se.fit=F, deriv, tr, what, get.data=F, f3d=F)
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
\item{deriv}{
vector to control prediction and plotting of derivatives. E.g.
\code{deriv=1} for the local slope; \code{deriv=c(1,1)} for local curvature.
}
\item{tr}{
Transformation for likelihood models. Default is the inverse of
the link function.
}
\item{what}{
What to predict/plot. Default is \code{coef} (the local fit). Other choices
are \code{infl} (influence function); \code{nlx} (length of weight diagram);
\code{band} (bandwidth); \code{degr} (local polynomial degree).
}
\item{get.data}{
If \code{TRUE}, the original data is attached to the returned object, and
added to the plot.
}
\item{f3d}{
If \code{TRUE}, forces a trellis plot.
}}
\value{
An object with class \code{"locfitpred"}, containing the predicted
values and additional information used to construct the plot.
}
\seealso{
\code{\link{locfit}}, \code{\link{predict.locfit}}, \code{\link{plot.locfit}}.
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
