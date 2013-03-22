\name{fitted.locfit}
\alias{fitted.locfit}
\title{
Fitted values for a `"locfit"' object.
}
\usage{
fitted(object, data, what, se.fit, cv. tr)
}
\arguments{
\item{object}{
\code{"locfit"} object.
}
\item{data}{
The data frame for the original fit. Usually, shouldn't be needed.
}
\item{what}{
What to compute fitted values of. See \code{preplot.locfit()}.
}
\item{se.fit}{
If \code{TRUE}, standard errors of fit are also computed.
}
\item{cv}{
If \code{TRUE}, leave-one-out cross validated fitted values are approximated.
}
\item{tr}{
Back transformation for likelihood models.
}}
\value{
A numeric vector of the fitted values.
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
