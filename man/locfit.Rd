\name{locfit}
\alias{locfit}
\title{
Local Regression, Likelihood and Density Estimation.
}
\usage{
locfit(formula, data, weights, cens, base, subset, geth, \dots, lfproc)
}
\arguments{
\item{formula}{
Model Formula; e.g. \code{y~x} for a regression model; \code{~x} for a
density estimation model
}
\item{data}{
Data Frame.
}
\item{weights}{
Prior weights (or sample sizes) for individual observations. This is
typically used where observations have unequal variance.
}
\item{cens}{
Censoring indicator. \code{1} (or \code{TRUE}) denotes a censored observation.
\code{0} (or \code{FALSE}) denotes uncensored.
}
\item{base}{
Baseline for local fitting. For local regression models, specifying
a \code{base} is equivalent to using \code{y-base} as the reponse. But \code{base}
also works for local likelihood.
}
\item{subset}{
Subset observations in the data frame.
}
\item{geth}{
Don't use.
}
\item{...}{
Other arguments to \code{\link{locfit.raw}()} (or the \code{lfproc}).
}
\item{lfproc}{
A processing function to compute the local fit. Default is
\code{locfit.raw()}. Other choices include \code{locfit.robust()},
\code{locfit.censor()} and \code{locfit.quasi()}.
}}
\value{
An object with class \code{"locfit"}. A standard set of methods for printing,
ploting, etc. these objects is provided.
\code{locfit()} is a model-formula based front end to \code{\link{locfit.raw}()}.
}
\seealso{
\code{\link{locfit.raw}}
}
\examples{
data(ethanol)
locfit(NOx~E,data=ethanol)
}
\references{
Consult the Web page \url{http://cm.bell-labs.com/stat/project/locfit/}.
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
