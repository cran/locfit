\name{residuals.locfit}
\alias{residuals.locfit}
\title{
Residuals from a "locfit" object.
}
\usage{
residuals(object, data, type, studentize)
}
\arguments{
\item{object}{
Fitted object from \code{\link{locfit}()}.
}
\item{data}{
The data frame for the original data. Usually this isn't needed.
}
\item{type}{
Type of residuals for local likelihood. Default \code{"deviance"}.
Other choices include \code{"response"}; \code{"pearson"};
\code{"ldot"} (likelihood derivative).
}
\item{studentize}{
Studentize the residuals.
}
\value{
A vector of residuals.
}
\examples{
data(ethanol)
fit <- locfit(NOx~E,data=ethanol)
residuals(fit)
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
