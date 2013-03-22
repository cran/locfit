\name{rv}
\alias{rv}
\title{
Residual variance from a "locfit" object.
}
\usage{
rv(fit)
}
\arguments{
\item{fit}{
\code{"locfit"} object.
}}
\examples{
data(ethanol)
fit <- locfit(NOx~E,data=ethanol)
rv(fit)
}
\seealso{
\link{locfit}, \link{rv<-}
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a5.
