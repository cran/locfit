\name{rva}
\alias{rv<-}
\title{
Substitute variance estimate on a "locfit" object.
}
\description{
By default, Locfit uses the normalized residual sum of squares as the
variance estimate when constructing confidence intervals.
In some cases, the user may like to use alternative variance
estimates; this function allows the default value to be changed.
}
\usage{
rv(fit) <- 1.2345
}
\arguments{
\item{fit}{
\code{"locfit"} object.
}}
\seealso{
\link{locfit}(), \link{rv}(), \link{plot.locfit}()
}
\keyword{locfit}
