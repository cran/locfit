\name{ang}
\alias{ang}
\title{
Angular Term for a Locfit model.
}
\usage{
locfit(y~ang(x), scale=1/(2*pi))
}
\description{
  The \code{ang()} function is used in a locfit model formula
  to specify that a variable should be treated as an angular
  or periodic term. The \code{scale} argument of
  \code{\link{locfit}()} is used to set the period.
}

\arguments{
  \item{x}{numeric variable to be treated periodically.}
  \item{scale}{Locfit's \code{scale} argument is used to specify
    the period divided by \eqn{2\pi} of the term. The default is
    \code{sacle=1}, giving a period of  \eqn{2\pi}.}
}

\examples{
# generate an x variable, and a response with period 0.2
x <- seq(0,1,length=200)
y <- sin(10*pi*x)+rnorm(200)/5

# compute the periodic local fit. Note the scale argument is period/(2pi)
fit <- locfit(y~ang(x),scale=0.2/(2*pi))

# plot the fit over a single period
plot(fit)

# plot the fit over the full range of the data
plot(fit,xlim=c(0,1))
}

\references{
Loader, C. (1999). Local Regression and Likelihood. Springer, NY (Section 6.2).
}

\seealso{
  \code{\link{locfit}}.
}

\keyword{models}
% Converted by Sd2Rd version 0.2-a5.
