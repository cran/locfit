\name{locfit.raw}
\alias{locfit.raw}
\title{
Local Regression, Likelihood and Density Estimation.
}
\usage{
locfit.raw(x, y, weights, \dots)
}
\arguments{
\item{x}{
Independent variable (or matrix).
}
\item{y}{
Response variable.
}
\item{weights}{
Prior weights for observations (reciprocal of variance, or sample size).  
}
\item{cens}{
Censoring indicators for hazard rate or censored regression.  
}
\item{base}{
Baseline parameter estimate.  
}
\item{xlim}{
For density estimation, optional vector of lower and upper
bounds for variables.  
}
\item{flim}{
A vector of lower and upper bounds for the evaluation structure.
Defaults to the data range.  
}
\item{scale}{
A scale to apply to each variable. Effectively, the data is
transformed before fitting, by dividing each component of the 
independent variable by the corresponding component of \code{scale}.  
}
\item{alpha}{
Smoothing parameter. A single number (e.g. \code{alpha=0.7})
is interpreted as a nearest neighbor fraction. With two
componentes (e.g. \code{alpha=c(0.7,1.2)}), the first component
is a nearest neighbor fraction, and the second component is
a fixed component. A third component is the penalty term in locally
adaptive smoothing.
}
\item{ev}{
Evaluation Structure, default = \code{"tree"}. Also available are 
\code{"phull"}, \code{"data"}, \code{"grid"}, \code{"kdtree"}, \code{"kdcenter"}, \code{"crossval"}.
}
\item{deg}{
Degree of local polynomial. Default: 2 (local quadratic).  
}
\item{family}{
Local likelihood family; \code{"gaussian"}; 
\code{"binomial"}; \code{"poisson"}; \code{"gamma"} and \code{"geom"}.
Density and rate estimation families are \code{"dens"}, \code{"rate"} and 
\code{"hazard"} (hazard rate).  A default is selected based on whether or 
not \code{"y"} and \code{"cens"} arguments are given; specifying a 
\code{family} argument overrides the default.  
}
\item{link}{
Link function for local likelihood fitting. Depending on the family, 
choices may be \code{"ident"}, \code{"log"}, \code{"logit"}, \code{"inverse"}, 
\code{"sqrt"}.  
}
\item{maxk}{
Controls space assignment for evaluation structures, default \code{50} 
}
\item{kern}{
Weight function, default = \code{"tcub"}. Others are \code{"rect"}, \code{"trwt"},
\code{"tria"}, \code{"epan"}, \code{"bisq"} and \code{"gauss"}. Choices may be restricted
when derivatives are required; e.g. for confidence bands and some
bandwidth selectors.
}
\item{kt}{
Kernel type, \code{"sph"} (default); \code{"prod"}.
In multivariate problems, \code{"prod"\} uses a 
simplified product model which speeds up computations.
}
\item{itype}{
Integration type for density estimation.
}
\item{mint}{
Points for numerical integration rules. Default 20.
}
\item{maxit}{
Maximum iterations for local likelihood estimation. Default 20.
}
\item{cut}{
Refinement parameter for adaptive partitions. Default 0.8; smaller 
values result in more refined partitions.  
}
\item{dc}{
Derivative adjustment.  
}
\item{geth}{
Don't use! 
}
\item{mg}{
Margin size for grids.  
}
\item{deriv}{
Don't use! 


}}
\value{
An object with class "locfit". A standard set of methods for printing, ploting, etc. these objects is provided.
}
\references{
Consult the Web page \url{http://cm.bell-labs.com/stat/project/locfit/}.
}
\keyword{locfit}
% Converted by Sd2Rd version 0.2-a4.
