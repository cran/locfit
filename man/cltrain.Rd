\name{cltrain}
\title{Training dataset for classification}
\usage{data(cltrain)}
\alias{cltrain}
\format{
Data Frame. Three variables x1, x2 and y. The latter indicates
class membership.
}
\description{
200 observations from a 2 population model. Under population 0,
\(x_{1,i}\) has a standard normal distribution, and
\(x_{2,i} = (2-x_{1,i}^2+z_i)/3\), where \(z_i\) is also standard normal.
Under population 1, \(x_{2,i} = -(2-x_{1,i}^2+z_i)/3\).
The optimal classification regions form a checkerboard pattern,
with horizontal boundary at \(x_2=0\), vertical boundaries at
\(x_1 = \pm \sqrt{2}\).

This is the same model as the cltest dataset.
}
\keyword{datasets}
