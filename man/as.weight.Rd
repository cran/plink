\name{as.weight}
\alias{as.weight}
\title{Calibration Weights}
\description{
This function facilitates the creation of either a set of quadrature weights or
weights based on a set of theta values for use in the function 
\code{\link[=plink-methods]{plink} }
}
\usage{
as.weight(n=40, theta1 = NULL, weight1 = NULL, theta2 = NULL, weight2 = NULL, 
  normal.wt = FALSE, ...)
}

\arguments{
  \item{n}{number of quadrature points to use}

  \item{theta1}{vector of theta values for the set of \code{TO} item parameters
    (i.e. the item parameters that will not be transformed)}

  \item{weight1}{vector of weights corresponding to \code{theta1}.  They must
    sum to one.}

  \item{theta2}{vector of theta values for the set of \code{FROM} item parameters
    (i.e. the item parameters that will be transformed)}

  \item{weight2}{vector of weights corresponding to \code{theta2}.  They must
    sum to one.}

  \item{normal.wt}{if \code{TRUE} and \code{weight1} is \code{NULL}, the weights for 
    \code{theta1} and \code{theta2} will be computed to correspond to the
    densities from a normal distribution.}

  \item{...}{further arguments passed to other methods }
}

\details{
If no theta values are supplied, \code{n} normal quadrature points and weights will be 
computed using the \code{\link[statmod:gauss.quad.prob]{gauss.quad.prob}} function from 
the \code{statmod} package. Additional arguments can be passed to this function to change
the distribution used to compute the quadrature weights. 

If values are only supplied for \code{theta1}, the values for \code{theta2} will equal those of
\code{theta1}.  If values are supplied for \code{theta1} and \code{theta2}, the length
of the two vectors must be equal.

If no weights are supplied, and \code{normal.wt}=FALSE, the values for \code{weight1} 
and \code{weight2} will be equal to 1/m for m theta values. If no weights are supplied, 
and \code{normal.wt}=TRUE, the values for \code{weight1} and \code{weight2} will be equal 
to the normal densities corresponding to the theta values.  

If values are only supplied for \code{weight1}, the values for \code{weight2} will equal 
those of \code{weight1}.  If values are supplied for \code{weight1} and/or \code{weight2}, 
the length of the vector(s) should be equal to the length of the theta vector(s). If values 
are supplied for \code{weight1} or \code{weight2}, \code{normal.wt} will be ignored.

See Kolen & Brennan (2004) for more information on calibration weights.
}

\value{
Returns a list of length two.  The first list element is an n x 2 matrix of theta values 
corresponding to \code{theta1} and \code{theta2}.  The second list element is an
n x 2 matrix of weights corresponding to \code{weight1} and \code{weight2}.
}

\references{ 
  Kolen, M. J., & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking.}
  New York: Springer
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\note{When \code{normal.wt} equals \code{TRUE}, \code{theta1} and \code{theta2} are
standardized to determine normal density weights}

\seealso{\code{\link[=plink-methods]{plink}}, \code{\link[statmod:gauss.quad.prob]{gauss.quad.prob}} }

\examples{
# Create a set of 40 (default) theta values and weights using quadrature
wt <- as.weight()

# Create a set of 20 theta values and weights using quadrature
wt <- as.weight(20)

# Create a set of 40 quadrature points and weights using a uniform distribution
wt <- as.weight(dist="uniform")

# Create weights for a set of random normal theta values
wt <- as.weight(theta1=rnorm(100))

# Create an object with equal interval theta values and normal density weights
wt <- as.weight(theta1=seq(-4,4,0.05), normal.wt=TRUE)

}

\keyword{ misc }
