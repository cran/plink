\name{as.weight}
\alias{as.weight}
\title{Calibration Weights}
\description{
This function facilitates the creation of either a set of quadrature weights or
weights based on a set of theta values for use in the function 
\code{\link[=plink-methods]{plink} }
}
\usage{
as.weight(n, theta, weight, quadrature = FALSE, normal.wt = FALSE, 
  dimensions = 1, ...)
}

\arguments{
  \item{n}{numeric value or vector identifying the number of theta values
  to use for each dimension. If only one value is supplied but \code{dimensions} is
  greater than one, \code{n} will be the same for all dimensions.}

  \item{theta}{vector or list of theta values. If \code{dimensions} is greater than
  one, list elements should correspond to each dimension, otherwise the single vector
  will be repeated for all dimensions.}

  \item{weight}{vector or list of weights. If \code{dimensions} is greater than
  one, list elements should correspond to each dimension, otherwise the single vector
  will be repeated for all dimensions.}

  \item{quadrature}{if \code{TRUE}, quadrature points and weights will be 
    computed for the corresponding number of points \code{n} for each dimension
    respectively. See \code{\link[statmod:gauss.quad.prob]{gauss.quad.prob}} 
    for more information on creating quadrature points and weights.}

  \item{normal.wt}{if \code{TRUE} and \code{weight} is \code{missing}, the weights 
    for \code{theta} will be computed to correspond to the densities from a 
    normal distribution. }

  \item{dimensions}{number of dimensions for which the weights should be created}

  \item{...}{further arguments passed to other methods }
}

\details{
When weighting expected response probabilities at different theta values using
characteristic curve linking methods, there are a variety of approaches one can take.
These range from uniform weights to normal weights, to quadrature weights corresponding
to some a priori specified distribution. The purpose of this function is to facilitate 
the creation of these weights for use in \code{\link[=plink-methods]{plink} }. 

For all approaches, when more than one dimension is specified, the weights for each combined 
set of theta values will be a multiplicative weight. For example, if there are two dimensions 
and the specified weights corresponding to two specific theta values on each dimension
respectively are 0.8 and 1.2, the final weight for this pair of theta values will be 0.96.

\describe{
  \item{Uniform Weights}{Five general approaches can be used to create uniform weights.}
    \itemize{
      \item{}{If no arguments are supplied, a set of weights (all equal to one) will be 
      returned for a single dimension, for 40 equal interval theta values ranging from 
      -4 to 4. If \code{dimensions} is greater than one, seven equal interval theta values
      ranging from -4 to 4 will be specified for each dimension. For instance, for two
      dimensions, there will be weights for 7 x 7 (49) points.}

      \item{}{If only a value for \code{n} is supplied, uniform weights (all equal to one)
      will be created for \code{n} points ranging from -4 to 4 for each dimension specified.}

      \item{}{If values are only specified for \code{theta}, uniform weights (all equal to one) 
      will be created for each of these values for each dimension specified.}
      
      \item{}{If values are only specified for \code{weight} where the values are all equal. In 
      this case, equal interval theta values will be selected from -4 to 4 to correspond to
      the number of specified weights.}
      
      \item{}{If values are specified for \code{theta} and uniform values are specified for 
      \code{weight}.}
  }

  \item{Non-Uniform Weights}{Four general approaches can be used to create non-uniform weights.
    \itemize{
      \item{}{If values are only specified for \code{weight} where the values are not equal. In 
      this case, equal interval theta values will be selected from -4 to 4 to correspond to
      the number of specified weights. } 

      \item{}{If values are specified for \code{theta} and varying values are specified for 
      \code{weight}. }

      \item{}{If \code{quadrature} is equal to \code{TRUE}, and no other arguments are specified
      \code{n} will default to 40. if \code{dimensions} is greater than one, \code{n} will default
      to seven for each dimension. In either case \code{n} quadrature points and weights will be 
      selected from a standard normal distribution. To change the mean and/or standard deviation of 
      this distribution, values for \code{mu} and \code{sigma} respectively should be specified. See
      \code{\link[statmod:gauss.quad.prob]{gauss.quad.prob}} for more information. Different means
      and/or SDs can be supplied for each dimension. If values are specified for \code{theta} or 
      \code{weight}, the quadrature approach will not be used. }

      \item{}{If \code{quadrature} equals \code{TRUE} other distributions can be specified for
      \code{n} points. See \code{\link[statmod:gauss.quad.prob]{gauss.quad.prob}} for more 
      information.}

      \item{}{If \code{normal.wt} equals \code{TRUE}, normal densities will be created for the 
      specified \code{theta} values (if supplied) or equal interval values. The default distribution 
      will be standard normal, but different means and/or standard deviations can be specified by
      passing arguments for \code{mean} and \code{sd} respectively. Different means and/or SDs can be 
      supplied for each dimension. If no values are included for \code{theta}, equal interval theta 
      values will be created for the range of three SDs above and below the mean. If values are 
      specified for \code{weight}, the \code{normal.wt} argument will be ignored. }
    }
  }
}


See Kolen & Brennan (2004) for more information on calibration weights.
}

\value{
Returns a list of length two.  The first list element is an n x m matrix of theta values 
for m dimensions.  The second list element is a vector of weights corresponding to the 
theta values in the first list element.
}

\references{ 
  Kolen, M. J., & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking.}
  New York: Springer
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{\code{\link[=plink-methods]{plink}}, \code{\link[statmod:gauss.quad.prob]{gauss.quad.prob}} }

\examples{
###### Unidimensional Examples ######
# Create a set of 40 (default) theta values and uniform weights
wt <- as.weight()

# Create a set of 40 quadrature points and weights using a standard 
# normal distribution
wt <- as.weight(quadrature=TRUE)

# Create a set of 30 quadrature points and weights using a normal 
# distribution with a mean of 0.5
wt <- as.weight(n=30, quadrature=TRUE, mu=0.5)

# Create weights for a set of random normal theta values
wt <- as.weight(theta=rnorm(100))

# Create an object with equal interval theta values and standard 
# normal density weights
wt <- as.weight(theta=seq(-4,4,0.05), normal.wt=TRUE)

# Create an object with equal interval theta values and normal 
# density weights with a mean of 0.5 and SD equal to .92
wt <- as.weight(theta=seq(-4,4,0.05), normal.wt=TRUE, mean=0.5, sd=0.92)

###### Multidimensional Examples ######
# Create a set of 49 theta values and uniform weights 
# (based on seven points for each dimension)
wt <- as.weight(dimensions=2)

# Create a set of 100 quadrature points and weights using a normal 
# distribution with a means of 0 and 0.5 for the two dimensions respectively
wt <- as.weight(n=10, quadrature=TRUE, mu=c(0,0.5), dimensions=2)

# Create an object with equal interval theta values and standard 
# normal density weights for three dimensions
wt <- as.weight(theta=seq(-3,3), normal.wt=TRUE, dimensions=3)

# Create an object with two sets of equal interval theta values for 
# two dimensions
wt <- as.weight(theta=list(seq(-4,4),seq(-3,3)), dimensions=2)

# Create an object with two sets of random normal theta values and 
# standard normal density weights for two dimensions
wt <- as.weight(theta=list(rnorm(10),rnorm(10)), normal.wt=TRUE, dimensions=2)
}

\keyword{ misc }
