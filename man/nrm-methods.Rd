\name{nrm-methods}
\docType{methods}
\alias{nrm}
\alias{nrm-methods}
\alias{nrm,matrix,numeric-method}
\alias{nrm,data.frame,numeric-method}
\alias{nrm,list,numeric-method}
\alias{nrm,irt.pars,ANY-method}
\alias{nrm,sep.pars,ANY-method}
\title{Nominal Response Model Probabilities}
\description{
This function computes the probability of responding in a specific category for one or more 
items for a given set of theta values using the nominal response model or multidimensional
nominal response model.
}

\usage{
nrm(x, cat, theta, dimensions = 1, items, ...)

## S4 method for signature 'matrix', 'numeric'
nrm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'data.frame', 'numeric'
nrm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'list', 'numeric'
nrm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'irt.pars', 'ANY'
nrm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'sep.pars', 'ANY'
nrm(x, cat, theta, dimensions, items, ...)
}

\arguments{
  \item{x}{Object containing item parameters. See below for more details. }

  \item{cat}{vector identifying the number of response categories for each item } 

  \item{theta}{vector, matrix, or list of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.5. See details below for more information. }

  \item{dimensions}{number of modeled dimensions }

  \item{items}{numeric vector identifying the items for which probabilities should be computed}

  \item{...}{further arguments passed to or from other methods }
}


\section{Methods}{
  \describe{

    \item{x = "matrix", cat = "numeric"}{This method allows one to specify an n x (m x 2k) 
      matrix for n items, m dimensions, and k equal to the maximum number of response 
      categories across items.  The first (m x k) columns are for category slope parameters
      and the remaining columns are for the category difficulty parameters.  For any items 
      with fewer categories than the maximum, the remaining cells in each block of (m x k) 
      columns should be \code{NA}.  

      \describe{
        \item{Unidimensional Specification:}{Say we have one four category item and one 
        five category item, the first four columns of the four response item would include 
        the slope parameters. The fifth column for this item would be \code{NA}.  The next 
        four columns would include the category difficulty values, and the last column would 
        be \code{NA}. }

        \item{Multidimensional Specification:}{In the multidimensional case, the columns for 
        the slope and difficulty parameters should be grouped first by dimension and then by 
        category.  Using the same example for the two items with two dimensions there will be 
        20 columns.  The first four columns for the four category item would include the slope 
        parameters associated with the first dimension for each of the four categories 
        respectively.  Columns 9-10 would be \code{NA}.  Columns 11-14 would include the 
        category difficulties associated with the first dimension and columns 19-20 would 
        be \code{NA}. }
    }
}

    \item{x = "data.frame", cat = "numeric"}{See the method for x = "matrix" }

    \item{x = "list", cat = "numeric"}{This method is for a list with two elements. The 
	first element is an n x (m x k) matrix of category slope values for n items, m 
      dimensions, and k equal to the maximum number of response categories across items.
      The second list element is an n x (m x k) matrix of category difficulty parameters. 
      For either element, for items with fewer categories than the maximum, the remaining 
      cells in the rows should be \code{NA} (see the examples for method x = "matrix" for
      specification details). }

    \item{x = "irt.pars", cat = "ANY"}{This method can be used to compute probabilities 
	for the nrm items in an object of class \code{"irt.pars"}. If \code{x} contains 
	dichotomous items or items associated with another polytomous model, a warning 
	will be displayed stating that probabilities will be computed for the nrm
	items only. If \code{x} contains parameters for multiple groups, a list of 
	\code{"irt.prob"} objects will be returned.}

    \item{x = "sep.pars", cat = "ANY"}{This method can be used to compute probabilities 
	for the mcm items in an object of class \code{sep.pars}. If \code{x} contains 
	dichotomous items or items associated with another polytomous model, a warning 
	will be displayed stating that probabilities will be computed for the nrm
	items only.}
  }
}

\details{
\code{theta} can be specified as a vector, matrix, or list.  For the unidimensional case, \code{theta}
should be a vector. If a matrix or list of values is supplied, they will be converted to a single vector
of theta values.  For the multidimensional case, if a vector of values is supplied it will be assumed
that this same set of values should be used for each dimension.  Probabilities will be computed for each
combination of theta values.  Similarly, if a list is supplied, probabilities will be computed for each
combination of theta values.  In instances where probabilities are desired for specific combinations of 
theta values, a j x m matrix should be specified for j ability points and m dimensions where the columns 
are ordered from dimension 1 to m.
}

\value{
  Returns an object of class \code{\linkS4class{irt.prob}}
}

\references{ 
  Bock, R.D. (1972)
  Estimating item parameters and latent ability when responses are scored in two or more 
  nominal categories.
  \emph{Psychometrika}, \emph{37}(1), 29-51.

  Bock, R.D. (1996)
  The nominal categories model.
  In W.J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory}. New York: Springer-Verlag

  Bolt, D. M. & Johnson, T. J. (in press)
  Applications of a MIRT model to self-report measures: 
  Addressing score bias and DIF due to individual differences in response style. 
  \emph{Applied Psychological Measurement}.

  Kolen, M. J., & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking.}
  New York: Springer

  Takane, Y., & De Leeuw, J. (1987)
  On the relationship between item response theory and factor analysis of discretized variables. 
  \emph{Psychometrika}, \emph{52}(3), 393-408.

  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
###### Unidimensional Example ######
## Item parameters from Bock (1972, p. 46,47)
a <- matrix(c(.905, .522, -.469, -.959, NA, 
  .828, .375, -.357, -.079, -.817), 2,5,byrow=TRUE)
c <- matrix(c(.126, -.206, -.257, .336, NA, 
  .565, .865, -1.186, -1.199, .993), 2,5,byrow=TRUE)
pars <- cbind(a,c)
x <- nrm(pars, c(4,5))
plot(x,auto.key=list(space="right"))

###### Multidimensional Example ######
# From Bolt & Johnson (in press)
pars <- matrix(c(-1.28, -1.029, -0.537, 0.015, 0.519, 0.969, 1.343,
1.473, -0.585, -0.561, -0.445, -0.741, -0.584, 1.444,
0.29, 0.01, 0.04, 0.34, 0, -0.04, -0.63), 1,21)
x <- nrm(pars, cat=7, dimensions=2)
# Plot separated surfaces
plot(x,separate=TRUE,drape=TRUE)
}

\keyword{methods}
\keyword{distribution}
