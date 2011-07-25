\name{drm-methods}
\docType{methods}
\alias{drm}
\alias{drm-methods}
\alias{drm,numeric-method}
\alias{drm,matrix-method}
\alias{drm,data.frame-method}
\alias{drm,list-method}
\alias{drm,irt.pars-method}
\alias{drm,sep.pars-method}
\title{Dichotomous Response Model Probabilities}
\description{
This function computes the probability of correct responses (and optionally, incorrect 
responses) for one or more items for a given set of theta values using the 1PL, 2PL, 3PL,
M1PL, M2PL or M3PL model, depending on the included item parameters and the specified
number of dimensions.
}

\usage{
drm(x, theta, dimensions = 1, D = 1, incorrect = FALSE, 
  print.mod = FALSE, items, ...)

## S4 method for signature 'numeric'
drm(x, theta, dimensions, D, incorrect, print.mod, items, ...)

## S4 method for signature 'matrix'
drm(x, theta, dimensions, D, incorrect, print.mod, items, ...)

## S4 method for signature 'data.frame'
drm(x, theta, dimensions, D, incorrect, print.mod, items, ...)

## S4 method for signature 'list'
drm(x, theta, dimensions, D, incorrect, print.mod, items, ...)

## S4 method for signature 'irt.pars'
drm(x, theta, dimensions, D, incorrect, print.mod, items, ...)

## S4 method for signature 'sep.pars'
drm(x, theta, dimensions, D, incorrect, print.mod, items, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters}

  \item{theta}{vector, matrix, or list of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.5. See details below for more information. }

  \item{dimensions}{number of modeled dimensions }

  \item{D}{scaling constant. The default value assumes that the parameters are already in the 
  desired metric. If the parameters are in the logistic metric, they can be transformed to a 
  normal metric by setting D = 1.7}

  \item{incorrect}{if \code{TRUE}, compute probabilities for the incorrect response category) }

  \item{print.mod}{if \code{TRUE}, print the model that was used (i.e., 1PL, 2PL, 3PL, M1PL, M2PL, or M3PL) }

  \item{items}{numeric vector identifying the items for which probabilities should be computed}

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
  \describe{

    \item{x = "numeric"}{This method should only be used for the Rasch model where \code{x} is a vector of 
      difficulty parameters (or parameters related to item difficulty in the multidimensional case).  Under 
      this method the slopes and lower asymptote values for all items will default to one and zero respectively.  
      This is true for both the unidimensional and multidimensional case.}

    \item{x = "matrix"}{This method allows one to specify an n x k matrix for n items, for m
      dimensions, and k equal to m, m+1, or m+2.  The first m columns generally correspond to the 
      discrimination/slope parameters, and the last two columns correspond to the difficulty, and
      lower asymptote parameters, although this may not be the case depending on the supplied
      parameters.
        \describe{
          \item{1PL/M1PL:}{For the 1PL and M1PL model with discriminations equal to 1 (Rasch Model), 
          an n x 1 matrix of item difficulties can be supplied.  An n x (m+1) matrix can also be used with
          all the values in the first m columns equal to 1 and difficulty parameters in the last
          column.  For discrimination values other than 1, \code{x} should include at least m+1
          columns where the first m columns contain the item discriminations (identical for all items) 
          and the last column is for the item difficulties.  The lower asymptote defaults to zero for 
          all items; however, m+2 columns can be specified where the values in the last column all 
          equal zero.}

          \item{2PL/M2PL:}{\code{x} should include at least m+1 columns where the first m columns contain 
          the item discriminations and the last column is for the item difficulties.  The lower asymptote 
          defaults to zero for all items; however, m+2 columns can be specified where the values in the 
          last column all equal zero.}

          \item{3PL/M3PL:}{\code{x} should include m+2 columns where the first m columns contain the item 
          discriminations, the m+1 column is for item difficulties, and the last column is for lower 
          asymptote values.}
	}
    }

    \item{x = "data.frame"}{See the method for x = "matrix" }

    \item{x = "list"}{This method can include a list with 1, 2, or 3 elements.  In general, these 
      elements correspond to discrimination, difficulty, and lower asymptote parameters, although 
      this may not be the case depending on the supplied parameters.  For the unidimensional case,
      each list element can be a vector of length n or an n x 1 matrix for n items.  For the 
      multidimensional case, the list element for the discrimination parameters should be an n x m
      matrix for m dimensions.  The other list elements can be vectors of length n or an n x 1 matrix 
      for n items.
        \describe{
          \item{1PL/M1PL:}{For the 1PL and M1PL model with discriminations equal to 1 (Rasch Model), 
          one element with item difficulties can be supplied.  Alternatively, two elements can be used 
          with the first list element containing a matrix/vector of ones and difficulty parameters in 
          the second list element.  For discrimination values other than 1, \code{x} should contain 
          at  least two list elements where the first contains the item discriminations (identical 
          for all items) and the second is for item difficulties.  The lower asymptote defaults to 
          zero for all items; however, a third element with a vector/matrix of zeros can be included.}

          \item{2PL/M2PL:}{\code{x} should contain at least two list elements where the first element 
          contains the item discriminations and the second element includes the item difficulties.  
          The lower asymptote defaults to zero for all items; however, a third element with a 
          vector/matrix of zeros can be included.}

          \item{3PL/M3PL:}{\code{x} should include three list elements where the first element contains the 
          item discriminations, the second element includes the item difficulties, and the third element
          contains the lower asymptote values.}
	}
    }

    \item{x = "irt.pars"}{This method can be used to compute probabilities for the dichotomous
      items in an object of class \code{"irt.pars"}.  If \code{x} contains polytomous items, a
      warning will be displayed stating that probabilities will be computed for the 
      dichotomous items only.  If \code{x} contains parameters for multiple groups, a list of 
      \code{"irt.prob"} objects will be returned.  The argument \code{dimensions} does not need
      to be included for this method.}

    \item{x = "sep.pars"}{This method can be used to compute probabilities for the dichotomous
      items in an object of class \code{sep.pars}.  If \code{x} contains polytomous items, a
      warning will be displayed stating that probabilities will be computed for the 
      dichotomous items only.  The argument \code{dimensions} does not need to be included for 
      this method.}
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
  Adams, R. J., Wilson, M., & Wang, W. (1997).
  The multidimensional random coefficients multinomial logit model.
  \emph{Applied Psychological Measurement}, \emph{21}(1), 1-23.

  Birnbaum, A. (1968).
  Some latent trait models and their use in inferring an examinee's ability
  In F. M. Lord & M. R. Novick (Eds.)
  \emph{Statistical Theories of Mental Test Scores.}
  Reading, MA: Addison-Wesley.

  Rasch, G. (1960).
  \emph{Probabilistic Models for Some Intelligence and Attainment Tests}
  Copenhagen, Denmark: Danish Institute for Educational Research.

  Reckase, M. D. (1985). 
  The difficulty of items that measure more than one ability.
  \emph{Applied Psychological Measurement}, \emph{9}(4), 401-412. 

  Reckase, M. D. (1997). 
  A linear logistic multidimensional model for dichotomous item response data. 
  In W. J. van der Linden & R. K. Hambleton (Eds.)
  \emph{Handbook of Modern Item Response Theory} (pp. 271-286). 
  New York: Springer-Verlag.

  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\note{ 
The indices 0 and 1 are used in the column labels for incorrect and correct responses 
respectively (e.g. the label d1.0 indicates the incorrect response column for item 1, d1.1 
indicates the correct response column for the same item).  If \code{incorrect} = \code{FALSE}, 
all column labels will end with .1
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
###### Unidimensional Examples ######
## 1PL
# A single item with a difficulty at 1 logit
x <- drm(1) 
plot(x)

# Three items with a difficulties at -2.2, -1.7, and 0.54 logits 
# respectively. Assume a constant discrimination of 1 for each item
b <- c(-2.2,-.7,1.54)
x <- drm(b) 
# Plot the expected probabilites
plot(x)

# Five items with a constant discrimination equal to 1.2
# Use the scaling constant for the normal metric
a <- rep(1.2,5)
b <- c(-2,-1,0,1,2)
pars <- list(a,b)
x <- drm(pars,D=1.7)
# Extract the expected probabilities
get.prob(x)
# Plot the expected probabilities
plot(x,combine=5,item.names="Items 1-5")


## 2PL
# Compute probabilities for five items at a theta value of 1.5 logits
# Print the model
a <- c(1,.8,.6,1.2,.9)
b <- c(-1.7,2.2,.4,-1,-.4)
pars <- cbind(a,b)
drm(pars,theta=1.5,print.mod=TRUE)


# Include a vector of zeros for the lower asymptote
# Compute probabilities for incorrect responses
a <- c(.63,1.14,.89,1.01,1.51,.58)
b <- c(-1.9,.08,1.6,-1.4,.5,-2.3)
c <- rep(0,6)
pars <- cbind(a,b,c)
x <- drm(pars,incorrect=TRUE)
plot(x)


## 3PL
a <- c(1,.8,.4,1.5,.9)
b <- c(-2,-1,0,1,2)
c <- c(.2,.25,.18,.2,.22)
pars <- list(a,b,c)
x <- drm(pars)
plot(x)

# Use theta values from -3 to 3 with an increment of 0.2
a <- c(.71,.96,.36,1.05,1.76,.64)
b <- c(-.16,1.18,2.03,1.09,.82,-1.56)
c <- c(.22,.17,.24,.08,.20,.13)
theta <- seq(-3,3,.2)
pars <- cbind(a,b,c)
x <- drm(pars,theta)
plot(x,combine=6,item.names="Items 1-6",auto.key=list(space="right"))

###### Multidimensional Examples ######
## M1PL
# A single item with a parameter related to item difficulty at 1 logit
x <- drm(1, dimensions=2) 
plot(x)

# Three items with a difficulties at -2.2, -1.7, and 0.54 logits 
# respectively. Assume a constant discrimination of 1 for each item
d <- c(-2.2,-.7,1.54)
x <- drm(d, dimensions=2) 
plot(x, drape=TRUE)


## M2PL
# Items 27-30 from Reckase (1985)
a <- matrix(c(1.66,1.72,.69,.19,.88,1.12,.68,1.21),4,2,byrow=TRUE)
d <- c(-.38,-.68,-.91,-1.08)
pars <- list(a,d)
x <- drm(pars,dimensions=2)
plot(x, drape=TRUE, item.names=paste("Item",27:30))

# Create contourplots for these items
plot(x,type="contourplot",cuts=10)


## M3PL
# Single item from Reckase (1997, p. 274)
pars <- t(c(.8,1.4,-2,.2))
x <- drm(pars, dimensions=2)
plot(x, default.scales=list(arrows=FALSE),drape=TRUE)

# Compute the probabilities of an incorrect response
# for the Reckase (1997) item
x <- drm(pars, dimensions=2, incorrect=TRUE)
plot(x, screen=list(z=-40,x=-60), auto.key=list(space="right"))

# Four items from the included example for BMIRT (Yao, 2003)
# modeled using four dimensions
pars <- matrix(c(0.5038,2.1910,1.1317,0.2493,0.5240,0.1763, 
0.2252,1.1999,0.5997,0.2087, 2.1730,0.4576,
0.2167,0.2487,1.4009,0.3865,-1.5270,0.3507, 
2.3428,1.1530,0.3577,0.4240,-1.4971,0.3641),4,6,byrow=TRUE)
colnames(pars) <- c("a1","a2","a3","a4","d","c")
x <- drm(pars, dimensions=4, print.mod=TRUE)
# Plot the item response surfaces for item 4
\dontrun{
plot(x,items=4,item.names="Item 4",drape=TRUE)
}
}

\keyword{methods}
\keyword{distribution}
