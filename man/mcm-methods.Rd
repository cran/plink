\name{mcm-methods}
\docType{methods}
\alias{mcm}
\alias{mcm-methods}
\alias{mcm,matrix,numeric-method}
\alias{mcm,data.frame,numeric-method}
\alias{mcm,list,numeric-method}
\alias{mcm,irt.pars,ANY-method}
\alias{mcm,sep.pars,ANY-method}
\title{Multiple-Choice Model Response Probabilities}
\description{
This function computes the probability of responding in a specific category for one or more 
items for a given set of theta values using the multiple-choice model or multidimensional
multiple-choice model.}

\usage{
mcm(x, cat, theta, dimensions = 1, items, ...)

## S4 method for signature 'matrix', 'numeric'
mcm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'data.frame', 'numeric'
mcm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'list', 'numeric'
mcm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'irt.pars', 'ANY'
mcm(x, cat, theta, dimensions, items, ...)

## S4 method for signature 'sep.pars', 'ANY'
mcm(x, cat, theta, dimensions, items, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters}

  \item{cat}{vector identifying the number of response categories plus one for each item 
  (the additional category is for 'do not know') }

  \item{theta}{vector, matrix, or list of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.5. See details below for more information. }

  \item{dimensions}{number of modeled dimensions }

  \item{items}{numeric vector identifying the items for which probabilities should be computed}

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
  \describe{

    \item{x = "matrix", cat = "numeric"}{This method allows one to specify an n x h matrix
      for n items and h equal to [m x 2k + (k-1)] where m is the number of modeled dimensions
      and k is equal to the maximum number of response categories (including the 'do not know' 
      category) across items.  The first (m x k) columns are for category slope parameters,
      the next block of (m x k) columns are for category difficulty parameters, and the 
      remaining columns are for the lower asymptote (guessing) parameters.  For any items 
      with fewer categories than the maximum, the remaining cells in each block of (m x k) 
      columns or the last k columns should be \code{NA}.  

      \describe{
        \item{Unidimensional Specification:}{Say we have one item with four actual response
        categories and one item with five response categories.  There will be 17 columns.  
        The first six columns are for the category slope parameters.  The first column should 
        contain the parameters for the 'do not know' category.  Column six for the four category 
        item should be \code{NA}.  The next six columns (7-12) are for category difficulty 
        parameters.  The first column of this subset of columns (column 7) should contain the 
        category difficulties for the 'do not know' category.  Similar to the block of columns 
        containing the slopes, the last column in this subset of columns (column 12) for the 
        four category item should be \code{NA}.  The remaining five columns are for the lower
        asymptote (guessing) parameters.  The last column for the four category item would 
        be \code{NA}. }

        \item{Multidimensional Specification:}{In the multidimensional case, the columns for 
        the slope and difficulty parameters should be grouped first by dimension and then by 
        category.  Using the same example for the two items with two dimensions there will be 
        29 columns.  The first five columns for the four category item would include the slope 
        parameters associated with the first dimension for the 'do not know' category and each 
        of the four actual categories respectively.  Columns 11-12 would be \code{NA}.  Columns 
        13-17 would include the category difficulties associated with the first dimension (again
        the parameters for the 'do not know' category should be in the first column of this 
        block of columns) and columns 23-24 would be \code{NA}.  The remaining five columns are 
        for the lower asymptote (guessing) parameters.  The last column for the four category 
        item would be \code{NA}. }
    }}

    \item{x = "data.frame", cat = "numeric"}{See the method for x = "matrix" }

    \item{x = "list", cat = "numeric"}{This method is for a list with three elements. The 
	first element is an n x (m x k) matrix of category slope values for n items, m 
      dimensions, and k equal to the maximum number of response categories across items
      (including the 'do not know' category).  The second list element is an n x (m x k) 
      matrix of category difficulty parameters and the last element is an n x (k-1) matrix
      of lower asymptote (guessing) parameters.  For any list element, for items with fewer 
      categories than the maximum, the remaining cells in the rows should be \code{NA} (see 
      the examples for method x = "matrix" for specification details). }

    \item{x = "irt.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the mcm items in an object of class \code{"irt.pars"}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the mcm items only.  If \code{x} 
      contains parameters for multiple groups, a list of \code{"irt.prob"} objects will be 
      returned.  The argument \code{dimensions} does not need to be included for this method.}

    \item{x = "sep.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the mcm items in an object of class \code{sep.pars}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the mcm items only.  The argument 
      \code{dimensions} does not need to be included for this method.}
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
  Bolt, D. M. & Johnson, T. J. (in press) 
  Applications of a MIRT model to self-report measures: 
  Addressing score bias and DIF due to individual differences in response style. 
  \emph{Applied Psychological Measurement}.

  Thissen, D., & Steinberg, L. (1984). 
  A response model for multiple choice items.
  \emph{Psychometrika}, \emph{49}(4), 501-519.

  Thissen, D., & Steinberg, L. (1996)
  A response model for multiple choice items.
  In W.J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory}. New York: Springer-Verlag

  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\note{No multidimensional extension of the multiple-choice model has officially been presented in
the literature; however, this model is consistent with how the 3PL and nominal response model
were extended to the multidimensional context.
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
###### Unidimensional Examples ######
## Item parameters from Thissen & Steinberg (1984, p. 510)
## Items R,S,T,U for the whole test
a <- matrix(c(-1.7, -1, 1.1, .3, 1.9, 
  -2.1, -.6, 1.2, 2.3, -.8, 
  -1.3, -.9, -.2, 1.9, .5, 
  -1.9, -.5, 0, -.6, 1.9),4,5,byrow=TRUE)
c <- matrix(c(.3, -2.3, 2.4, -2.5, 2.1, 
  2.1, .05, -3, -.6, 1, 
  -.9, -2.5, -.1, 1.8, 1.6, 
  -.1, -2, .5, .8, .8),4,5,byrow=TRUE)
d <- matrix(c(.25, .25, .25, .25, 
  .2, .2, .4, .2, 
  .2, .2, .4, .2, 
  .25, .25, .25, .25),
  4,4,byrow=TRUE)
pars <- cbind(a,c,d)
x <- mcm(pars, rep(5,4))
plot(x,item.names=paste("Item",c("R","S","T","U")), 
  auto.key=list(space="right"))

## Item parameters from Thissen & Steinberg (1984, p. 511)
## Items W,X,Y,Z for the 
pars <- vector("list",3)
pars[[1]] <- matrix(c(-2.3, -.2, 2, .9, -.3, 
  -.8, .6, -.5, 1.1, -.4, 
  -.5, -.2, 2, -1.2, 0, 
  -1.5, -.7, -.2, .1, 2.3),4,5,byrow=TRUE)
pars[[2]] <- matrix(c(.5, .7, -.5, -1.9, 1.1, 
  1.6, -2.8, 1.5, 0, -.3, 
  -.3, .7, -1, .7, 0, 
  .4, .4, -.5, .5, -.8),4,5,byrow=TRUE)
pars[[3]] <- matrix(c(.2, .4, .2, .2, 
  .2, .2, .4, .2, 
  .2, .4, .2, .2, 
  .2, .2, .2, .4), 4,4,byrow=TRUE)
x <- mcm(pars, rep(5,4))
plot(x,item.names=paste("Item",c("W","X","Y","Z")), 
  auto.key=list(space="right"))


###### Multidimensional Example ######
## Discrimination and category parameters from Bolt & Johnson (in press)
pars <- matrix(c(-1.28, -1.029, -0.537, 0.015, 0.519, 0.969, 1.343,
1.473, -0.585, -0.561, -0.445, -0.741, -0.584, 1.444,
0.29, 0.01, 0.04, 0.34, 0, -0.04, -0.63,
0.01, 0.09, 0.09, 0.28, 0.22, 0.31),1,27)
x <- mcm(pars, cat=7, dimensions=2)
# Plot separated surfaces
plot(x,separate=TRUE,drape=TRUE)
}

\keyword{methods}
\keyword{distribution}
