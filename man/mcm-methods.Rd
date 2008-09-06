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
items for a given set of theta values using the multiple-choice model.}

\usage{
mcm(x, cat, theta = seq(-4,4,0.05), ...)

## S4 method for signature 'matrix', 'numeric'
mcm(x, cat, theta, ...)

## S4 method for signature 'data.frame', 'numeric'
mcm(x, cat, theta, ...)

## S4 method for signature 'list', 'numeric'
mcm(x, cat, theta, ...)

## S4 method for signature 'irt.pars', 'ANY'
mcm(x, cat, theta, ...)

## S4 method for signature 'sep.pars', 'ANY'
mcm(x, cat, theta, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters}

  \item{cat}{vector identifying the number of response categories plus one for each item 
  (the additional category is for 'do not know') }

  \item{theta}{vector of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.05 }

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
  \describe{

    \item{x = "matrix", cat = "numeric"}{This method allows one to specify an n x k matrix
      for n items and k equal to three times the maximum number of response categories 
      (including the 'do not know' category) across items minus one (e.g. if some items 
      have four response categories and others have five categories, k will equal 14). 

      \describe{
        \item{Discrimination Parameters:}{The first (k+1)/3 columns are for category discrimination 
        parameters.  The first column should contain the parameters for the 'do not know' category. 
        For any items with fewer categories than the maximum, the remaining cells in the (k+1)/3 
        row should be \code{NA}.}

        \item{Category Difficulty Parameters:}{The second set of (k+1)/3 columns are for category 
        difficulty parameters.  The first column of this subset of columns should contain the 
        category difficulties for the 'do not know'  category.  Similar to the discrimination 
        columns, items with fewer response categories than the maximum should have \code{NA}s 
        in the remaining cells of the (k+1)/3 subset of columns.} 

        \item{Guessing  Probabilities:}{The remaining ((k+1)/3)-1 columns are for the guessing
        probabilities.  As with the discrimination and category difficulties, the remaining cells 
        for items with fewer categories than the maximum should be filled with \code{NA}s.}
      } 

      For example, for one four category and one five category item, the first column includes 
      the discrimination parameters for the 'do not know' categories.  Columns 2-5 contain the
      discrimination parameters for the other categories.  The fifth column for the four 
      response item should be \code{NA}.  The sixth column includes the category difficulties
      for the 'do not know' category.  Columns 7-10 contain the difficulty parameters for the
      other categories.  The tenth column for the four response item should be \code{NA}.  The
      remaining four columns include the guessing  probabilities, although the last column
      for the four response item should be \code{NA}.}

    \item{x = "data.frame", cat = "numeric"}{See the method for x = "matrix" }

    \item{x = "list", cat = "numeric"}{This method is for a list with three elements.  The first 
      list element is for discimination parametes, the second is for category difficulty
      parameters, and the third is for alternative response probability parameters.
      \describe{
        \item{Discrimination Parameters:}{The first list element should conatain an n x k matrix of 
        category discrimination values for n items and k equal to the maximum number of response 
        categories (including the 'do not know' category) across all items.  The first column of this 
        matrix should include the parameters for the 'do not know' category.  For any items with fewer 
        categories than the maximum, the remaining cells in the row should be \code{NA}.}
    
       \item{Category Difficulty Parameters:}{The second list element is an n x k matrix of category 
        difficulty parameters.  The first element of this matrix should include the category 
        difficulties for the 'do not know' category.  Similar to the discrimination list element, 
        items with fewer response categories than the maximum should have \code{NA}s in the remaining 
        cells.} 

      \item{Guessing  Probabilities:}{The third list element is an n x (k-1) matrix of guessing 
        probabilities (i.e. there is no column for the 'do not know' category). As with the 
        discrimination and category difficulties, the remaining cells for items with fewer categories 
        than the maximum should be filled with \code{NA}s.}
      } 

      For example, for one four category and one five category item, the first four columns for 
      the four response item in the first list element would include the discrimination parameters. 
      The fifth column for this item would be \code{NA}.  The first four columns for the four 
      response item in the second list element would include the category difficulty parameters. 
      The fifth column for this item would be \code{NA}.  The first three columns for the four
      response item in the third list element would include the guessing probabilities.  The fourth
      column for this item would be \code{NA}.}

    \item{x = "irt.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the mcm items in an object of class \code{"irt.pars"}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the mcm items only.  If \code{x} 
      contains parameters for multiple groups, a list of \code{"irt.prob"} objects will be returned.}

    \item{x = "sep.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the mcm items in an object of class \code{sep.pars}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the mcm items only.}
  }
}

\value{
  Returns an object of class \code{\linkS4class{irt.prob}}
}

\references{ 
  Thissen, D., & Steinberg, L. (1984). 
  A response model for multiple choice items.
  \emph{Psychometrika}, \emph{49}(4), 501-519.

  Thissen, D., & Steinberg, L. (1996)
  A response model for multiple choice items.
  In W.J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory}. New York: Springer-Verlag
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com}} }

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
## Item parameters from Thissen & Steinberg (1984, p. 510)
## Items R,S,T,U for the whole test
a <- matrix(c(-1.7, -1, 1.1, .3, 1.9, -2.1, -.6, 1.2, 2.3, -.8, -1.3, -.9, -.2, 1.9, .5, 
  -1.9, -.5, 0, -.6, 1.9),4,5,byrow=TRUE)
c <- matrix(c(.3, -2.3, 2.4, -2.5, 2.1, 2.1, .05, -3, -.6, 1, -.9, -2.5, -.1, 1.8, 1.6, 
  -.1, -2, .5, .8, .8),4,5,byrow=TRUE)
d <- matrix(c(.25, .25, .25, .25, .2, .2, .4, .2, .2, .2, .4, .2, .25, .25, .25, .25),
  4,4,byrow=TRUE)
pars <- cbind(a,c,d)
x <- mcm(pars, rep(5,4))
plot(x,item.names=paste("Item",c("R","S","T","U")), auto.key=list(space="right"))

## Item parameters from Thissen & Steinberg (1984, p. 511)
## Items W,X,Y,Z for the 
pars <- vector("list",3)
pars[[1]] <- matrix(c(-2.3, -.2, 2, .9, -.3, -.8, .6, -.5, 1.1, -.4, -.5, -.2, 2, -1.2, 0, 
  -1.5, -.7, -.2, .1, 2.3),4,5,byrow=TRUE)
pars[[2]] <- matrix(c(.5, .7, -.5, -1.9, 1.1, 1.6, -2.8, 1.5, 0, -.3, -.3, .7, -1, .7, 0, 
  .4, .4, -.5, .5, -.8),4,5,byrow=TRUE)
pars[[3]] <- matrix(c(.2, .4, .2, .2, .2, .2, .4, .2, .2, .4, .2, .2, .2, .2, .2, .4),
  4,4,byrow=TRUE)
x <- mcm(pars, rep(5,4))
plot(x,item.names=paste("Item",c("W","X","Y","Z")), auto.key=list(space="right"))
}

\keyword{methods}
\keyword{distribution}
