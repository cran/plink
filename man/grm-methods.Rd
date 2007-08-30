\name{grm-methods}
\docType{methods}
\alias{grm}
\alias{grm-methods}
\alias{grm,matrix,numeric-method}
\alias{grm,data.frame,numeric-method}
\alias{grm,list,numeric-method}
\alias{grm,irt.pars,ANY-method}
\alias{grm,sep.pars,ANY-method}
\title{Graded Response Model Probabilities}
\description{
This function computes the cumulative probability of responding within or above a certain 
category or the probability of responding in a specific category for one or more items for 
a given set of theta values using the graded response model.
}

\usage{
grm(x, cat, theta = seq(-4,4,0.05), catprob = FALSE, D = 1.7, location = FALSE, ...)

## S4 method for signature 'matrix', 'numeric'
grm(x, cat, theta, catprob, D, location, ...)

## S4 method for signature 'data.frame', 'numeric'
grm(x, cat, theta, catprob, D, location, ...)

## S4 method for signature 'list', 'numeric'
grm(x, cat, theta, catprob, D, location, ...)

## S4 method for signature 'irt.pars', 'ANY'
grm(x, cat, theta, catprob, D, location, ...)

## S4 method for signature 'sep.pars', 'ANY'
grm(x, cat, theta, catprob, D, location, ...)
}

\arguments{
  \item{x}{ans \code{R} object containing item parameters.}

  \item{cat}{vector identifying the number of response categories (not the number of step
    parameters) for each item.}

  \item{theta}{vector of theta values for which probabilities will be computed. 
    If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
    with an increment of 0.05 }

  \item{catprob}{if \code{TRUE}, compute category probabilities instead of cumulative
    probabilities }

  \item{D}{scaling constant (for the normal metric \code{D} = 1.7, for the logistic metric \code{D} = 1) }

  \item{location}{if \code{TRUE}, the step parameters are deviations from a difficulty
    parameter }

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
  \describe{

    \item{x = "matrix", cat = "numeric"}{ This method allows one to specify an n x k 
      matrix for n items.  The number of columns can vary depending on whether a location 
      parameter is included.  The first column is for item discriminations and the remaining 
      columns for threshold/threshold deviation parameters.
	
      Threshold parameters can be characterized in two ways: the actual thresholds or 
      deviations from an overall item difficulty (location).  In the deviation scenario the 
      \code{location} argument should be \code{TRUE} and the location parameters should be 
      in the second column.  The columns for the threshold/threshold deviation parameters 
      will always follow the discrimination and, if applicable, location columns.  The number 
      of threshold/threshold deviation parameters can vary for each item.  In these instances, 
      all cells with missing values should be filled with \code{NA}s.  For example, for a model 
      with no location parameter, if one item has five categories (four threshold parameters) 
      and another item has three categories (two threshol parameters), there should be five 
      columns.  The first column includes the discrimination parameters and columns 2-5 include 
      the threshold parameters.  The values in the last two columns for the item with three 
      categories should be \code{NA}.}

    \item{x = "data.frame", cat = "numeric"}{See the method for x = "matrix" }

    \item{x = "list", cat = "numeric"}{This method is for a list with two elements.  The first 
      list element is for item discriminations and the second for threshold/threshold deviation 
      parameters.  The discrimination values should be a vector of length n or an n x 1 
      matrix for n items.  The threshold/threshold deviation parameters should be formatted 
      as an n x k matrix.  If the thresholds are deviations from a location parameter, the 
      argument \code{location} should be \code{TRUE} and the location parameters should be in 
      the first column.  The number of threshold/threshold deviation parameters can vary for 
      each item.  In these instances, all cells with missing values should be filled with 
      \code{NA}s (See the example in the method for x = "matrix"). }

    \item{x = "irt.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the grm items in an object of class \code{"irt.pars"}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the grm items only. If \code{x} 
      contains parameters for multiple groups, a list of \code{"irt.prob"} objects will be returned.}

    \item{x = "sep.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the grm items in an object of class \code{sep.pars}. If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the grm items only.}
  }
}

\value{
  Returns an object of class \code{\linkS4class{irt.prob}}
}

\references{ 
  Embretson, S. E., & Reise, S. P. (2000)
  \emph{Item Response Theory for Psychologists.}
  Mahwah, New Jersey: Lawrence Erlbaum Associates

  Kolen, M. J. & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking (2nd ed.).}
  New York: Springer

  Samejima, F. (1969)
  Estimation of latent ability using a response pattern of graded scores.
  \emph{Psychometrika Monograph}, No. 17

  Samejima, F. (1996)
  The graded response model.
  In W. J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory} (pp. 85-100). 
  New York: Springer-Verlag
}

\author{ Jonathan Weeks \email{weeksjp@gmail.com} }

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
## Item parameters from Embretson & Reise (2000, p. 101) items 9-11
## Cumulative probabilities 
# Includes a location parameter
a <- c(2.09,1.18,1.69)
b <- matrix(c(-1.07,-1.03,.39,.86,1.78,1.85,-.87,-.39,.31,NA,-1.4,-.42,NA,NA,-1.74),3,5) 
pars <- cbind(a,b)
rownames(pars) <- paste("Item",9:11,sep=" ")
colnames(pars) <- c("a","b",paste("c",1:4,sep=""))
cat <- c(3,4,5)
x <- grm(pars,cat,location=TRUE)
plot(x)

## Category probabilities
# Does not include a location parameter
a <- c(2.09,1.18,1.69)
b <- matrix(c(-1.93,-2.81,-1.46,-.2,-.64,.08,NA,.37,.81,NA,NA,2.13),3,4) 
pars <- cbind(a,b)
rownames(pars) <- paste("Item",9:11,sep=" ")
colnames(pars) <- c("a",paste("b",1:4,sep=""))
cat <- c(3,4,5)
x <- grm(pars,cat,catprob=TRUE)
plot(x)

## Item parameters from Kolen & Brennan (2004, p. 210)
## Cumulative probabilities for a single item
# Use theta values from -3 to 3 with an increment of 0.5
pars <- t(c(1.2,-.5,.6,1.1,1.3))
x <- grm(pars,5,theta=seq(-3,3,.05))
plot(x,item.lab=FALSE)

## Category probabilities for a single item
pars <- t(c(1.2,-.5,.6,1.1,1.3))
x <- grm(pars,5,seq(-3,3,.05),catprob=TRUE)
plot(x,item.lab=FALSE)

## Create sep.pars object then compute category probabilities
a <- c(2.09,1.18,1.69)
b <- matrix(c(-1.93,-2.81,-1.46,-.2,-.64,.08,NA,.37,.81,NA,NA,2.13),3,4) 
pars <- cbind(a,b)
cat <- c(3,4,5)
pm <- as.poly.mod(3,"grm")
out <- sep.pars(pars,cat,pm)
x <- grm(out,catprob=TRUE)
plot(x)
}

\keyword{methods}
\keyword{distribution}
