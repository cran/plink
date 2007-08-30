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
items for a given set of theta values using the nominal response model.
}

\usage{
nrm(x, cat, theta = seq(-4,4,0.05), ...)

## S4 method for signature 'matrix', 'numeric'
nrm(x, cat, theta, ...)

## S4 method for signature 'data.frame', 'numeric'
nrm(x, cat, theta, ...)

## S4 method for signature 'list', 'numeric'
nrm(x, cat, theta, ...)

## S4 method for signature 'irt.pars', 'ANY'
nrm(x, cat, theta, ...)

## S4 method for signature 'sep.pars', 'ANY'
nrm(x, cat, theta, ...)
}

\arguments{
  \item{x}{Object containing item parameters. See below for more details. }

  \item{cat}{vector identifying the number of response categories for each item } 

  \item{theta}{vector of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.05 }

  \item{...}{further arguments passed to or from other methods }
}


\section{Methods}{
  \describe{

    \item{x = "matrix", cat = "numeric"}{This method allows one to specify an n x k matrix
	for n items and k equal two times the maximum number of response categories across
	items (e.g. if some items have four response categories and others have five 
	categories, k will equal 10). The first k/2 columns are for category discrimination 
	parameters. For any items with fewer categories than the maximum, the remaining
	cells in the k/2 row should be \code{NA}. The remaining columns are for the
	category difficulty parameters. Similarly, items with fewer response categories
	than the maximum should have \code{NA}s in the remaining cells. For example, for 
	one four category and one five category item, the first four columns for the four
	response item would include the discrimination parameters. The fifth column for
	this item would be \code{NA}. The next four columns would include the category
	difficulty values, and the last column would be \code{NA}. }

    \item{x = "data.frame", cat = "numeric"}{See the method for x = "matrix" }

    \item{x = "list", cat = "numeric"}{This method is for a list with two elements. The 
	first element is and n x k matrix of category discrimination values for n items 
	and k equal to the maximum number of response categories across items. The second 
	list element is an n x k matrix of category difficulty parameters. For either 
	element, for items with fewer categories than the maximum, the remaining cells in
	the rows should be \code{NA}. For example, for one four category and one five 
	category item, the first four columns for the four response item in the first list
	element would include the discrimination parameters. The fifth column for
	this item would be \code{NA}.The first four columns for the four response item in 
	the second list element would include the category difficulty parameters. The fifth 
	column for this item would be \code{NA}. }

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

  Kolen, M. J., & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking.}
  New York: Springer
}

\author{ Jonathan Weeks \email{weeksjp@gmail.com}} }

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
## Item parameters from Bock (1972, p. 46,47)
a <- matrix(c(.905, .522, -.469, -.959, NA, .828, .375, -.357, -.079, -.817),
  2,5,byrow=TRUE)
c <- matrix(c(.126, -.206, -.257, .336, NA, .565, .865, -1.186, -1.199, .993),
  2,5,byrow=TRUE)
pars <- cbind(a,c)
x <- nrm(pars, c(4,5))
plot(x,auto.key=list(space="right"))


}

\keyword{methods}
\keyword{distribution}
