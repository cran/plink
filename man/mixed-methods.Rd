\name{mixed-methods}
\docType{methods}
\alias{mixed}
\alias{mixed-methods}
\alias{mixed,numeric,numeric-method}
\alias{mixed,matrix,numeric-method}
\alias{mixed,data.frame,numeric-method}
\alias{mixed,list,numeric-method}
\alias{mixed,irt.pars,ANY-method}
\alias{mixed,sep.pars,ANY-method}
\title{Mixed-Format Response Probabilities}
\description{
This function computes the probability of responding in a specific category for one or more 
items for a given set of theta values when the items are from a mixed-format test.
}

\usage{
mixed(x, cat, poly.mod, theta = seq(-4,4,0.05), ...)

## S4 method for signature 'numeric', 'numeric'
mixed(x, cat, poly.mod, theta, ...)

## S4 method for signature 'matrix', 'numeric'
mixed(x, cat, poly.mod, theta, ...)

## S4 method for signature 'data.frame', 'numeric'
mixed(x, cat, poly.mod, theta, ...)

## S4 method for signature 'list', 'numeric'
mixed(x, cat, poly.mod, theta, ...)

## S4 method for signature 'irt.pars', 'ANY'
mixed(x, cat, poly.mod, theta, ...)

## S4 method for signature 'sep.pars', 'ANY'
mixed(x, cat, poly.mod, theta, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters}

  \item{cat}{vector identifying the number of response categories for each item.  If 
    multiple-choice model items are included, \code{cat} for these items should equal 
    the number of response categories plus one (the additional category is for 
    'do not know')}

  \item{poly.mod}{object of class \code{\link[=poly.mod-class]{poly.mod}} identifying 
    the items associated with each IRT model}

  \item{theta}{vector of theta values for which probabilities will be computed.  If 
    \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
    with an increment of 0.05 }

  \item{...}{further arguments passed to or from other methods }
}

\details{
The item parameters supplied to this method can be associated with a single IRT model or
multiple models.  When the parameters are tied to only one model, the format of \code{x}
should follow the conventions in \code{\link{drm}} for dichotomous response models (i.e.
1PL, 2PL, 3PL), \code{\link{gpcm}} for the partial credit model and generalized partial
credit model, \code{\link{grm}} for the graded response model, \code{\link{mcm}} for the
multiple-choice model, and \code{\link{nrm}} for the nominal response model.  When the 
parameters are associated with two or more models, the parameters need to be combined.
See \code{\link{as.irt.pars}} or \code{\link{sep.pars}} for a complete description of 
how the parameters from different models can be combined.

Additional arguments for the above models can be passed to this method as well.
}

\value{
  Returns an object of class \code{\linkS4class{irt.prob}}
}

\author{ Jonathan Weeks \email{weeksjp@gmail.com}} }

\seealso{ 
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
# Compute probabilites for three dichotomous (3PL) items and two polytomous (gpcm) items
# without a location parameter
dichot <- matrix(c(1.2, .8, .9, 2.3, -1.1, -.2, .24, .19, .13),3,3)
poly <- matrix(c(.64, -1.8, -.73, .45, NA, .88, .06, 1.4, 1.9, 2.6),2,5,byrow=TRUE)
pars <- rbind(cbind(dichot,matrix(NA,3,2)),poly)
cat <- c(2,2,2,4,5)
pm <- as.poly.mod(5, c("drm","gpcm"), list(1:3,4:5))
x <- mixed(pars, cat, pm) 
plot(x)


# Compute probabilities for three dichotomous (3PL) items, four polytomous items,
# two gpcm items and two nrm items. Include a location parameter for the 
# gpcm items
a <- matrix(c(
  1.2, rep(NA,4),
  .8, rep(NA,4),
  .9, rep(NA,4),
  .64, rep(NA,4),
  .88, rep(NA,4),
  .905, .522, -.469, -.959, NA, 
  .828, .375, -.357, -.079, -.817),7,5,byrow=TRUE)
b <- matrix(c(
  2.3, rep(NA,4),
  -1.1, rep(NA,4),
  -.2, rep(NA,4),
  -.69, -1.11, -.04, 1.14, NA,
  1.49, -1.43, -.09, .41, 1.11,
  .126, -.206, -.257, .336, NA, 
  .565, .865, -1.186, -1.199, .993),7,5,byrow=TRUE)
c <- c(.14, .19, .26, rep(NA,4))
pars <- list(a,b,c)
cat <- c(2,2,2,4,5,4,5)
pm <- as.poly.mod(7, c("drm","gpcm","nrm"), list(1:3,4:5,6:7))
x <- mixed(pars, cat, pm, location=TRUE)
plot(x)
}


\keyword{methods}
\keyword{distribution}
