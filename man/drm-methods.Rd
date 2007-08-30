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
responses) for one or more items for a given set of theta values using the 1PL, 2PL, 
or 3PL model, depending on the included item parameters.
}

\usage{
drm(x, theta = seq(-4,4,0.05), D = 1.7, incorrect = FALSE, print.mod = FALSE, ...)

## S4 method for signature 'numeric'
drm(x, theta, D, incorrect, print.mod, ...)

## S4 method for signature 'matrix'
drm(x, theta, D, incorrect, print.mod, ...)

## S4 method for signature 'data.frame'
drm(x, theta, D, incorrect, print.mod, ...)

## S4 method for signature 'list'
drm(x, theta, D, incorrect, print.mod, ...)

## S4 method for signature 'irt.pars'
drm(x, theta, D, incorrect, print.mod, ...)

## S4 method for signature 'sep.pars'
drm(x, theta, D, incorrect, print.mod, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters}

  \item{theta}{vector of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.05 }

  \item{D}{scaling constant (for the normal metric \code{D} = 1.7, for the logistic metric \code{D} = 1) }

  \item{incorrect}{if \code{TRUE}, compute probabilities for the incorrect response category) }

  \item{print.mod}{if \code{TRUE}, print the model that was used (i.e. 1PL, 2PL, or 3PL) }

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
  \describe{

    \item{x = "numeric"}{This method should only be used for the Rasch model where \code{x} is 
      a vector of difficulty parameters.  Under this method the discrimination values will 
      default to 1 for all items. }

    \item{x = "matrix"}{This method allows one to specify an n x k matrix for n items and k
      equal to 1, 2, or 3.  The columns generally correspond to discrimination, difficulty, 
      and lower asymptote parameters, although this may not be the case depending on the 
      supplied parameters.
        \describe{
          \item{1PL:}{For the 1PL model with discriminations equal to 1 (Rasch Model), an n x 1 
          matrix of item difficulties can be supplied.  An n x 2 matrix can also be used with
          all the values in the first column equal to 1 and difficulty parameters in the second
          column.  For discrimination values other than 1, \code{x} should contain at least two 
          columns, the first for item discriminations (identical for all items) and the second 
          for item difficulties.  The lower asymptote defaults to zero for all items; however, a 
          third column of zeros can be included.}

          \item{2PL:}{\code{x} should include at lease two columns, the first for item 
          discriminations and the second for item difficulties.  The lower asymptote defaults to 
          zero for all items; however, a third column of zeros can be included.}

          \item{3PL:}{\code{x} should include three columns, the first for item discriminations, 
          the second for item difficulties, and the third for lower asymptote values.}
	}
    }

    \item{x = "data.frame"}{See the method for x = "matrix" }

    \item{x = "list"}{This method can include a list with 1, 2, or 3 elements.  In general, these 
      elements correspond to discrimination, difficulty, and lower asymptote parameters, although 
      this may not be the case depending on the supplied parameters.  Each list element can be a 
      vector of length n or an n x 1 matrix for n items. 
        \describe{
          \item{1PL:}{For the 1PL model with discriminations equal to 1 (Rasch Model), one element
          with item difficulties can be supplied.  Alternatively, two elements can be used with
          the first list element containing a matrix/vector of ones and difficulty parameters in 
          the second list element.  For discrimination values other than 1, \code{x} should contain 
          at  least two list elements, the first for item discriminations (identical for all items) 
          and the second for item difficulties.  The lower asymptote defaults to zero for all items; 
          however, a third element with a vector/matrix of zeros can be included.}

          \item{2PL:}{\code{x} should include at lease two list elements, the first for item 
          discriminations and the second for item difficulties.  The lower asymptote defaults to 
          zero for all items; however, a third element with a vector/matrix of zeros can be included.}

          \item{3PL:}{\code{x} should include three list elements, the first for item discriminations, 
          the second for item difficulties, and the third for lower asymptote values.}
	}
    }

    \item{x = "irt.pars"}{This method can be used to compute probabilities for the dichotomous
      items in an object of class \code{"irt.pars"}.  If \code{x} contains polytomous items, a
      warning will be displayed stating that probabilities will be computed for the 
      dichotomous items only.  If \code{x} contains parameters for multiple groups, a list of 
      \code{"irt.prob"} objects will be returned.}

    \item{x = "sep.pars"}{This method can be used to compute probabilities for the dichotomous
      items in an object of class \code{sep.pars}.  If \code{x} contains polytomous items, a
      warning will be displayed stating that probabilities will be computed for the 
      dichotomous items only. }
  }
}

\value{
  Returns an object of class \code{\linkS4class{irt.prob}}
}

\references{ 
  Birnbaum, A. (1968)
  Some latent trait models and their use in inferring an examinee's ability
  In F. M. Lord & M. R. Novick (Eds.)
  \emph{Statistical Theories of Mental Test Scores.}
  Reading, MA: Addison-Wesley

  Rasch, G. (1960)
  \emph{Probabilistic Models for Some Intelligence and Attainment Tests}
  Copenhagen, Denmark: Danish Institute for Educational Research
}

\author{ Jonathan Weeks \email{weeksjp@gmail.com} }

\note{ 
The indices 0 and 1 are used in the column labels for incorrect and correct responses 
respectively (e.g. the label d1.0 indicates the incorrect response column for item 1, d1.1 
indicates the correct response column for the same item).  If \code{incorrect} = \code{FALSE}, 
all column labels will end with .1
}

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
## 1PL
# A single item with a difficulty at 1 logit
x <- drm(1) 
plot(x)

# Three items with a difficulties at -2.2, -1.7, and 0.54 logits respectively
# Assume a constant discrimination of 1 for each item
b <- c(-2.2,-.7,1.54)
x <- drm(b) 
plot(x)

# Five items with a constant discrimination equal to 1.2
# Do not use the scaling constant
a <- rep(1.2,5)
b <- c(-2,-1,0,1,2)
pars <- list(a,b)
x <- drm(pars,D=1)
plot(x,combine=5,item.names="Items 1-5")


## 2PL
# Compute probabilities for five items at a theta value of 1.5 logits
# Print the model
a <- c(1,.8,.6,1.2,.9)
b <- c(-1.7,2.2,.4,-1,-.4)
pars <- cbind(a,b)
drm(pars,theta=1.5,print.mod=TRUE)@prob

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
}

\keyword{methods}
\keyword{distribution}
