\name{gpcm-methods}
\docType{methods}
\alias{gpcm}
\alias{gpcm-methods}
\alias{gpcm,matrix,numeric-method}
\alias{gpcm,data.frame,numeric-method}
\alias{gpcm,list,numeric-method}
\alias{gpcm,irt.pars,ANY-method}
\alias{gpcm,sep.pars,ANY-method}
\title{Partial Credit/Generalized Partial Credit Model Response Probabilities}
\description{
This function computes the probability of responding in a specific category for one or 
more items for a given set of theta values using the partial credit model or the 
generalized partial credit model, depending on the included item parameters.
}

\usage{
gpcm(x, cat, theta = seq(-4,4,0.05), D = 1.7, location = FALSE, print.mod = FALSE, ...)

## S4 method for signature 'matrix', 'numeric'
gpcm(x, cat, theta, D, location, print.mod, ...)

## S4 method for signature 'data.frame', 'numeric'
gpcm(x, cat, theta, D, location, print.mod, ...)

## S4 method for signature 'list', 'numeric'
gpcm(x, cat, theta, D, location, print.mod, ...)

## S4 method for signature 'irt.pars', 'ANY'
gpcm(x, cat, theta, D, location, print.mod, ...)

## S4 method for signature 'sep.pars', 'ANY'
gpcm(x, cat, theta, D, location, print.mod, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters }

  \item{cat}{vector identifying the number of response categories (not the number of step
  parameters) for each item. }

  \item{theta}{vector of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.05 }

  \item{D}{scaling constant (for the normal metric \code{D} = 1.7, for the logistic metric \code{D} = 1) }

  \item{location}{if \code{TRUE}, the step parameters are deviations from a difficulty
    parameter }

  \item{print.mod}{if \code{TRUE}, print the model that was used (i.e. Partial Credit
    Model or Generalized Partial Credit Model) }

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
  \describe{
    \item{x = "matrix", cat = "numeric"}{This method allows one to specify an n x k 
      matrix for n items.  The number of columns can vary depending on the model 
      (partial credit or generalized partial credit model) and whether a location 
      parameter is included.  Generally, the first column is for item discriminations
      and the remaining columns are for step parameters.
        \describe{
          \item{Discrimination Parameters:}{The partial credit model is typically 
          specified with all discrimination values equal to 1.  For this model it is
          unnecessary (although optional) to include a ones in the first column.  For 
          discrimination values other than one (equal for all items) or for the 
          generalized partial credit model, discrimination parameters should be 
          included in the first column.}
	
          \item{Step/Step Deviation Parameters:}{Step parameters can be characterized in
          two ways: as the actual steps or deviations from an overall item difficulty
          (location).  In the deviation scenario the \code{location} argument should be 
          equal to \code{TRUE}.  If a column is included for the discrimination parameters, 
          the location parameters should be in the second column; otherwise, they should be 
          in the first column.  The columns for the step/step deviation parameters will always
          follow the discrimination and/or location columns (or they may potententially 
          start in the first column).

          The number of step/step deviation parameters can vary for each item.  In these
          instances, all cells with missing values should be filled with \code{NA}s.  For
          example, for the generalized partial credit model with no location parameter, if 
          one item has five categories (four step parameters) and another item has three 
          categories (two step parameters), there should be five columns.  The first column
          includes the discrimination parameters and columns 2-5 include the step parameters.
          The values in the last two columns for the item with three categories should be \code{NA}.}
	}
    }
 
    \item{x = "data.frame", cat = "numeric"}{See the method for x = "matrix" }

    \item{x = "list", cat = "numeric"}{This method can include a list with one or two
        elements. Generally, the first element is for item discriminations and the second
        for step/step deviation parameters. 
        \describe{
          \item{Discrimination Parameters:}{For the partial credit model with all 
          discrimination values equal to 1 it is unnecessary (although optional) to include
          a list element for discriminations.  If no discrimination values are included, the 
          first element would contain the step/deviation step parameters.  For discrimination 
          values other than 1 (equal for all items) or for the generalized partial credit model, 
          discrimination values should be included in the first list element.  These values should 
          be a vector with length equal to the number of items or an n x 1 matrix for n items.}

          \item{Step/Step Deviation Parameters:}{The step/step deviation parameters should be 
          formatted as an n x k matrix for n items.  If the steps are deviations from a location 
          parameter, the argument \code{location} should equal \code{TRUE} and the location 
          parameters should be in the first column.  The number of step/step deviation parameters 
          can vary for each item.  In these instances, all cells with missing values should be 
          filled with \code{NA}s (See the example in the method for x = "matrix"). }
        }
    }

    \item{x = "irt.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the gpcm items in an object of class \code{"irt.pars"}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning 
      will be displayed stating that probabilities will be computed for the gpcm items only. 
      If \code{x} contains parameters for multiple groups, a list of \code{"irt.prob"} objects 
      will be returned.}

    \item{x = "sep.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the gpcm items in an object of class \code{sep.pars}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the gpcm items only.}
  }
}

\value{
  Returns an object of class \code{\linkS4class{irt.prob}}
}

\references{ 
  Embretson, S. E., & Reise, S. P. (2000)
  \emph{Item Response Theory for Psychologists.}
  Mahwah, New Jersey: Lawrence Erlbaum Associates

  Masters, G. N. (1982)
  A rasch model for partial credit scoring.
  \emph{Psychometrika}, \emph{47}, 149-174.

  Masters, G. N. & Wright, B. D. (1996)
  The partial credit model.
  In W. J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory} (pp. 101-121). 
  New York: Springer-Verlag

  Muraki, E. (1992)
  A generalized partial credit model: Application of an EM algorithm.
  \emph{Applied Psychological Measurement}, \emph{16}, 159-176.

  Muraki, E. (1996)
  A generalized partial credit model.
  In W. J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory} (pp. 153-164). 
  New York: Springer-Verlag
}

\author{ Jonathan Weeks \email{weeksjp@gmail.com} }

\note{ 
The determination of the model (partial credit or generalized partial credit) is based on
the number of non-NA columns for each item in \code{x} and the corresponding values in
\code{cat}.
}

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
## Partial Credit Model
## Item parameters from Embretson & Reise (2000, p. 108) item 5
b <- t(c(-2.519,-.063,.17,2.055))
x <- gpcm(b,5)
plot(x) 

## Generalized Partial Credit Model
## Item parameters from Embretson & Reise (2000, p. 112) items 5-7
a <- c(.683,1.073,.583)
b <- matrix(c(-3.513,-.041,.182,NA,-.873,.358,-.226,1.547,-4.493,-.004,NA,NA),3,4,byrow=TRUE)
pars <- cbind(a,b) # Does not include a location parameter
rownames(pars) <- paste("Item",5:7,sep="") 
colnames(pars) <- c("a",paste("b",1:4,sep=""))
cat <- c(4,5,3)
x <- gpcm(pars,cat,seq(-3,3,.05))
plot(x)

## Item parameters from Muraki (1992, p. 154)
a <- c(1,.5)
b <- matrix(c(.25,-1.75,1.75,.75,-1.25,1.25),2,3,byrow=TRUE)
pars <- cbind(a,b) # Include a location parameter
rownames(pars) <- paste("Item",1:2,sep="")
colnames(pars) <- c("a","b",paste("d",1:2,sep=""))
cat <- c(3,3)
x <- gpcm(pars,cat,location=TRUE,print.mod=TRUE)
# Plot category curves for two items
matplot(seq(-4,4,.05),x@prob[,2:4],xlab="Theta",ylab="Probability",ylim=c(0,1),lty=1,
  type="l",col="black") 
par(new=TRUE)
matplot(seq(-4,4,.05),x@prob[,5:7],xlab="Theta",ylab="Probability",ylim=c(0,1),lty=3,
  type="l",col="black") 
}

\keyword{methods}
\keyword{distribution}
