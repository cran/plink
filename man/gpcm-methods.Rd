\name{gpcm-methods}
\docType{methods}
\alias{gpcm}
\alias{gpcm-methods}
\alias{gpcm,matrix,numeric-method}
\alias{gpcm,data.frame,numeric-method}
\alias{gpcm,list,numeric-method}
\alias{gpcm,irt.pars,ANY-method}
\alias{gpcm,sep.pars,ANY-method}
\title{Generalized Partial Credit Model Response Probabilities}
\description{
This function computes the probability of responding in a specific category for one or 
more items for a given set of theta values using the partial credit model, generalized 
partial credit model, or multidimensional extension of these models, depending 
on the included item parameters and the specified number of dimensions.
}

\usage{
gpcm(x, cat, theta, dimensions = 1, D = 1, location = FALSE, 
  print.mod = FALSE, items, ...)

## S4 method for signature 'matrix', 'numeric'
gpcm(x, cat, theta, dimensions, D, location, print.mod, items, ...)

## S4 method for signature 'data.frame', 'numeric'
gpcm(x, cat, theta, dimensions, D, location, print.mod, items, ...)

## S4 method for signature 'list', 'numeric'
gpcm(x, cat, theta, dimensions, D, location, print.mod, items, ...)

## S4 method for signature 'irt.pars', 'ANY'
gpcm(x, cat, theta, dimensions, D, location, print.mod, items, ...)

## S4 method for signature 'sep.pars', 'ANY'
gpcm(x, cat, theta, dimensions, D, location, print.mod, items, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters }

  \item{cat}{vector identifying the number of response categories (not the number of step
  parameters) for each item. }

  \item{theta}{vector, matrix, or list of theta values for which probabilities will be computed. 
  If \code{theta} is not specified, an equal interval range of values from -4 to 4 is used 
  with an increment of 0.5. See details below for more information. }

  \item{dimensions}{number of modeled dimensions }

  \item{D}{scaling constant. The default value assumes that the parameters are already in the 
  desired metric. If the parameters are in the logistic metric, they can be transformed to a 
  normal metric by setting D = 1.7}

  \item{location}{if \code{TRUE}, the step parameters are deviations from a difficulty
    parameter }

  \item{print.mod}{if \code{TRUE}, print the model that was used (i.e., Partial Credit
    Model, Generalized Partial Credit Model, Multidimensional Partial Credit
    Model or Multidimensional Generalized Partial Credit Model) }

  \item{items}{numeric vector identifying the items for which probabilities should be computed}

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
In the following description, references to the partial credit model and generalized
partial credit model should be thought of as encompassing both the unidimensional and 
multidimensional models.
  \describe{
    \item{x = "matrix", cat = "numeric"}{This method allows one to specify an n x k 
      matrix for n items.  The number of columns can vary depending on the model 
      (partial credit or generalized partial credit model), number of dimensions, and 
      whether a location parameter is included.  Generally, the first m columns, for m
      dimensions, are for item slopes and the remaining columns are for step parameters.
        \describe{
          \item{Slope Parameters:}{The partial credit model is typically specified with 
          all slopes equal to 1.  For this model it is unnecessary (although optional) to 
          include ones in the first m columns.  For slope values other than one (equal for 
          all items) or for the generalized partial credit model, slope parameters should 
          be included in the first m columns.}
	
          \item{Step/Step Deviation Parameters:}{Step parameters can be characterized in
          two ways: as the actual steps or deviations from an overall item difficulty
          (location).  In the deviation scenario the \code{location} argument should be 
          equal to \code{TRUE}.  If column(s) are included for the slope parameters, 
          the location parameters should be in the m+1 column; otherwise, they should be 
          in the first column.  The columns for the step/step deviation parameters will always
          follow the slope and/or location columns (or they may potentially 
          start in the first column for the partial credit model with no location parameter).

          The number of step/step deviation parameters can vary for each item.  In these
          instances, all cells with missing values should be filled with \code{NA}s.  For
          example, for a unidimensional generalized partial credit model with no location 
          parameter, if one item has five categories (four step parameters) and another item 
          has three categories (two step parameters), there should be five columns.  The first 
          column includes the slope parameters and columns 2-5 include the step parameters.
          The values in the last two columns for the item with three categories should be \code{NA}.}
	}
    }
 
    \item{x = "data.frame", cat = "numeric"}{See the method for x = "matrix" }

    \item{x = "list", cat = "numeric"}{This method can include a list with one or two
        elements. Generally, the first element is for item slopes and the second is
        for step/step deviation parameters. 
        \describe{
          \item{Slope Parameters:}{For the partial credit model with all slopes equal to 1 it is 
          unnecessary (although optional) to include a list element for the item slopes.  If no 
          slope values are included, the first element would contain the step/deviation step 
          parameters.  For slopes other than 1 (equal for all items) or for the generalized partial 
          credit model, slope values should be included in the first list element.  For the 
          unidimensional case, these values should be a vector of length n or an n x 1 matrix for 
          n items. For the multidimensional case, an n x m matrix of values for m dimensions 
          should be supplied}

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
      will be returned. The argument \code{dimensions} does not need to be included for this method.}

    \item{x = "sep.pars", cat = "ANY"}{This method can be used to compute probabilities 
      for the gpcm items in an object of class \code{sep.pars}.  If \code{x} contains 
      dichotomous items or items associated with another polytomous model, a warning will be 
      displayed stating that probabilities will be computed for the gpcm items only. The argument 
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
  Adams, R. J., Wilson, M., & Wang, W. (1997).
  The multidimensional random coefficients multinomial logit model.
  \emph{Applied Psychological Measurement}, \emph{21}(1), 1-23.

  Embretson, S. E., & Reise, S. P. (2000).
  \emph{Item Response Theory for Psychologists.}
  Mahwah, New Jersey: Lawrence Erlbaum Associates.

  Masters, G. N. (1982).
  A rasch model for partial credit scoring.
  \emph{Psychometrika}, \emph{47}(2), 149-174.

  Masters, G. N. & Wright, B. D. (1996)
  The partial credit model.
  In W. J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory} (pp. 101-121). 
  New York: Springer-Verlag.

  Muraki, E. (1992)
  A generalized partial credit model: Application of an EM algorithm.
  \emph{Applied Psychological Measurement}, \emph{16}(2), 159-176.

  Muraki, E. (1996)
  A generalized partial credit model.
  In W. J. van der Linden & Hambleton, R. K. (Eds.) 
  \emph{Handbook of Modern Item Response Theory} (pp. 153-164). 
  New York: Springer-Verlag.

  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}

  Yao, L. (2003). 
  BMIRT: Bayesian multivariate item response theory [Computer Program]. 
  Monterey, CA: CTB/McGraw-Hill.	

  Yao, L., & Schwarz, R. D. (2006). 
  A multidimensional partial credit model with associated item and test statistics: 
  An application to mixed-format tests. 
  \emph{Applied Psychological Measurement}, \emph{30}(6), 469-492.
}

\note{ 
The determination of the model (partial credit or generalized partial credit) is based on
the number of non-NA columns for each item in \code{x} and the corresponding values in
\code{cat}.
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{ 
\code{\link[=mixed-methods]{mixed}:} compute probabilities for mixed-format items\cr
\code{\link[=plot.irt.prob]{plot}:} plot item characteristic/category curves\cr
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}:} classes
}

\examples{
###### Unidimensional Examples ######
## Partial Credit Model
## Item parameters from Embretson & Reise (2000, p. 108) item 5
b <- t(c(-2.519,-.063,.17,2.055))
x <- gpcm(b,5)
plot(x) 

## Generalized Partial Credit Model
## Item parameters from Embretson & Reise (2000, p. 112) items 5-7
a <- c(.683,1.073,.583)
b <- matrix(c(-3.513,-.041,.182,NA,-.873,.358,-.226,
  1.547,-4.493,-.004,NA,NA),3,4,byrow=TRUE)
pars <- cbind(a,b) # Does not include a location parameter
rownames(pars) <- paste("Item",5:7,sep="") 
colnames(pars) <- c("a",paste("b",1:4,sep=""))
cat <- c(4,5,3)
x <- gpcm(pars,cat,seq(-3,3,.05))
plot(x)

## Item parameters from Muraki (1996, p. 154)
a <- c(1,.5)
b <- matrix(c(.25,-1.75,1.75,.75,-1.25,1.25),2,3,byrow=TRUE)
pars <- cbind(a,b) # Include a location parameter
rownames(pars) <- paste("Item",1:2,sep="")
colnames(pars) <- c("a","b",paste("d",1:2,sep=""))
cat <- c(3,3)
x <- gpcm(pars,cat,location=TRUE,print.mod=TRUE, D=1.7)
# Plot category curves for two items
matplot(x@prob$theta,x@prob[,2:4],xlab="Theta",ylab="Probability",
  ylim=c(0,1),lty=1,type="l",col="black") 
par(new=TRUE)
matplot(x@prob$theta,x@prob[,5:7],xlab="Theta",ylab="Probability",
  ylim=c(0,1),lty=3,type="l",col="black") 


###### Multidimensional Examples ######
## Multidimensional Partial Credit Model
pars <- matrix(c(2.4207,0.245,-1.1041,NA,
2.173,-0.4576,NA,NA,
2.1103,-0.8227,.4504,NA,
3.2023,1.0251,-.7837,-1.3062),4,4,byrow=TRUE)
cat <- c(4,3,4,5)
x <- gpcm(pars,cat,dimensions=2,print.mod=TRUE)
# plot combined item category surfaces
# The screen argument adjusts the orientation of the axes
plot(x,screen=list(z=-60,x=-70))

## Multidimensional Generalized Partial Credit Model
a <- matrix(c(.873, .226, .516, .380, .613, .286 ),3,2,byrow=TRUE)
b <- matrix(c(2.255, 1.334, -.503, -2.051, -3.082,
1.917, 1.074, -.497, -1.521, -2.589,
1.624, .994, -.656, -1.978, NA),3,5,byrow=TRUE)
pars <- cbind(a,b)
cat <- c(6,6,5)
x <- gpcm(pars,cat,dimensions=2,print.mod=TRUE)

# plot combined item category surfaces
plot(x,screen=list(z=-40,x=-60), auto.key=list(space="right"))

# plot separated item category surfaces for item two
plot(x,items=2,separate=TRUE,drape=TRUE,panels=1)

# Compute response probabilities for a single three-category item with
# three dimensions. Plot the response surfaces for the first two 
# dimensions conditional on each theta value on the third dimension
pars <- matrix(c(1.1999,0.5997,0.8087,2.1730,-1.4576),1,5)
x <- gpcm(pars,3,dimensions=3,theta=-4:4)
plot(x, screen=list(z=-30,x=-60))
}

\keyword{methods}
\keyword{distribution}
