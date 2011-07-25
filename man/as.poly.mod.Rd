\name{as.poly.mod}
\alias{as.poly.mod}
\title{poly.mod objects}
\description{
This function attempts to turn the given values into a \code{\linkS4class{poly.mod}} object that
associates each item with a specific unidimensional or multidimensional item response model.
}
\usage{
as.poly.mod(n, model = "drm", items = NULL)
}

\arguments{
  \item{n}{total number of items}

  \item{model}{character vector identifying the IRT models used to estimate
    the item parameters.  The only acceptable models are 
    \code{drm, gpcm, grm, mcm, and nrm}.  See below for an explanation of
    the codes. }

  \item{items}{list identifying the item numbers from a set of parameters
	that correspond to the given model in \code{model}.}
}

\details{
When creating a \code{poly.mod} object, there is no difference in the specification for
unidimensional versus multidimensional item response models.  If all the items are dichotomous, 
it is only necessary to specify a value for \code{n}.  If all the items correspond to a 
single model (other than \code{drm}), only \code{n} and \code{model} need to be specified.

The IRT models associated with the codes:
  \describe{
    \item{\code{drm}:}{dichotomous response models (includes the 1PL, 2PL, 3PL, M1PL, 
      M2PL, and M3PL)}
    \item{\code{gpcm}:}{partial credit model, generalized partial credit model,
      multidimensional partial credit model, and multidimensional generalized partial credit model}
    \item{\code{grm}:}{graded response model and multidimensional graded response model}
    \item{\code{mcm}:}{multiple-choice model and multidimensional multiple-choice model}
    \item{\code{nrm}:}{nominal response model and multidimensional nominal response model}
  }
}

\value{
Returns an object of class \code{\linkS4class{poly.mod}}
}

\references{
  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{ \code{\linkS4class{poly.mod}} }

\examples{
# Ten dichotomous items
as.poly.mod(10)

# The first ten items in the set of associated (not present here) item 
# parameters are dichotomous and the last five were estimated using the 
# generalized partial credit model
as.poly.mod(15, c("drm", "gpcm"), list(1:10,11:15) )

# Ten multidimensional graded response model items
# Note: This same specification would be used for a unidimensional
# graded response model
as.poly.mod(10, "grm")
}
\keyword{ misc }
