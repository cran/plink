\name{as.poly.mod}
\alias{as.poly.mod}
\title{poly.mod objects}
\description{
This function attempts to turn the given values into a poly.mod object
}
\usage{
as.poly.mod(n, model = "drm", items = NULL)
}

\arguments{
  \item{n}{total number of items}

  \item{model}{character vector identifying the IRT models used to estimate
    the item parameters.  The only acceptable models are 
    \code{drm, gpcm, grm, mcm, and nrm}.  See below for an explantion of
    the codes. }

  \item{items}{list identifying the item numbers from a set of parameters
	that correspond to the given model in \code{model}.}
}

\details{
If all the items are dichotomous, it is only necessary to specify a value for \code{n}.
If all the items correspond to a single model (other than \code{drm}), only \code{n}
and \code{model} need to be specified.

The IRT models associted with the codes:
  \describe{
    \item{\code{drm}:}{dichotomous reponse models (includes 1PL, 2PL, and 3PL)}
    \item{\code{gpcm}:}{generalized partial credit model (includes the
      partial credit model)}
    \item{\code{grm}:}{graded response model}
    \item{\code{mcm}:}{multiple-choice model}
    \item{\code{nrm}:}{nominal response model}
  }
}

\value{
Returns an object of class \code{\linkS4class{poly.mod}}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{ \code{\linkS4class{poly.mod}} }

\examples{
# Ten dichotomous items
as.poly.mod(10)

# The first ten items in the set of associated (not present here) item parameters
# are dichotomous and the last five were estimated using the generalized
# partial credit model
as.poly.mod(15, c("drm", "gpcm"), list(1:10,11:15) )
}
\keyword{ misc }
