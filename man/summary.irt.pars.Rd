\name{summary.irt.pars}
\alias{summary.irt.pars}
\alias{summary.sep.pars}
\alias{summary.list}
\alias{summary.link}
\title{Summarize Item Parameters/plink Output }
\description{
This function summarizes the item parameters for a given object.  If the object is
of class \code{\linkS4class{link}}, the function prints the linking constants
}

\usage{
\method{summary}{irt.pars}(object, ..., descrip = FALSE)

\method{summary}{sep.pars}(object, ..., descrip = FALSE)

\method{summary}{list}(object, ..., descrip = FALSE)

\method{summary}{link}(object, ..., descrip = FALSE)
}

\arguments{

  \item{object}{an \code{R} object. See details below}

  \item{\dots}{further arguments passed to or from other methods}

  \item{descrip}{if \code{TRUE},print descriptive statistics for the item parameters in 
    \code{object}. }
}

\details{
The method for objects of class \code{irt.pars} summarizes the parameters for a single 
group or multiple groups.  If multiple groups are included in the object, parameters will 
be summarized for each group separately.

The method for objects of class \code{sep.pars} summarizes the parameters for a single
\code{sep.pars} object.

Objects of class \code{list} can include any combination of \code{irt.pars} or
\code{sep.pars} objects.  The parameters will be summarized for each group separately.

The method for objects of class \code{link} prints the linking constants for n-1 groups.
If \code{descrip} = \code{TRUE}, the method summarizes the parameters for all groups 
separately.  The labels are specified in the following manner "group1/group2", meaning
the group1 parameters were transformed to the group2 test using the associated constants.
The base group is indicated by an asterisk.
}

\references{
  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{
\code{\link{as.irt.pars}}, \code{\linkS4class{irt.pars}}, \code{\link[=sep.pars-methods]{sep.pars}},
\code{\linkS4class{sep.pars}}, \code{\linkS4class{link}}
}

\keyword{ utilities }
