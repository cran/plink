\name{combine.pars}
\alias{combine.pars}
\title{Combine Item Parameters for Multiple Groups}
\description{
This function combines objects of class \code{\linkS4class{irt.pars}} and/or 
\code{\link{sep.pars}} into a single \code{irt.pars} object}

\usage{
combine.pars(x, common, grp.names)
}

\arguments{
  \item{x}{an ordered list containing two or more \code{irt.pars} objects, two or more 
      \code{sep.pars} objects, or a combination of \code{irt.pars} and \code{sep.pars} 
      objects.}

  \item{common}{an n x 2 matrix or list of matrices identifying the common items between 
      adjacent pairs of objects in \code{x}. See below for more details. }

  \item{grp.names}{character vector of names for all the groups in the returned object }
 }

\details{
Although many of the methods in this package allow for lists containing \code{irt.pars} 
and \code{sep.pars} objects, it may be helpful to combine the item parameters for multiple 
groups into a single object.  \code{x} can include a combination of \code{irt.pars} and 
\code{sep.pars} objects.  The \code{irt.pars} objects can contain information for a 
single group or multiple groups.  The list elements in \code{x} should be ordered such that 
adjacent elements correspond to adjacent groups.  If an \code{irt.pars} object is included 
with multiple groups, the list element following this object should contain information for a 
group that is adjacent to the last group in the \code{irt.pars} object.

If \code{x} contains only two elements, \code{common} should be a matrix.  If \code{x} contains
more than two elements, \code{common} should be a list.  In any of the \code{common} matrices
the first column identifies the common items for the first group of two adjacent list elements
in \code{x}.  The second column in \code{common} identifies the corresponding set of common 
items from the next list element in \code{x}. For example, if \code{x} contains only two list 
elements, a single set of common items links them together.  If item 4 in group one (row 4 in 
slot \code{pars}) is the same as item 6 in group two, the first row of \code{common} would be \code{"4,6"}.

If an \code{irt.pars} object is included with multiple groups, \code{common} should correspond
to the set of common items between the last group in the \code{irt.pars} object and the
group in the adjacent element in \code{x}.
}

\value{
Returns an object of class \code{\linkS4class{irt.pars}}
}

\references{
  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{\code{\linkS4class{irt.pars}}, \code{\link{as.irt.pars}}, \code{\linkS4class{sep.pars}},
\code{\link[=sep.pars-methods]{sep.pars}} }

\keyword{ misc }