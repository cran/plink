\name{act.nrm}
\alias{act.nrm}
\docType{data}
\title{ACT Mathematics Test 1997-1998 Nominal Response Model}
\description{
This (unidimensional )dataset includes nominal response model item parameter estimates 
for two groups.  There are 60 items, all of which are common items.  The estimates
are from the 1997 and 1998 ACT mathematics test.  The item parameters were
obtained from the examples for the software \code{nreq} written by Bradley
Hanson.
}
\usage{act.nrm}
\format{
A list of length two. The matrices in the two list elements contain
nominal response model item parameters for 1997 and 1998 respectively. The first
six columns in each matrix are the slope parameters and the next six columns are
the category difficulty parameters.
}

\source{
  Kim, J.-S. & Hanson, B. A. (2002). 
  Test Equating Under the Multiple-Choice Model.
  \emph{Applied Psychological Measurement}, \emph{26}(3), 255-270.

The item parameters can be found here
\url{http://www.b-a-h.com/software/mcmequate/index.html}
}

\keyword{datasets}
