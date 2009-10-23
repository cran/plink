\name{link.con}
\alias{link.con}
\alias{link.pars}
\alias{link.ability}
\alias{get.prob}
\title{Extract Output Information}
\description{
These functions extract information from an object of class \code{\linkS4class{irt.pars}},
typically created using \code{\link[=plink-methods]{plink}}, or an object of 
class \code{\linkS4class{irt.prob}} containing expected probabilities.
}

\usage{
link.con(x, method = "ALL")
link.pars(x)
link.ability(x)
get.prob(x)
}

\arguments{
  \item{x}{an \code{R} object}

  \item{method}{character vector identifying the linking methods for which
  constants should be returned. The only acceptable values are 
    \code{ALL, MM, MS, HB, SL, and RM}. See below for details. }
}

\details{
\code{link.con} extracts the linking constants,
\code{link.pars} extracts the rescaled item parameters (if present), and 
\code{link.ability} extracts the rescaled ability estimates (if present)
\code{get.prob} extracts expected probabilities for an \code{irt.prob} object

For \code{link.con}, the following values can be included for \code{method}.
For multidimensional constants, although matrices \code{A}, \code{K}, and \code{T} are
included in the object of class "\code{link}" for the Haebara and Stocking-Lord
methods, only the matrix \code{A} will be returned for these methods.

\describe{
   \item{ALL}{Returns the constants for all linking methods}
   \item{\code{MM}:}{Mean/Mean}
   \item{\code{MS}:}{Mean/Sigma}
   \item{\code{HB}:}{Haebara}
   \item{\code{SL}:}{Stocking-Lord}
   \item{\code{RM}:}{Reckase-Martineau (for multidimensional constants only)}
}
}

\value{
Returns a matrix or list of linking constants
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{\code{\linkS4class{link}}, \code{\link[=plink-methods]{plink}} }

\keyword{ misc }