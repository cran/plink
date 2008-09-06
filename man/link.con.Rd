\name{link.con}
\alias{link.con}
\alias{link.pars}
\alias{link.ability}
\title{Extract Linking Information}
\description{
These functions extract information from an object of class "\code{link}", typically
created using \code{\link[=plink-methods]{plink}}. }

\usage{
link.con(x)
link.pars(x)
link.ability(x)
}

\arguments{
  \item{x}{an \code{R} object}
}

\value{
Returns a matrix or list
}

\details{
\code{link.con} extracts the linking constants,
\code{link.pars} extracts the rescaled item parameters (if present), and 
\code{link.ability} extracts the rescaled ability estimates (if present)
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{\code{\linkS4class{link}}, \code{\link[=plink-methods]{plink}} }

\keyword{ misc }