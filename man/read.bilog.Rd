\name{read.bilog}
\alias{read.bilog}
\alias{read.parscale}
\title{Import Parameters from IRT Software}
\description{
This function imports item and/or ability parameters from BILOG-MG and PARSCALE.
}

\usage{
read.bilog(file, ability=FALSE, pars.only=TRUE, as.irt.pars=FALSE)

read.parscale(file, ability=FALSE, pars.only=TRUE, as.irt.pars=FALSE, loc.out=TRUE)
}

\arguments{
  \item{file}{filename of \code{.PAR} or \code{.SCO} file to be imported}

  \item{ability}{if \code{TRUE}, \code{file} is an \code{.SCO} file}

  \item{pars.only}{if \code{TRUE}, only the item parameters will be imported (i.e.,
	any other information like standard errors will be dropped)}

  \item{as.irt.pars}{if \code{TRUE}, the parameters will be output as an 
	\code{\linkS4class{irt.pars}} object (this is only applicable to item
	parameters)}

  \item{loc.out}{if \code{TRUE}, the step/threshold parameters will be reformated to 
	be deviations from a location parameter }
}

\value{
Returns a data.frame or an object of class \code{\linkS4class{irt.pars}} if 
\code{as.irt.pars} = \code{TRUE}.
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\note{These functions are currently unable to handle output generated when subtests are used.}

\keyword{ utilities }
