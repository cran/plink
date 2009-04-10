\name{read.bilog}
\alias{read.bilog}
\alias{read.parscale}
\alias{read.testfact}
\title{Import Parameters from IRT Software}
\description{
This function imports item and/or ability parameters from BILOG-MG, PARSCALE, 
and TESTFACT.
}

\usage{
read.bilog(file, ability = FALSE, pars.only = TRUE, as.irt.pars = FALSE)

read.parscale(file, ability = FALSE, pars.only = TRUE, as.irt.pars = FALSE, loc.out = TRUE)

read.testfact(file, ability = FALSE, dimensions = 2, guessing = FALSE, bifactor = FALSE, 
    as.irt.pars = FALSE)
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

  \item{dimensions}{number of modeled dimensions}

  \item{guessing}{logical value indicating whether a guessing parameter was modeled}

  \item{bifactor}{logical value indicating whether the bifactor model was used to
      estimate the item/ability parameters}
}

\value{
Returns a data.frame or an object of class \code{\linkS4class{irt.pars}} if 
\code{as.irt.pars} = \code{TRUE}.
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\note{These functions are currently unable to handle output generated when subtests are used.}

\keyword{ utilities }
