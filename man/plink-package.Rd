\name{plink-package}
\alias{plink-package}
\docType{package}
\title{IRT Separate Calibration Linking Methods}
\description{This package uses item response theory methods to compute linking constants and 
  conduct chain linking of tests for multiple groups under a nonequivalent groups common item 
  design. It allows for mean/mean, mean/sigma, Haebara, and Stocking-Lord calibrations of 
  dichotomous (1PL, 2PL, and 3PL) and/or polytomous (graded response, partial credit/generalized 
  partial credit, nominal, and multiple-choice model) common items.
}

\details{
Running the separate calibration is typically a two-step process. The first step is to
format the item parameters for processing with \code{\link{plink}}. Parameters should be formatted
as a single \code{\linkS4class{irt.pars}} object with multiple groups, a set of \code{irt.pars}
objects, or a set of \code{\linkS4class{sep.pars}} object.  Once in this format, response
probabilities can be computed using \code{drm, gpcm, grm, mcm, or nrm} or linking constants
can be computed using \code{plink}. 

The functions \code{\link{as.irt.pars}}, \code{\link{sep.pars}}, and \code{\link{combine.pars}}
can be used to create the \code{irt.pars} and \code{sep.pars} objects. The function \code{summary}
can be used to summarize item parameters (including descriptive statistics) and linking constants,
and \code{plot} can be used to create item/category characteristic curves.
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com}} }

\keyword{ package }