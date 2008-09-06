\name{link-class}
\docType{class}
\alias{link-class}
\title{Class "link"}
\description{The formal S4 class for link.  This class compiles the linking constants,
item parameter descriptive statistics, and iteration/convergence information for
the moment methods and characteristic curve methods of a separate calibration linking
procedure.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("link", ...)}, but this is 
not encouraged.  Use the function \code{link[=plink-methods{plink}} instead.
}

\section{Slots}{
  \describe{
    \item{\code{MM}:}{vector for intercept and slope linking constants under the
      Mean/Mean method}

    \item{\code{MS}:}{vector for intercept and slope linking constants under the
      Mean/Sigma method}

    \item{\code{HB}:}{vector for intercept and slope linking constants under the
      Haebara method}

    \item{\code{SL}:}{vector for intercept and slope linking constants under the
      Stocking-Lord method. The value for \code{SL} can also be \code{NULL}.}

    \item{\code{descriptives}:}{data.frame or list containing the item parameter
      descriptive statistics}

    \item{\code{iterations}:}{vector of the number of optimization iterations for
      the Haebara and Stocking-Lord methods}

    \item{\code{convergence}:}{character vector identifying the type of convergence
      reached under the Haebara or Stocking-Lord optimizations.
      (see \url{http://netlib.bell-labs.com/cm/cs/cstr/153.pdf} for more
      information on the output values)}

    \item{\code{base.grp}:}{numeric value indicating the base group for the calibration}

    \item{\code{include.mcm.nrm}:}{logical vector indicating whether multiple-choice
      model or nominal response model items were included in the computation of
      the linking constants for the moment methods and characteristic curve methods
      respectively.}
  }
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{
\code{\link[=plink-methods]{plink}}
}

\keyword{classes}
