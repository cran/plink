\name{link-class}
\docType{class}
\alias{link-class}
\title{Class "link"}
\description{The formal S4 class for link.  This class compiles the linking constants,
item parameter descriptive statistics, and iteration/objective/convergence information
for the moment methods and characteristic curve methods from a separate calibration linking
procedure. It is also able to store rescaled item and ability parameters.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("link", ...)}, but this is 
not encouraged.  Use the function \code{\link[=plink-methods]{plink}} instead.
}

\section{Slots}{
  \describe{
    \item{\code{constants}:}{list of linking methods where the list elements include
      the linking constants.  For unidimensional models, each element is a vector of
      length two containing a slope and intercept to adjust the SD and mean respectively. 
      For multidimensional models, each element is a list with a matrix (or set of 
      matrices) to resolve rotational indeterminacy and scale, and a vector of 
      translation constants to adjust the mean of each dimension. }

    \item{\code{descriptives}:}{data.frame or list containing the item parameter
      descriptive statistics}

    \item{\code{iterations}:}{vector of the number of optimization iterations for
      the Haebara and Stocking-Lord methods}

    \item{\code{objective}:}{vector of criterion values for HBcrit and SLcrit at the point
      of convergence for the Haebara and Stocking-Lord methods}

    \item{\code{convergence}:}{character vector identifying the type of convergence
      reached under the Haebara or Stocking-Lord optimizations.
      (see \url{http://netlib.bell-labs.com/cm/cs/cstr/153.pdf} for more
      information on the output values)}

    \item{\code{base.grp}:}{numeric value indicating the base group for the calibration}

    \item{\code{n}:}{vector identifying the total number of common items, the number of
      dichotomous common items, and the number of polytomous common items}

    \item{\code{grp.names}:}{character vector of group names}

    \item{\code{mod.lab}:}{character vector identifying the dichotomous and/or polytomous
      models used to model the item responses}

    \item{\code{dilation}:}{character value identifying the dilation approach used when
      estimating multidimensional linking constants}
  }
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{ 
\code{\link[=plink-methods]{plink}}
}

\keyword{classes}
