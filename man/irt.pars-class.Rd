\name{irt.pars-class}
\docType{class}
\alias{irt.pars-class}
\alias{initialize,irt.pars-method}
\title{Class "irt.pars"}
\description{The formal S4 class for irt.pars. This class contains the item parameters
and characteristics of the item parameters for one or more groups.  When parameters
are included for two or more groups, the common items between the different groups are also
included.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("irt.pars", ...)}, but this is 
not encouraged.  Use the function \code{\link{as.irt.pars}} instead. 
}
\section{Slots}{
  \describe{
    \item{\code{pars}:}{matrix of item parameters for a single group or a list of matrices 
      containing item parameters for two or more groups}

    \item{\code{cat}:}{vector with length equal to the number of items, identifying the number 
      of response categories for each item for a single group or a list of response category 
      vectors for two or more groups}

    \item{\code{poly.mod}:}{a \code{\linkS4class{poly.mod}} object for one group or a list of 
      \code{poly.mod} objects identifying the items associated with each IRT model (see class 
      \code{\linkS4class{poly.mod}} for more information)}

    \item{\code{common}:}{a j x 2 matrix for j common items or a list of matrices identifying the 
      common items between pairs of item parameters}

    \item{\code{location}:}{logical vector identifying whether a given set of item parameters 
      (i.e., for a given group) includes a location parameter }

    \item{\code{groups}:}{the number of groups (i.e. sets of item parameters)}

    \item{\code{dimensions}:}{numeric vector identifying the number of modeled dimensions}
  }
}

\details{
Objects of class \code{irt.pars} contain all the information necessary to produce
expected response probabilities, linking constants (when data for two or more groups are 
included), and to plot item response curves/surfaces, vector plots, and parameter 
comparison plots.
}

\section{Extends}{
Class \code{\linkS4class{poly.mod}}, directly.
Class \code{\linkS4class{list.poly}}, by class \code{poly.mod}, distance 2.
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{
\code{\link{as.irt.pars}}, \code{\link{as.poly.mod}}
}

\keyword{classes}
