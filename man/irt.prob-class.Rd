\name{irt.prob-class}
\docType{class}
\alias{irt.prob-class}
\title{Class "irt.prob"}
\description{The formal S4 class for irt.prob.  This class contains the expected 
probabilities of responding in a given category for a set of items and theta values 
under the specified IRT models.  The class also includes characteristics of the items.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("irt.prob", ...)}, but this is
not encouraged.  Use one of the functions \code{drm, gpcm, grm, mcm, nrm, or mixed} instead.
}

\section{Slots}{
  \describe{
      \item{\code{prob}:}{data.frame of item probabilities with dimensions n x j+1 for n theta 
        values, for j items (the first column contains theta values)}

      \item{\code{p.cat}:}{vector identifying the number of categories for each item for which 
        probabilities were computed}

      \item{\code{mod.lab}:}{character vector of labels for the model(s).}

      \item{\code{model}:}{character vector identifying all the models used to compute the 
        probabilities in \code{prob}.  The only acceptable models are 
        \code{drm, gpcm, grm, mcm, and nrm} (see class \code{\linkS4class{poly.mod}} for more 
        information).}

    \item{\code{items}:}{list with the same length as \code{model}, where each element 
	identifies the items associated with the model(s) specified in \code{model}.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{poly.mod}"}, directly.\cr
Class \code{"\linkS4class{list.poly}"}, by class "poly.mod", distance 2.
}

\author{ Jonathan Weeks \email{weeksjp@gmail.com} }

\seealso{
\code{\link[=drm-methods]{drm}:} for dichotomous response models (1PL, 2PL, and 3PL)\cr
\code{\link[=gpcm-methods]{gpcm}:} for the partial credit/generalized partial credit model\cr
\code{\link[=grm-methods]{grm}:} for the graded response model\cr
\code{\link[=mcm-methods]{mcm}:} for the multiple-choice model\cr
\code{\link[=nrm-methods]{nrm}:} for the nominal response model\cr
\code{\link[=mixed-methods]{mixed}:} for mixed-format items\cr
\code{\linkS4class{irt.pars}}
}

\keyword{classes}
