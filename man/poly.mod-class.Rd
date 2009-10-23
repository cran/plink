\name{poly.mod-class}
\docType{class}
\alias{poly.mod-class}
\title{Class "poly.mod"}
\description{The formal S4 class for poly.mod.  This class characterizes the models and 
associated items for a set of item parameters.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("poly.mod", ...)}, but this is 
not encouraged. Use the function \code{\link{as.poly.mod}} instead.
}

\section{Slots}{
  \describe{
    \item{\code{model}:}{character vector identifying all the models associated with the 
	corresponding set of item parameters.  The only acceptable models are 
	\code{drm}, \code{gpcm}, \code{grm}, \code{mcm}, and \code{nrm}.  See below for more 
      details.}

    \item{\code{items}:}{list with the same length as \code{model}, where each element 
	identifies the items(rows) in the corresponding set of item parameters associated 
	with the model(s) identified in \code{model}.}
  }
}

\details{ 
    The IRT models associated with the codes:
  \describe{
    \item{\code{drm}:}{dichotomous response models (includes 1PL, 2PL, and 3PL)}
    \item{\code{gpcm}:}{generalized partial credit model (includes the
	partial credit model)}
    \item{\code{grm}:}{graded response model}
    \item{\code{mcm}:}{multiple-choice model}
    \item{\code{nrm}:}{nominal response model}
  }
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\note{
The names of the list elements for \code{items} must correspond to the elements in \code{model}. 
For example, if the poly.mod object is \code{pm} and the first element in \code{pm@model} is 
\code{drm},one should be able to reference the associated items by \code{pm@items$drm}. If 
the list elements are unnamed, some functions may not work properly.
}

\seealso{
\code{\link{as.poly.mod}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}},
\code{\linkS4class{irt.prob}}
}

\keyword{classes}
