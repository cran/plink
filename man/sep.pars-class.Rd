\name{sep.pars-class}
\docType{class}
\alias{sep.pars-class}
\title{Class "sep.pars"}
\description{The formal S4 class for sep.pars. This class stores a set of separated item
parameters and characteristics of these parameters.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sep.pars", ...)}, but this is 
not encouraged. Use the function \code{\link{sep.pars}} instead.
}

\section{Slots}{
  \describe{
    \item{\code{a}:}{matrix of discrimination/slope parameters}

    \item{\code{b}:}{matrix of difficulty, threshold, step, and category parameters 
	(depending on the associated IRT model)}

    \item{\code{c}:}{matrix of lower asymptote parameters and category proportions for the 
	dichotomous response models and multiple-choice model respectively. \code{c} is equal to 
	\code{NA} for each item for all other models.}

    \item{\code{cat}:}{vector identifying the number of categories associated with each item. 
	All dichotomous items will have \code{cat} values equal to 2. Graded response model and 
	partial credit/generalized partial credit model items will have \code{cat} values equal 
	to the number of step/threshold parameters plus one. Nominal response model items will 
	have \code{cat} values equal to the number of categories, and multiple-choice model items 
	will have \code{cat} values equal to the number of categories plus one (the 'do not know' 
	category). }

    \item{\code{n}:}{vector identifying the total number of items, the total number of dichotomous 
	and polytomous items, and the number of items associated with each polytomous model.}

    \item{\code{mod.lab}:}{character vector of labels for the model(s).} 

    \item{\code{location}:}{logical value. If \code{TRUE}, the step and/or threshold parameters in 
	slot \code{b} for the graded response model and generalized partial credit model include a 
	location parameter. }

    \item{\code{dimensions}:}{numeric value identifying the number of modeled dimensions.}

    \item{\code{model}:}{character vector identifying all the models associated with the 
	corresponding set of item parameters. The only acceptable models are 
	\code{drm}, \code{gpcm}, \code{grm}, \code{mcm}, and \code{nrm} (see class \code{poly.mod} 
      for more information).}

    \item{\code{items}:}{list with the same length as \code{model}, where each element identifies
	 the items(rows) in the corresponding set of item parameters associated with the model(s) 
	identified in \code{model}.}
  }
}

\details{
\code{a} will be an n x m matrix for n items and m dimensions if there are no nominal response model 
or multiple-choice model items.  Otherwise, if nrm or mcm items are included, \code{a} will be an 
n x (k x m) matrix with k equal to the maximum number of response categories across items. If nrm or 
mcm items are included, the discrimination/slope parameters for the dichotomous response models, the 
graded response model, partial credit model and generalized partial credit model (and the multidimensional
extensions of these models) are listed in the first m columns with all other columns filled with \code{NA}s.

\code{b} is an n x (k x m) matrix of difficulty, threshold, step, or category parameters (depending on 
the corresponding model) for n items and m dimensions with k equal to the maximum number of \code{b} 
parameters across all items. For items with \code{b} less than k x m, the row is right-filled with \code{NA}s.

\code{c} will be an n x 1 matrix for n items if there are no multiple-choice model items. Otherwise, 
if mcm items are included, \code{c} will be an n x k matrix with k equal to the maximum number of 
response categories across items minus one. If mcm items are included, the lower asymptote parameters 
for the dichotomous response models, are listed in the first column with all other columns filled with 
\code{NA}s. The \code{c} values for the 1PL, 2PL, M1PL, and M2PL equal zero.
}

\section{Extends}{
Class \code{\linkS4class{poly.mod}}, directly.\cr
Class \code{\linkS4class{list.poly}}, by class \code{poly.mod}, distance 2.
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{
\code{\link{sep.pars}}, \code{\linkS4class{irt.pars}}
}
\keyword{classes}
