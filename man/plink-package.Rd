\name{plink-package}
\alias{plink-package}
\docType{package}
\title{IRT Separate Calibration Linking Methods}
\description{This package uses unidimensional and multidimensional item response theory methods 
  to compute linking constants and conduct chain linking of tests for multiple groups under a 
  nonequivalent groups common item design.}

\details{
The package consists of three types of functions: 

\itemize{
  \item{}{Response Probability Functions: These functions compute response probabilities for specified 
  theta values (for as many dimensions and theta values as allowed by your computer's memory). The package 
  can estimate probabilities for the item response models listed below. For all of these models, I have 
  tried to include examples from the corresponding journal articles as a way to show that the equations are correct.
  \itemize{
    \item{}{1PL}
    \item{}{2PL}
    \item{}{3PL}
    \item{}{Graded Response Model* (cumulative or category probabilities)}
    \item{}{Partial Credit Model*}
    \item{}{Generalized Partial Credit Model*}
    \item{}{Nominal Response Model}
    \item{}{Multiple-Choice Model}

    \item{}{M1PL}
    \item{}{M2PL}
    \item{}{M3PL}
    \item{}{MD Graded Response Model* (cumulative or category probabilities)}
    \item{}{MD Partial Credit Model*}
    \item{}{MD Generalized Partial Credit Model*}
    \item{}{MD Nominal Response Model}
    \item{}{MD Multiple-Choice Model (this model has not formally been presented in the literature)}
  }
  * These models can be specified using a location parameter
  }

  \item{}{Linking Function: There is only one linking function, but it includes a variety of linking methods. 
  The most notable feature of this function (as compared to other software) is that it allows you to input 
  parameters for more than two groups and then chain link all of the tests together in a single run. 
  Below are several options.
  \itemize{
    \item{}{Symmetric or Non-Symmetric linking (as originally presented by Haebara)}
    \item{}{Specification of weights (can include uniform, quadrature, and normal density 
    weights for default or specified theta points)}
    \item{}{Choice of base group}
    \item{}{Choice of method to rescale item parameters and/or ability estimates (if included)}
  }

  When the item parameters correspond to a unidimensional model, linking constants can be estimated 
  using the following methods
  \itemize{
    \item{}{Mean/Mean}
    \item{}{Mean/Sigma}
    \item{}{Haebara}
    \item{}{Stocking-Lord}
  }

  When the item parameters correspond to a multidimensional model, linking constants can be estimated 
  using the following methods (with different options for the first two)
  \itemize{
    \item{}{Multidimensional extension of Haebara}
    \item{}{Multidimensional extension of Stocking-Lord}
    \item{}{Reckase-Martineau (based on the oblique procrustes rotation and least squares method 
    for estimating the translation vector \code{m}.}
  }

  For the first two approaches there is a "dilation" argument where
  \itemize{
    \item{}{all of the elements of the rotation/scaling matrix \code{A} and the translation vector \code{m} 
    are estimated via the optimization routine (see Oshima, Davey, & Lee, 2000)}
    \item{}{an orthongonal procrustes rotation to resolve the rotational indeterminacy and 
    a single dilation parameter is estimated (see Li & Lissitz, 2000)}
    \item{}{an orthongonal procrustes rotation to resolve the rotational indeterminacy and different 
    parameters are estimated to adjust the scale of each dimension (see Min, 2003). }
  }

  For all of these approaches an optional matrix can be specified to identify the ordering of dimensions. 
  It allows you to specify an overall factor structure that can take into account construct shift or 
  different orderings of the factor loadings.}

  \item{}{Utility Functions:}
  \itemize{
    \item{}{Import item parameters and/or ability estimates from BILOG-MG 3, PARSCALE 4, MULTILOG 7, TESTFACT 4,
    ICL, BMIRT, and the \code{eRm} and \code{ltm} packages.}

    \item{}{Separate item parameters from matrices or lists into the slope, difficulty, category, 
    etc. values for use in computing response probabilities, descriptive statistics, etc.}

    \item{}{Combine multiple sets of parameters into a single object for use in the linking function 
    (this essentially creates a blueprint of all the items including the IRT model used for each item, 
    the number of response categories, and the mapping of common items between groups)}

    \item{}{Plot item characteristic/category curves (for unidimensional items), item response/category 
    surfaces, contour plots, level plots, and three types of vector plots (for multidimensional items).
    The first three types of multidimensional plots include functionality for creating "conditional" 
    surfaces when there are more than two dimensions).}

    \item{}{Summarize the item parameters (unique and/or common) for each item response model separately, and overall}
  }
}

Running the separate calibration is typically a two-step process.  The first step is to format the
item parameters and the second step is to run the function \code{\link[=plink-methods]{plink}}. In the
simplest scenario, the parameters should be formatted as a single \code{\linkS4class{irt.pars}} 
object with multiple groups. Refer to the function \code{\link{as.irt.pars}} for specific details. Once 
in this format, the linking constants can be computed using \code{plink}.  The 
\code{\link[=summary.irt.pars]{summary}} function can be used to summarize the common item parameters 
(including descriptive statistics) and the linking constants.

To compute response probabilities for a given model, the following functions can be used: 
\code{\link[=drm-methods]{drm}}, \code{\link[=gpcm-methods]{gpcm}}, \code{\link[=grm-methods]{grm}},
\code{\link[=mcm-methods]{mcm}}, or \code{\link[=nrm-methods]{nrm}}. The \code{\link[=plot.irt.prob]{plot}} 
function can be used to create item/category characteristic curves, item/category response surfaces, 
and vector plots.
}


\references{
  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com}}

\keyword{ package }
