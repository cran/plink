\name{plink-package}
\alias{plink-package}
\docType{package}
\title{Separate Calibration Linking Methods}
\description{This package uses unidimensional item response theory methods to compute 
  linking constants and conduct chain linking of tests for multiple groups under a 
  nonequivalent groups common item design.  It allows for mean/mean, mean/sigma, Haebara, 
  and Stocking-Lord calibrations of single-format or mixed-format dichotomous (1PL, 2PL, 
  and 3PL) and polytomous (graded response partial credit/generalized partial credit, 
  nominal, and multiple-choice model) common items.
}

\details{
The primary purpose of this package is to use IRT methods to place scores and item parameters
for two or more groups onto a common scale. These methods can be used for horizontal equating
or vertical scaling. The linking methodology used in this package is commonly referred to as 
separate calibration.  The package allows for chain linking of multiple-groups with the base
scale specified for any group.  Unique and common items can can be included and rescaled using
any of the linking constants (if desired). Ability estimates can also be included and rescaled.
Symmetric or non-symmetric optimizations can be specified, various scoring functions can be
used for the Stocking-Lord method, and there is great flexibility in specifying thetas and 
theta weights to be integrated over in the characteristic curve methods.

In addition to computing linking constants and the option to rescale ability and item parameters,
the methods in the package can be used to compute item/category response probabilities and 
create plots of item/category characteristic curves. 

The package is designed to allow a variety of formats for the item parameters including vectors,
lists, matrices, and other objects (\code{irt.pars} and \code{sep.pars}) available in the 
package. Item parameters and calibration output can be summarized, and descriptive statistics
for the item parameters can be displayed as well.

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

\author{ Jonathan Weeks \email{weeksjp@gmail.com}} }

\references{ 
  Haebara, T. (1980). 
  Equating logistic ability scales by a weighted least squares method. 
  \emph{Japanese Psychological Research}, /emph{22}(3), 144-149.

  Kim, S. & Lee, W.-C. (2006). 
  An Extension of Four IRT Linking Methods for Mixed-Format Tests. 
  \emph{Journal of Educational Measurement}, \emph{43}(1), 53-76.

  Kolen, M. J. & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking} (2nd ed.).
  New York: Springer

  Loyd, B. H. & Hoover, H. D. (1980). 
  Vertical Equating Using the Rasch Model. 
  \emph{Journal of Educational Measurement}, \emph{17}(3), 179-193.

  Marco, G. L. (1977). 
  Item Characteristic Curve Solutions to Three Intractable Testing Problems.
  \emph{Journal of Educational Measurement}, \emph{14}(2), 139-160.

  Stocking, M. L. & Lord, F. M. (1983). 
  Developing a common metric in item response theory. 
  \emph{Applied Psychological Measurement}, \emph{7}(2), 201-210.
}

\keyword{ package }