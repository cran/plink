\name{dgn}
\alias{dgn}
\docType{data}
\title{Unidimensional Mixed-Format Item Parameters}
\description{
This (unidimensional) dataset includes item parameters for two groups estimated under 
the three parameter logistic model, the generalized partial credit model, and 
the nominal response model. There are 15 common items (7 drm, 3 gpcm, 5 nrm).
This dataset is for illustrative purposes to show how items from different
response models can be mixed together and common items can be in different
positions across groups.
}

\usage{dgn}
\format{
A list of length four.  The first element is a list of length two with item
parameter estimates for groups 1 and 2 respectively.  The gpcm items do not
include a location parameter.  The second list element is a list identifying 
the number of response categories for the two groups.  The third element, 
specifies which items correspond to the different item response models for each 
group respectively.  The last element is a matrix of common items between the two 
groups.
}

\references{
The item parameters are loosely based on parameters listed in the following articles

  Kim, J.-S. (2006). 
  Using the Distractor Categories of Multiple-Choice Items to Improve IRT Linking. 
  \emph{Journal of Educational Measurement}, \emph{43}(3), 193-213.

  Kim, S. & Lee, W.-C. (2006). 
  An Extension of Four IRT Linking Methods for Mixed-Format Tests. 
  \emph{Journal of Educational Measurement}, \emph{43}(1), 53-76.
}

\keyword{datasets}
