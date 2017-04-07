\name{reckase9}
\alias{reckase9}
\docType{data}
\title{Reckase MIRT Multidimensional Item Parameters}
\description{
This dataset includes multidimensional two parameter logistic model (M2PL) item parameter 
estimates for two dimensions and two groups.  The Form 1 items are from a "new" form and 
the Form 2 items are from an "old" form.  There are 16 common items between the two forms.  
}
\usage{reckase9}
\format{
A list of length two.  The first element \code{pars} is a list of length two. 
Each of the list elements in \code{pars} is a 80 x 3 matrix of item parameters.  Element 
one is for \code{Form 1}, and element two is for \code{Form 2}.  The first two columns 
in each matrix are for the slope parameters and the third column is for the parameter 
related to item difficulty.  

The second list element in \code{reckase9} is a matrix identifying the common items
between the two sets of parameters in \code{pars}.  That is, the first column in
\code{common} identifies the rows in \code{pars$form1} that are common items.
The second column identifies the corresponding set of common items in \code{pars$form2}.
}

\source{
  Reckase, M. D. (2009).  Multidimensional item response theory.  New York: Springer.
}

\keyword{datasets}
