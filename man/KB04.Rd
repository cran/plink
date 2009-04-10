\name{KB04}
\alias{KB04}
\docType{data}
\title{Kolen - Brennan (2004)}
\description{
This (unidimensional) dataset includes three parameter logistic model (3PL) item parameter 
estimates for two groups.  The Form X items are from a "new" form and the Form Y items are 
from an "old" form.  The data are listed in Table 6.5 in \emph{Test Equating, Scaling, and 
Linking} (2nd ed.).  There are 12 common items between the two forms.  
}
\usage{KB04}
\format{
A list of length two.  The first element \code{pars} is a list of length two. 
Each of the list elements in \code{pars} is a 36 x 3 matrix of item parameters.
Element one is for \code{Form X}, and element two is for \code{Form Y}.

The second list element in \code{KB04} is a matrix identifying the common items
between the two sets of parameters in \code{pars}.  That is, the first column in
\code{common} identifies the rows in \code{pars$form.x} that are common items.
The second column identifies the corresponding set of common items in \code{pars$form.y}.
}

\note{In the book, the authors transform the Form X parameters to the Form Y scale.}

\source{
  Kolen, M. J., & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking.} (2nd ed.), p. 192
  New York: Springer
}

\keyword{datasets}
