\name{TK07}
\alias{TK07}
\docType{data}
\title{Tong - Kolen (2007)}
\description{
This (unidimensional) dataset includes three parameter logistic model (3PL) item parameter 
estimates from simulated dataset based on items from the 1992 Iowa Test of Basic 
Skills vocabulary test.  The data are for six groups (grades 3-8) with varying 
numbers of common items between adjacent grades. 
}

\usage{TK07}
\format{
A list of length two.  The first element \code{pars} is a list of length six.  Each
of the list elements in \code{pars} is a n x 3 matrix of item parameters where n
differs for each grade.  The second list element in \code{TK07} is a list of length
five identifying the common items between adjacent groups in \code{pars}.  Each list
element in \code{common} is a j x 2 matrix for j common items.  For example, there
are 13 common items between grades 3 and 4. The first list element in \code{common} 
is for the common items between these grades.  The first column in the matrix show that
items 14-26 in grade 3 are the same as items 1-13 in grade 4.  Similar matrices are 
identified for all other grade pairs.
}

\source{
The parameters were obtained from Dr. Ye Tong at Pearson Educational Measurement
and were used in the following article:

Tong, Y. & Kolen, M. J. (2007)
Comparisons of methodologies and results in vertical scaling for educational 
achievement tests. \emph{Applied Measurement in Education}, \emph{20}(2), 227-253.
}

\keyword{datasets}
