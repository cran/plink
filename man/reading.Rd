\name{reading}
\alias{reading}
\docType{data}
\title{Multi-Group Mixed-Format Reading Assessment}
\description{
This (unidimensional) dataset includes item parameters from a large-scale reading 
assessment.  The parameters were estimated using a combination of the three parameter 
logistic model (3PL) and the generalized partial credit model (GPCM).  There are 
six sets of parameters which are based on tests administered in four grades over 
three years.  In particular, the data include a grade 3 and grade 4 test in year 0, 
a grade 4 and grade 5 test in year 1, and a grade 5 and grade 6 test in year 2. The 
label for the information related to the grade 3 test is grade3.0 where the value 
after the decimal indicates the year.  Similar labels are used for the other grade/year
combinations.  This dataset is used for illustrative purposes to show in a multi-group
scenario how items from different response models can be mixed together and common 
items can be in different positions across groups.
}

\usage{reading}
\format{
A list of length four.  The first element is a list of length six with item
parameter estimates for each grade/year.  There is no location parameter for the GPCM
items.  The second list element is a list identifying the number of response categories 
for the six grade/year combinations.  The third element specifies which items correspond 
to the different item response models for each grade/year respectively.  The last element 
is a list of common item matrices for each adjacent grade/year group.
}

\source{
  Briggs, D. C. & Weeks, J. P. (In Press)
  The Impact of Vertical Scaling Decisions on Growth Interpretations.
  \emph{Educational Measurement: Issues and Practices}.
}

\keyword{datasets}
