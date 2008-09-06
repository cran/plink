\name{sep.pars-methods}
\docType{methods}
\alias{sep.pars}
\alias{sep.pars-methods}
\alias{sep.pars,numeric-method}
\alias{sep.pars,data.frame-method}
\alias{sep.pars,irt.pars-method}
\alias{sep.pars,matrix-method}
\alias{sep.pars,list-method}
\title{Separate Item Parameters}
\description{
This function splits the item parameters in the specified object into
discrimination parameters, difficulty/step/threshold/category parameters,
and lower asymptote/category probability parameters.
}

\usage{
sep.pars(x, cat, poly.mod, location = FALSE, loc.out = FALSE, ...)

## S4 method for signature 'numeric'
sep.pars(x, cat, poly.mod, location, loc.out, ...)

## S4 method for signature 'matrix'
sep.pars(x, cat, poly.mod, location, loc.out, ...)

## S4 method for signature 'data.frame'
sep.pars(x, cat, poly.mod, location, loc.out, ...)

## S4 method for signature 'irt.pars'
sep.pars(x, cat, poly.mod, location, loc.out, ...)

## S4 method for signature 'list'
sep.pars(x, cat, poly.mod, location, loc.out, ...)
}

\arguments{
  \item{x}{Object containing item parameters. See below for more details. }

  \item{cat}{write stuff here }

  \item{poly.mod}{write stuff here}

  \item{location}{if \code{TRUE}, the step parameters are deviations from a difficulty
  	parameter }

  \item{loc.out}{write stuff here }

  \item{...}{further arguments passed to or from other methods }
}

\section{Methods}{
\describe{

\item{x = "numeric"}{This method only applies to a single set of item parameters that 
  correspond to the Rasch model.  \code{x} is assumed to be a vector of item difficulties.  
  Discrimination parameters are set to one and the lower asymptote values are set to zero 
  for all items.}

\item{x = "matrix"}{\code{x} can include item parameters from multiple models.  The general 
  format for structuring \code{x} is an additive column approach.  That is, the left-most 
  columns are typically for discrimination parameters, the next column, if applicable, is 
  for location parameters, the next set of columns is for difficulty/threshold/step/category 
  parameters, and the final set of columns is for guessing parameters.  When multiple models 
  are included, or models with differing numbers of response categories, not all cells in 
  \code{x} will have data.  In these instances, cells with no data should be \code{NA}.
  \describe{
    \item{1PL:}{For the 1PL model with discriminations equal to 1 (Rasch Model), a single 
      column of item difficulties can be supplied.  Two columns can also be used with
      all the values in the first column equal to 1 and difficulty parameters in the second
      column.  For discrimination values other than 1, \code{x} should contain at least two 
      columns, the first for item discriminations (identical for all items) and the second 
      for item difficulties.  The lower asymptote defaults to zero for all items; however, a 
      third column of zeros can be included.}

    \item{2PL:}{\code{x} should include at lease two columns, the first for item 
      discriminations and the second for item difficulties.  The lower asymptote defaults to 
      zero for all items; however, a third column of zeros can be included.}

    \item{3PL:}{\code{x} should include three columns, the first for item discriminations, 
      the second for item difficulties, and the third for lower asymptote values.}
    
    \item{Partial Credit Model:}{For a model with discrimination parameters equal to one,
      a column with the discrimination values does not need to be included; however, if
      the discrimination values do not equal one, the first column should include
      discrimination values (equal for all items).  The next column (or the first column if
      no discrimination are included), if applicable, if for location parameters. The 
      last set of columns is for step or step deviation parameters.  If no discrimination
      parameters and no location parameters are included the step parameters will begin
      in column one.  If either discrimination or location parameters are include, the step
      or step deviation parameters will begin in the second column.  If both discrimination 
      and location parameters are included, the step/step deviation parameters will begin
      in the third column.}

    \item{Generalized Partial Credit Model:}{The first column is for discrimination
      parameters.  If the step parameters are deviations from a location parameter, the 
      location parameter should go in the second column and the deviation step parameters
      should begin in the third column.  If no location parameter is included, the step 
      parameters should begin in the second column.}

    \item{Graded Response Model:}{The first column is for discrimination parameters.  If 
      the threshold parameters are deviations from a location parameter, the location 
      parameter should go in the second column and the deviation thresholds should begin
      in the third column.  If no location parameter is included, the threshold parameters 
      should begin in the second column.}

    \item{Nominal Response Model:}{The first k columns are for the item discrimination
      parameters and the second k columns are for category difficulty parameters.  If
      items have different numbers of response categories, the category difficulties 
      will not necessarily start in the column following the discrimination parameters.
      For k equal to the maximum number of response categories across all NRM items,
      the difficulty parameters for all NRM items should begin in column k+1.  Missing
      cells for the discrimination parameters should be \code{NA}.}

    \item{Multiple-Choice Model:}{The first k columns are for the item discrimination
      parameters, the second k columns are for category difficulty parameters, and the
      last k-1 columns are for guessing probabilities.  If items have different numbers 
      of response categories, the category difficulties and guessing probabilities
      will not necessarily start in the column following the discrimination or
      difficulty parameters.  For k equal to the maximum number of response categories 
      across all MCM items, the difficulty parameters for all MCM items should begin 
      in column k+1, and guessing probabilities should begin in column 2k+1.  Missing
      cells for the discrimination parameters, difficulty parameters, and guessing
      probabilies should be \code{NA}.}
  }
}

\item{x = "data.frame"}{See the method for x = "matrix" }

\item{x = "irt.pars"}{\code{x} is an object of class \code{\linkS4class{irt.pars}}

\item{x = "list"}{This method can include a list with one, two, or three 
  elements.  In general, these elements correspond to discrimination, difficulty, and guessing
  parameters, although this may not be the case depending on the supplied parameters.  If a
  combination of models are used, the number of list elements should correspond with the 
  maximum number of elements across models.  For example, if the 3PL model (3 list elements) 
  and nominal response model (2 list elements) are used, the list should include three
  elements. 

  If the nominal response model or multiple-choice model are used in conjunction with 
  another model with only a single discrimination parameter, the first list element should
  be a matrix with discrimination values for the single discrimination value models in the
  first column and \code{NA}s in the remaining columns. Similarly, if the multiple-choice
  model is included, lower asymptote parameters for the 3PL model should be included in the
  first column of the matrix of guessing probabilities in the third list element with \code{NA}s
  for all other columns. For all other models, all the columns should be \code{NA}.
 
  \describe{
    \item{1PL:}{For the 1PL model with discriminations equal to 1 (Rasch Model), one element
      with item difficulties can be supplied.  Alternatively, two elements can be used with
      the first list element containing a matrix/vector of ones and difficulty parameters in 
      the second list element.  For discrimination values other than 1, \code{x} should contain 
      at  least two list elements, the first for item discriminations (identical for all items) 
      and the second for item difficulties.  The lower asymptote defaults to zero for all items; 
      however, a third element with a vector/matrix of zeros can be included.

      If \code{x} includes parameters for multiple models and the number of list elements is
      two, the first element should include a vector of ones.  If the number of list elements
      is three, the third element should include a vector of zeros.}

    \item{2PL:}{\code{x} should include at lease two list elements, the first for item 
      discriminations and the second for item difficulties.  The lower asymptote defaults to 
      zero for all items; however, a third element with a vector/matrix of zeros can be included.

      If \code{x} includes parameters for multiple models and the number of list elements 
      is three, the third element should include a vector of zeros.}

    \item{3PL:}{\code{x} should include three list elements, the first for item discriminations, 
      the second for item difficulties, and the third for lower asymptote values.}

    \item{Partial Credit Model:}{For a partial credit model with discrimination values equal
      to one, one list element with step/step deviation parameters can be supplied.  If the 
      discrimination is not one, two list elements should be included.  The first list element 
      should contain a vector of discrimination values (equal for all item) and the second 
      element should include the step/step deviation parameters.

      The list element containing the step or step deviation parameters should be a matrix.
      If the parameters are deviations from a location parameter, the location parameter
      must be in the first column with the remaining columns for the step deviation values.
      If no location parameter is included, the step parameters should begin in column one.
      (See the method for x = 'matrix' above for more information on the formatting for this
      object). 
        
      If \code{x} includes multiple models and the number of list elements is two, the first
      list element should include a vector of ones or a vector of constant discrimination
      values.  If there are three columns, the third list element should conatain a vector 
      of \code{NA} for all the PCM items.}

    \item{Generalized Partial Credit Model:}{This model should include at least two list
      elements.  The first element should include a vector of discrimination parameters and 
      the second element should include a matrix of step/deviation step parameters. If a
      location parameter is included, it must be in the first column with the remaining columns 
      containing the step deviation values. If no location parameter is included, the step 
      parameters should begin in column one. (See the method for x = 'matrix' above for more 
      information on the formatting for this object). 
        
      If \code{x} includes multiple models and there are three list elements, the third element 
      should conatain a vector of \code{NA} for all the GPCM items.}

    \item{Graded Response Model:}{This model should include at least two list elements.  The 
      first element should include a vector of discrimination parameters and the second element 
      should include a matrix of threshold/deviation threshold parameters.  If a location parameter 
      is included, it must be in the first column with the remaining columns containing the 
      threshold deviation values.  If no location parameter is included, the threshold
      parameters should begin in column one. (See the method for x = 'matrix' above for more 
      information on the formatting for this object).}
        
      If \code{x} includes multiple models and there are three list elements, the third element 
      should conatain a vector of \code{NA} for all the GRM items.}

    \item{Nominal Respons Model:}{There should be two list elements for this model. The 
      first element is a matrix of discrimination parameters and the second element is for
      category difficulty parameters (See the method for x = 'matrix' above for more 
      information on the formatting for these objects). If the multiple-choice model is 
      included, all the columns in the matrix for the third list element should be \code{NA}
      for the NRM items.}

    \item{Multiple-Choice Model:}{There should be three list elements for this model. The 
      first element is a matrix of discrimination parameters, the second element is for
      category difficulty parameters, and the third element is for guessing probabilities.
      (See the method for x = 'matrix' above for more information on the formatting for these 
      objects).}
  }
}
}}

\value{
  Returns an object of class \code{\linkS4class{sep.pars}}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\examples{
# Create object for three dichotomous (1PL) items with difficulties -1, 0, 1
x <- sep.pars(c(-1,0,1))


# Create object for three dichotomous (3PL) items and two polytomous (gpcm) items
# without a location parameter (use signature matrix, missing)
dichot <- matrix(c(1.2, .8, .9, 2.3, -1.1, -.2, .24, .19, .13),3,3)
poly <- matrix(c(.64, -1.8, -.73, .45, NA, .88, .06, 1.4, 1.9, 2.6),2,5,byrow=TRUE)
pars <- rbind(cbind(dichot,matrix(NA,3,2)),poly)
cat <- c(2,2,2,4,5)
pm <- as.poly.mod(5, c("drm","gpcm"), list(1:3,4:5))
x <- sep.pars(pars, cat, pm)
summary(x)


# Create object for three dichotomous (3PL) items and two polytomous (gpcm) items
# without a location parameter
a <- c(1.2, .8, .9, .64, .88)
b <- matrix(c(
  2.3, rep(NA,3),
  -1.1, rep(NA,3),
  -.2, rep(NA,3),
  -1.8, -.73, .45, NA,
  .06, 1.4, 1.9, 2.6),5,4,byrow=TRUE)
c <- c(1.4, 1.9, 2.6, NA, NA)
pars <- list(a,b,c)
cat <- c(2,2,2,4,5)
pm <- as.poly.mod(5, c("drm","gpcm"), list(1:3,4:5))
x <- sep.pars(pars, cat, pm)
summary(x)


# Create object for three dichotomous (3PL) items, four polytomous items,
# two gpcm items and two nrm items. Include a location parameter for the 
# gpcm items. Maintain the location parameter in the output.
a <- matrix(c(
  1.2, rep(NA,4),
  .8, rep(NA,4),
  .9, rep(NA,4),
  .64, rep(NA,4),
  .88, rep(NA,4),
  .905, .522, -.469, -.959, NA, 
  .828, .375, -.357, -.079, -.817),7,5,byrow=TRUE)
b <- matrix(c(
  2.3, rep(NA,4),
  -1.1, rep(NA,4),
  -.2, rep(NA,4),
  -.69, -1.11, -.04, 1.14, NA,
  1.49, -1.43, -.09, .41, 1.11,
  .126, -.206, -.257, .336, NA, 
  .565, .865, -1.186, -1.199, .993),7,5,byrow=TRUE)
c <- c(.14, .19, .26, rep(NA,4))
pars <- list(a,b,c)
cat <- c(2,2,2,4,5,4,5)
pm <- as.poly.mod(7, c("drm","gpcm","nrm"), list(1:3,4:5,6:7))
x <- sep.pars(pars, cat, pm, location=TRUE, loc.out=TRUE)
summary(x, TRUE)


# Create irt.pars object with two groups then run sep.pars
pm <- as.poly.mod(36)
x <- as.irt.pars(KB04$pars, KB04$common, cat=list(rep(2,36),rep(2,36)), 
  list(pm,pm), grp.names=c("form.x","form.y"))
out <- sep.pars(x)
summary(out, TRUE)
}

\keyword{methods}
