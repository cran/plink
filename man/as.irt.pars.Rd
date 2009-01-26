\name{as.irt.pars}
\alias{as.irt.pars}
\alias{as.irt.pars-methods}
\alias{as.irt.pars,numeric,missing-method}
\alias{as.irt.pars,data.frame,missing-method}
\alias{as.irt.pars,matrix,missing-method}
\alias{as.irt.pars,list,missing-method}
\alias{as.irt.pars,sep.pars,missing-method}
\alias{as.irt.pars,list,matrix-method}
\alias{as.irt.pars,list,list-method}
\title{irt.pars objects}
\description{
This function attempts to turn the given values into an irt.pars object
}
\usage{
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)

## S4 method for signature 'numeric', 'missing'
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)

## S4 method for signature 'matrix', 'missing'
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)

## S4 method for signature 'data.frame', 'missing'
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)

## S4 method for signature 'list', 'missing'
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)

## S4 method for signature 'sep.pars', 'missing'
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)

## S4 method for signature 'list', 'matrix'
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)

## S4 method for signature 'list', 'list'
as.irt.pars(x, common, cat, poly.mod, location = FALSE, grp.names = NULL, ...)
}

\arguments{
  \item{x}{an \code{R} object containing item parameters. See below for more details }

  \item{common}{j x 2 matrix or list of matrices identifying the common items between
    adjacent groups in \code{x}. This argument is only applicable when \code{x}
    includes two or more groups}

  \item{cat}{vector or list of the number of response categories for each item in \code{x}.
    if \code{x} is a list, \code{cat} should be a list as well. For multiple-choice
    model items, \code{cat} is the number of response categories plus one (the 
    additional category is for 'do not know')}

  \item{poly.mod}{a \code{\linkS4class{poly.mod}} object or a list of \code{poly.mod}
    objects (one for each group in \code{x})}

  \item{location}{logical vector identifying whether the parameters in \code{x}
    include a location parameter (i.e. for polytomous items)}

  \item{grp.names}{character vector of group names}

  \item{...}{further arguments passed to or from other methods}
}

\section{Methods}{
\describe{

\item{x = "numeric", common = "missing"}{This method only applies to a single set
  of item parameters that correspond to the Rasch model.  \code{x} is assumed to be
  a vector of item difficulties.  Discrimination parameters are set to one and the
  lower asymptote values are set to zero for all items.}

\item{x = "matrix", common = "missing"}{\code{x} can include item parameters from
  multiple models.  The general format for structuring \code{x} is an additive column 
  approach.  That is, the left-most columns are typically for discrimination parameters,
  the next column, if applicable, is for location parameters, the next set of columns
  is for difficulty/threshold/step/category parameters, and the final set of columns
  is for guessing parameters.  When multiple models are included, or models with 
  differing numbers of response categories, not all cells in \code{x} will have data.
  In these instances, cells with no data should be \code{NA}.
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

\item{x = "data.frame", common = "missing"}{See the method for x = "matrix" }

\item{x = "list", common = "missing"}{This method can include a list with one, two, or three 
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
      information on the formatting for this object).
        
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

\item{x = "sep.pars", common = "missing"}{\code{x} is an object of class \code{\linkS4class{sep.pars}}.
  The arguments \code{cat}, \code{poly.mod}, and \code{location} do not need to be included.}

\item{x = "list", common = "matix"}{This method is intended for the creation of an 
  object with only two groups.  Each list element should conform to one of the formats 
  listed above (i.e. a numeric vector, matrix/data.frame, list, or sep.pars object)
  or be an \code{\linkS4class{irt.pars}} object.  The format of each list element does
  not need to be the same.  That is, the first list element might be a \code{"sep.pars"} 
  object while the second element is a list of item parameters. If an object of class
  \code{irt.pars} is included, the object can contain information for a single group
  or for multiple groups.  If the \code{irt.pars} object includes multiple groups,
  \code{common} should identify the common items between the last group in \code{x[[1]]}
  and the first group in \code{x[[2]]}.
  
  For this method \code{common} is a j x 2 matrix where the first column identifies the 
  common items for the first group (\code{x[[1]]}).  The second column in \code{common} 
  identifies the corresponding set of common items from the next list element.  For example, 
  if item 4 in group one (row 4 in \code{x[[1]]}) is the same as item 6 in group two, the 
  first row of \code{common} would be \code{"4,6"}.

  If the objects in \code{x} are of class \code{sep.pars} or \code{irt.pars}, \code{cat}, 
  \code{poly.mod}, and \code{location} do not need to be included.}

\item{x = "list", common = "list"}{This method is intended for the creation of an object 
  with more than two groups; however, if there are two groups \code{common} can be 
  defined as a list with length one.  Each list element should conform to one of the formats 
  listed above (i.e. a numeric vector, matrix/data.frame, list, or sep.pars object)
  or be an \code{\linkS4class{irt.pars}} object.  The format of each list element does
  not need to be the same.  That is, the first list element might be a \code{"sep.pars"} 
  object while the second element is a list of item parameters.  If an object of class
  \code{irt.pars} is included, the object can contain information for a single group
  or for multiple groups.  The list elements in \code{x} should be ordered such that 
  adjacent elements correspond to adjacent groups.  If an \code{"irt.pars"} object is 
  included with multiple groups, the list element following this object should contain 
  information for a  group that is adjacent to the last group in the \code{"irt.pars"} object
  
  For this method \code{common} is a list of j x 2 matrices where the first column identifies 
  the identifies the common items for the first group of two adjacent list elements in
  The second column in \code{common} identifies the corresponding set of common items from the 
  next list element in \code{x}.  For example, if \code{x} contains only two list elements, a 
  single set of common items links them together.  If item 4 in group one (row 4 in slot \code{pars}) 
  is the same as item 6 in group two, the first row of \code{common} would be \code{"4,6"}.

  If all the objects in \code{x} are of class \code{sep.pars} or \code{irt.pars}, the arguments
  \code{cat}, \code{poly.mod}, and \code{location} do not need to be included.}
	}
}

\value{
Returns an object of class \code{\linkS4class{irt.pars}}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\seealso{\code{\linkS4class{irt.pars}}, \code{\link{as.poly.mod}}, \code{\linkS4class{poly.mod}}, 
\code{\linkS4class{sep.pars}} }

\examples{
# Create object for three dichotomous (1PL) items with difficulties -1, 0, 1
x <- as.irt.pars(c(-1,0,1))


# Create object for three dichotomous (3PL) items and two polytomous (gpcm) items
# without a location parameter (use signature matrix, missing)
dichot <- matrix(c(1.2, .8, .9, 2.3, -1.1, -.2, .24, .19, .13),3,3)
poly <- matrix(c(.64, -1.8, -.73, .45, NA, .88, .06, 1.4, 1.9, 2.6),2,5,byrow=TRUE)
pars <- rbind(cbind(dichot,matrix(NA,3,2)),poly)
cat <- c(2,2,2,4,5)
pm <- as.poly.mod(5, c("drm","gpcm"), list(1:3,4:5))
x <- as.irt.pars(pars, cat=cat, poly.mod=pm)
summary(x)


# Create object for three dichotomous (3PL) items and two polytomous (gpcm) items
# without a location parameter (use signature list, missing)
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
x <- as.irt.pars(pars, cat=cat, poly.mod=pm)
summary(x)


# Create object for three dichotomous (3PL) items, four polytomous items,
# two gpcm items and two nrm items. Include a location parameter for the 
# gpcm items (use signature list, missing)
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
x <- as.irt.pars(pars, cat=cat, poly.mod=pm, location=TRUE)
summary(x, TRUE)


# Create object with two groups (all dichotomous items)
pm <- as.poly.mod(36)
x <- as.irt.pars(KB04$pars, KB04$common, cat=list(rep(2,36),rep(2,36)), 
  list(pm,pm), grp.names=c("form.x","form.y"))
summary(x, TRUE)


# Create object with six groups (all dichotomous items)
pars <- TK07$pars
common <- TK07$common
cat <- list(rep(2,26),rep(2,34),rep(2,37),rep(2,40),rep(2,41),rep(2,43))
pm1 <- as.poly.mod(26)
pm2 <- as.poly.mod(34)
pm3 <- as.poly.mod(37)
pm4 <- as.poly.mod(40)
pm5 <- as.poly.mod(41)
pm6 <- as.poly.mod(43)
pm <- list(pm1, pm2, pm3, pm4, pm5, pm6)
x <- as.irt.pars(pars, common, cat, pm, grp.names=paste("grade",3:8,sep=""))

# Create an object with two groups using mixed-format items and 
# a mixed placement of common items. This example uses the dgn dataset.
pm1=as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group1)
pm2=as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group2)
x=as.irt.pars(dgn$pars,dgn$common,dgn$cat,list(pm1,pm2))
summary(x, TRUE)
}

\keyword{ utilities }
