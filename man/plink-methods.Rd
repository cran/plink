\name{plink-methods}
\docType{methods}
\alias{plink}
\alias{plink-methods}
\alias{plink,list,matrix-method}
\alias{plink,list,data.frame-method}
\alias{plink,list,list-method}
\alias{plink,irt.pars,ANY-method}
\title{Separate Calibration}
\description{
This function conducts separate calibration of unidimensional or multidimensional
IRT single-format or mixed-format item parameters for multiple groups.  The  
unidimensional methods include the Mean/Mean, Mean/Sigma, Haebara, and Stocking-Lord 
methods.  The multidimensional methods include a least squares approach (Li & Lissitz,
2000; Min, 2003; Reckase & Martineau, 2004) and extensions of the Haebara and 
Stocking-Lord method using a single dilation parameter (Li & Lissitz, 2000), multiple 
dilation parameters (Min, 2003), or the Oshima, Davey, & Lee (2000) approach.  The 
methods allow for symmetric and non-symmetric optimization and chain-linked rescaling 
of item and ability parameters.
}

\usage{
plink(x, common, rescale, ability, method, weights.t, weights.f, 
  startvals, exclude, score = 1, base.grp = 1, symmetric = FALSE, 
  rescale.com = TRUE, grp.names = NULL, dilation = "LL", 
  dim.order = NULL, ...)

## S4 method for signature 'list', 'matrix'
plink(x, common, rescale, ability, method, weights.t, weights.f, 
  startvals, exclude, score, base.grp, symmetric, rescale.com, 
  grp.names, dilation, dim.order, ...)

## S4 method for signature 'list', 'data.frame'
plink(x, common, rescale, ability, method, weights.t, weights.f, 
  startvals, exclude, score, base.grp, symmetric, rescale.com, 
  grp.names, dilation, dim.order, ...)

## S4 method for signature 'list', 'list'
plink(x, common, rescale, ability, method, weights.t, weights.f, 
  startvals, exclude, score, base.grp, symmetric, rescale.com, 
  grp.names, dilation, dim.order, ...)

## S4 method for signature 'irt.pars', 'ANY'
plink(x, common, rescale, ability, method, weights.t, weights.f, 
  startvals, exclude, score, base.grp, symmetric, rescale.com, 
  grp.names, dilation, dim.order, ...)
}

\arguments{
  \item{x}{an object of class \code{\linkS4class{irt.pars}} with multiple groups or
    a list of \code{irt.pars} and/or \code{\linkS4class{sep.pars}} objects. }

  \item{common}{matrix or list of common items. See below for more details. }

  \item{rescale}{if missing (default), the parameters in \code{x} will not be 
    transformed to the base scale. To transform the parameters use 
    \code{"MM","MS","HB","SL","LS"} for the Mean/Mean, Mean/Sigma, Haebara, 
    Stocking-Lord, and least squares linking constants respectively.}

  \item{ability}{list of theta values with length equal to the number of groups. 
    If supplied, these values will be transformed to the base scale using the
    constants identified in \code{rescale} or the Haebara constants if \code{rescale}
    is missing.}

  \item{method}{character vector identifying the linking methods to use. Values can
    include \code{"MM","MS","HB","SL","LS"} for the Mean/Mean, Mean/Sigma, Haebara, 
    Stocking-Lord, and least squares linking constants respectively, or if
    missing, constants will be estimated for all methods.}

  \item{weights.t}{list containing information about the theta values and weights
    on the \emph{To} scale for use with the characteristic curve methods. 
    See below for more details. }

  \item{weights.f}{list containing information about the theta values and weights
    on the \emph{From} scale for use with the characteristic curve methods. This
    argument will be ignored if \code{symmetric}=FALSE. See below for more details. }

  \item{startvals}{vector or list of slope and intercept starting values for the 
    characteristic curve methods or character value(s) of "MM", "MS", or "LS" to identify 
    that values from the Mean/Mean, Mean/Sigma, or least squares method should be used
    for the starting values. See below for more details.}

  \item{exclude}{character vector or list identifying common items that should be
    excluded when estimating the linking constants. See below for more details.}

  \item{score}{if \code{score} = 1, score responses for the Stocking-Lord
    method with zero for the lowest category and k-1 for the highest, k, category
    for each item.  If \code{score} = 2, score responses with one for the
    lowest category and k for the highest, k, category for each item.  A vector or 
    list of scoring weights for each response category can be supplied, but this is 
    only recommended for advanced users.}

  \item{base.grp}{integer identifying the group for the base scale}

  \item{symmetric}{if \code{TRUE} use symmetric minimization for the characteristic 
    curve methods.  See Kim and Lee (2006) for more information}

  \item{rescale.com}{if \code{TRUE} rescale the common item parameters using the 
    estimated linking constants; otherwise, insert the non-transformed common item
    parameters into the set of unique transformed item parameters}

  \item{grp.names}{character vector of group names }

  \item{dilation}{character value identifying whether the Oshima, Davey, & Lee (2000)
    approach \code{"ODL"}, the Li and Lissitz (2000) approach \code{"LL"}, or the Min (2003)
    approach \code{"MIN"} should be used to estimate the linking constants for the 
    multidimensional Haebara and Stocking-Lord methods.}

  \item{dim.order}{matrix for identifying the ordering of factors across groups.  See 
    below for details.}

  \item{...}{further arguments passed to or from other methods. See below for details.}
}

\details{
If \code{x} contains only two elements, \code{common} should be a matrix.  If \code{x} 
contains more than two elements, \code{common} should be a list.  In any of the \code{common} 
matrices the first column identifies the common items for the first group of two adjacent 
list elements in \code{x}.  The second column in \code{common} identifies the corresponding 
set of common items from the next list element in \code{x}. For example, if \code{x} 
contains only two list elements, a single set of common items links them together.  If 
item 4 in group one (row 4 in slot \code{pars}) is the same as item 6 in group two, the 
first row of \code{common} would be \code{(4,6)}.

\code{startvals} can be a vector or list of starting values for the slope(s) and intercept(s). 
For unidimensional methods, when there are only two groups, this argument should be a vector of 
length of two with the first value for the slope and the second value for the intercept or a
character value equal to "MM" or"MS". When there are more than two groups a vector of starting 
values or a character value can be supplied (the same numeric values, if a vector is supplied,
will be used for all pairs, or the corresponding MM/MS values for each pair of tests will be 
used) or a list of vectors/character values can be supplied with the number of list elements
equal to the number of groups minus one. For the multidimensional methods, the same general 
structure applies (a vector or character value for a single group or a vector, character value 
or list for multiple groups); however, the length of the vector will vary depending on the 
dilation approach used. If \code{dilation} is "ODL", the first m*m values in \code{startvals}, 
for m dimensions, should correspond to the values in the transformation matrix (starting with 
the value in the upper-left corner, then the next value in the column, ..., then the first value 
in the next column, etc.). The remaining m values should be for the translation vector. If 
\code{dilation} is "LL", the first value will be the slope parameter and the remaining m values 
will be for the translation vector. If \code{dilation} is "MIN", the first m values are the 
slopes for each dimension and the remaining m values are for the translation vector. 

\code{weights.t} can be a list or a list of lists.  The purpose of this object is to specify
the theta values on the \emph{To} scale to integrate over in the characteristic curve methods 
as well as any weights associated with the theta values.  See Kim and Lee (2006) or Kolen 
and Brennan (2004) for more information of these weights.  The function \code{\link{as.weight}} 
can be used to facilitate the creation of this object.  If \code{weights.t} is missing, the 
default--in the unidimensional case--is to use equally spaced theta values ranging from -4 to 4 
with an increment of 0.05 and theta weights equal to one for all theta values. In the
multidimensional case, the default is to use 1000 randomly sampled values from a multivariate
normal distribution with correlations equal to 0.6 for all dimensions.  The theta weights are
equal to the normal distribution weights at these points. 

To better understand the elements of \code{weights.t}, let us assume for a moment that \code{x}
has parameters for only two groups and that we are using non-symmetric linking.  In this instance, 
\code{weights.t} would be a single list with length two.  The first element should be a vector
of theta values corresponding to points on the \emph{To} scale.  The second list element should 
be a vector of weights corresponding the theta values.  If \code{x} contains multiple 
groups, a single \code{weights.t} object can be supplied, and the same set of thetas and weights
will be used for all adjacent groups.  However, a separate list of theta values and theta weights 
for each adjacent group in \code{x} can be supplied.

The specification of \code{weights.f} is the same as that for \code{weights.t}, although the 
theta values and weights for this object correspond to theta values on the \emph{From} scale.
This argument is only used when \code{symmetric}=TRUE. If \code{weights.f} is missing and 
\code{symmetric}=TRUE, the same theta values and weights used for \code{weights.t} will be
used for \code{weights.f}.

In general, all of the common items identified in \code{x} or \code{common} will be used
in estimating linking constants; however, there are instances where there is a need to exclude
certain common items (e.g., NRM or MCM items or items exhibiting parameter drift). Instead of
creating a new matrix or list of matrices for \code{common}, the \code{exclude} argument can
be used. \code{exclude} can be specified as a character vector or a list. In the former case,
a vector of model names (i.e., "drm", "grm", "gpcm", "nrm", "mcm") would be supplied, indicating
that any common item on any test associated with the given model(s) would be excluded from the
set of items used to estimate the linking constants. If the argument is specified as a list, 
\code{exclude} should have as many elements as groups in \code{x}. Each list element can include
model names and/or item numbers corresponding to the common items that should be excluded for
the given group. If no items need to be excluded for a given group, the list element should be
NULL or NA. For example, say we have two groups and we would like to exclude the NRM items and
item 23 from the first group, we would specify \code{exclude} as \code{exclude <- list(c("nrm",23),NA)}.
Notice that the item number corresponding item 23 in group 2 does not need to be specified.

The argument \code{dim.order} is a k x r matrix for k groups and r unique dimensions across
groups. This object identifies the common dimensions across groups. The elements in the 
matrix should correspond to the dimension (i.e., the column in the matrix of slope parameters)
for a given group. For example, say there are four unique dimensions across two groups,
each group only measures three dimensions, and there are only two common dimensions. We might
specify a matrix as follows dim.order <- matrix(c(1:3,NA,NA,1:3),2,4). In words, this means
that dimensions 2 and 3 in the first group correspond to dimensions 1 and 2 in the second
group respectively. If no \code{dim.order} is specified, it is assumed that all of the 
dimensions are common, or in instances with different numbers of factors, that the first
m dimensions for each group are common and the remaining r-m dimensions for a given group
are unique.

For the characteristic curve methods, response probabilities are computed using the functions
\code{drm}, \code{grm}, \code{gpcm}, \code{nrm}, and \code{mcm} for the associated models
respectively. Various arguments from these functions can be passed to \code{plink}. Specifically,
the argument \code{incorrect} can be passed to \code{drm} and \code{catprob} can be passed to
\code{grm}. In the functions \code{drm}, \code{grm}, and \code{gpcm} there is an argument \code{D}
for the value of a scaling constant. In \code{plink}, a single argument \code{D} can be passed
that will be applied to all applicable models, or arguments \code{D.drm}, \code{D.grm}, and 
\code{D.gpcm} can be specified for each model respectively. If an argument is specified for \code{D}
and, say \code{D.drm}, the values for \code{D.grm} and \code{D.gpcm} (if applicable) will be
set equal to \code{D}. If only \code{D.drm} is specified, the values for \code{D.grm} and 
\code{D.gpcm} (if applicable) will be set to 1.
}

\section{Methods}{
  \describe{

    \item{x = "list", common = "matrix"}{This method is used when \code{x} contains only
      two list elements. If either of the list elements is of class \code{irt.pars}, they
      can include multiple groups. \code{common} is the matrix of common items between
      the two groups in \code{x}.  See details for more information on \code{common}.}

    \item{x = "list", common = "data.frame"}{See the method for signature \code{x}="list", 
      \code{common}="matrix".}

    \item{x = "list", common = "list"}{This method is used when \code{x} includes two or
      more list elements. When \code{x} has length two, \code{common} (although a single
      matrix) should be a list with length one. If \code{x} has more than two list elements
      \code{common} identifies the common items between adjacent list elements.  If objects
      of class \code{irt.pars} are included with multiple groups, \code{common} should
      identify the common items between the first or last group in the \code{irt.pars} object,
      depending on its location in \code{x}, and the adjacent list element(s) in \code{x}.
      For example, if \code{x} has three elements: an \code{irt.pars} object with one group,
      an \code{irt.pars} object with four groups, and a \code{sep.pars} object, \code{common}
      will be a list with length two. The first element in \code{common} is a matrix 
      identifying the common items between the items in the first \code{irt.pars} object
      and the first group in the second \code{irt.pars} object. The second element in 
      \code{common} should identify the common items between the fourth group in the 
      second \code{irt.pars} object and the items in the \code{sep.pars} object.
      }

    \item{x = "irt.pars", common = "ANY"}{This method is intended for an \code{irt.pars}
      object with multiple groups.}
  }
}

\value{
  Returns an object of class \code{\linkS4class{link}}. The labels for the linking constants are 
  specified in the following manner "group1/group2", meaning the group1 parameters were transformed 
  to the group2 test. The base group is indicated by an asterisk.
}

\references{ 
  Haebara, T. (1980). 
  Equating logistic ability scales by a weighted least squares method. 
  \emph{Japanese Psychological Research}, \emph{22}(3), 144-149.

  Kim, S. & Lee, W.-C. (2006). 
  An Extension of Four IRT Linking Methods for Mixed-Format Tests. 
  \emph{Journal of Educational Measurement}, \emph{43}(1), 53-76.

  Kolen, M. J. & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking} (2nd ed.).
  New York: Springer

  Li, Y. H., & Lissitz, R. W. (2000). 
  An evaluation of the accuracy of multidimensional IRT linking. 
  \emph{Applied Psychological Measurement}, \emph{24}(2), 115-138.

  Loyd, B. H. & Hoover, H. D. (1980). 
  Vertical Equating Using the Rasch Model. 
  \emph{Journal of Educational Measurement}, \emph{17}(3), 179-193.

  Marco, G. L. (1977). 
  Item Characteristic Curve Solutions to Three Intractable Testing Problems.
  \emph{Journal of Educational Measurement}, \emph{14}(2), 139-160.

  Min, K. -S. (2007). 
  Evaluation of linking methods for multidimensional IRT calibrations. 
  \emph{Asia Pacific Education Review}, \emph{8}(1), 41-55. 

  Oshima, T. C., Davey, T., & Lee, K. (2000). 
  Multidimensional linking: Four practical approaches. 
  \emph{Journal of Educational Measurement}, \emph{37}(4), 357-373.

  Reckase, M. D., & Martineau, J. A. (2004). 
  \emph{The vertical scaling of science achievement tests}.
  Research report for the Center for Education and National Research Council.

  Stocking, M. L. & Lord, F. M. (1983). 
  Developing a common metric in item response theory. 
  \emph{Applied Psychological Measurement}, \emph{7}(2), 201-210.

  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\note{The translation vector for the multidimensional Stocking-Lord
method may converge to odd values. The least squares method has been
shown to produce more accurate constants than the characteristic curve
methods (in less time); however, if use of a characteristic curve 
approach is desired, it is recommended that the Haebara method be used 
and/or the relative tolerance for the optimization be lowered to 0.001 
using the argument \code{control=list(rel.tol=0.001)}.}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\examples{
###### Unidimensional Examples ######
# Create irt.pars object with two groups (all dichotomous items),
# rescale the item parameters using the Mean/Sigma linking constants,
# and exclude item 27 from the common item set
pm <- as.poly.mod(36)
x <- as.irt.pars(KB04$pars, KB04$common, cat=list(rep(2,36),rep(2,36)), 
  poly.mod=list(pm,pm), exclude=list(27,NA))
out <- plink(x, rescale="MS", base.grp=2, D=1.7)
summary(out, descrip=TRUE)
pars.out <- link.pars(out)


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
x <- as.irt.pars(pars, common, cat, pm, 
  grp.names=paste("grade",3:8,sep=""))
out <- plink(x)
summary(out)
constants <- link.con(out) # Extract linking constants

# Create an irt.pars object and a sep.pars object for two groups of
# nominal response model items. Compare non-symmetric and symmetric 
# minimization. Note: This example may take a minute or two to run
\dontrun{
pm <- as.poly.mod(60, "nrm", 1:60)
pars1 <- as.irt.pars(act.nrm$yr97, cat=rep(5,60), poly.mod=pm)
pars2 <- sep.pars(act.nrm$yr98, cat=rep(5,60), poly.mod=pm)
out <- plink(list(pars1, pars2), matrix(1:60,60,2))
out1 <- plink(list(pars1, pars2), matrix(1:60,60,2), symmetric=TRUE)
summary(out, descrip=TRUE)
summary(out1, descrip=TRUE)
}

# Compute linking constants for two groups with multiple-choice model
# item parameters. Rescale theta values and item parameters using
# the Haebara linking constants
# Note: This example may take a minute or two to run
\dontrun{
theta <- rnorm(100) # In practice, estimated theta values would be used
pm <- as.poly.mod(60, "mcm", 1:60)
x <- as.irt.pars(act.mcm, common=matrix(1:60,60,2), cat=list(rep(6,60),
  rep(6,60)), poly.mod=list(pm,pm))
out <- plink(x, ability=list(theta,theta), rescale="HB")
pars.out <- link.pars(out)
ability.out <- link.ability(out)
summary(out, descrip=TRUE)
}

# Compute linking constants for two groups using mixed-format items and 
# a mixed placement of common items. Compare calibrations with the
# inclusion or exclusion of NRM items. This example uses the dgn dataset.
pm1 <- as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group1)
pm2 <- as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group2)
x <- as.irt.pars(dgn$pars,dgn$common,dgn$cat,list(pm1,pm2))
# Run with the NRM common items included
out <- plink(x) 
# Run with the NRM common items excluded
out1 <- plink(x,exclude="nrm") 
summary(out)
summary(out1)


# Compute linking constants for six groups using mixed-format items and 
# a mixed placement of common items. This example uses the reading dataset.
# See the information on the dataset for an interpretation of the output.
pm1 <- as.poly.mod(41, c("drm", "gpcm"), reading$items[[1]])
pm2 <- as.poly.mod(70, c("drm", "gpcm"), reading$items[[2]])
pm3 <- as.poly.mod(70, c("drm", "gpcm"), reading$items[[3]])
pm4 <- as.poly.mod(70, c("drm", "gpcm"), reading$items[[4]])
pm5 <- as.poly.mod(72, c("drm", "gpcm"), reading$items[[5]])
pm6 <- as.poly.mod(71, c("drm", "gpcm"), reading$items[[6]])
pm <- list(pm1, pm2, pm3, pm4, pm5, pm6)
x <- as.irt.pars(reading$pars, reading$common, reading$cat, pm, base.grp=4)
out <- plink(x)
summary(out)


###### Multidimensional Examples ######
# Reckase Chapter 9
pm <- as.poly.mod(80, "drm", list(1:80))
x <- as.irt.pars(reckase9$pars, reckase9$common, 
  list(rep(2,80),rep(2,80)), list(pm,pm), dimensions=2)
# Compute constants using the least squares method and 
# the Min (2003) dilation. Rescale the item parameters 
# using the LS method
out <- plink(x, dilation="MIN", rescale="LS")
summary(out, descrip=TRUE)
# Extract the rescaled item parameters
pars.out <- link.pars(out)

# Compute constants using Oshima, Davey, Lee (2000)
# Display output and descriptives
out <- plink(x, dilation="ODL")
summary(out, descrip=TRUE)

# Compute constants using Li & Lissitz (2000)
# Rescale item parameters and ability estimates 
# using the "HB" and "LS" methods.
# Display the optimization iterations
ability <- matrix(rnorm(40),20,2)
out <- plink(x, method=c("HB","LS"), dilation="LL", 
   rescale="HB", ability=list(ability,ability),
   control=list(trace=1,rel.tol=0.001))
summary(out)
# Extract rescaled ability estimates
ability.out <- out$ability


# Compute linking constants for two groups using mixed-format items 
# modeled with the M3PL and MGPCM. Only compute constants using the 
# MIN dilation.
pm <- as.poly.mod(60,c("drm","gpcm"), list(c(1:60)[md.mixed$cat==2],
  c(1:60)[md.mixed$cat>2]))
x <- as.irt.pars(md.mixed$pars, matrix(1:60,60,2), 
  list(md.mixed$cat, md.mixed$cat), 
  list(pm, pm), dimensions=4, grp.names=c("Form.X","Form.Y"))
out <- plink(x,dilation="MIN")
summary(out, descrip=TRUE)


# Illustration of construct shift and how to specify common dimensions
pm <- as.poly.mod(80, "drm", list(1:80))
pars <- cbind(round(runif(80),2),reckase9$pars[[1]])
x <- as.irt.pars(list(pars,reckase9$pars[[2]]), reckase9$common, 
list(rep(2,80),rep(2,80)), list(pm,pm), dimensions=c(3,2))
dim.order <- matrix(c(1,2,3,NA,1,2),2,3,byrow=TRUE)
out <- plink(x, dilation="ODL", dim.order=dim.order, rescale="LS")
summary(out)
pars.out <- link.pars(out)

}

\keyword{methods}
