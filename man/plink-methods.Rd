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
This function conducts separate calibration of IRT single-format or 
mixed-format item parameters for multiple groups using the Mean/Mean, 
Mean/Sigma, Haebara, and Stocking-Lord methods. The function includes
symmetric and non-symmetric optimization and chain-linked rescaling
of item and ability parameters.
}

\usage{
plink(x, common, rescale, ability, weights, startvals, score=1, base.grp=1, 
	symmetric=TRUE, grp.names=NULL, mn=c(FALSE,FALSE), ...)

## S4 method for signature 'list', 'matrix'
plink(x, common, rescale, ability, weights, startvals, score, base.grp, 
	symmetric, grp.names, mn, ...)

## S4 method for signature 'list', 'data.frame'
plink(x, common, rescale, ability, weights, startvals, score, base.grp, 
	symmetric, grp.names, mn, ...)

## S4 method for signature 'list', 'list'
plink(x, common, rescale, ability, weights, startvals, score, base.grp, 
	symmetric, grp.names, mn, ...)

## S4 method for signature 'irt.pars', 'ANY'
plink(x, common, rescale, ability, weights, startvals, score, base.grp, 
	symmetric, grp.names, mn, ...)
}

\arguments{
  \item{x}{an object of class \code{\linkS4class{irt.pars}} with multiple groups or
    a list of \code{irt.pars} and/or \code{\linkS4class{sep.pars}} objects. }

  \item{common}{matrix or list of common items. See below for more details. }

  \item{rescale}{if missing (default), the parameters in \code{x} will not be 
    transformed to the base scale. To transform the parameters use \code{"MM","MS","HB","SL"}
    for the mean/mean, mean/sigma, Haebara, and Stocking-Lord linking constants
    respectively.}

  \item{ability}{list of theta values with length equal to the number of groups. 
    If supplied, these values will be transformed to the base scale using the
    constants identified in \code{rescale} or the Haebara constants if \code{rescale}
    is missing.}

  \item{weights}{list containing information about the theta values and weights
    to use for the characteristic curve methods. See below for more details. }

  \item{score}{if \code{score} = \code{1}, score responses for the Stocking-Lord
    method with zero for the lowest category and k-1 for the highest, k, category
    for each item. If \code{score} = \code{2}, score responses with one for the
    lowest category and k for the highest, k, category for each item. A vector of
    scores for each response category can be supplied, but this is only recommended
    for advanced users.}

  \item{startvals}{vector of length two of starting values for the slope and 
    intercept respectively for the characteristic curve methods. The default
    values (if \code{startvals} is missing) are 1 and 0.}

  \item{base.grp}{integer identifying the group for the base scale}

  \item{symmetric}{if \code{TRUE} use symmetric minimization for the characteristic 
    curve methods. See Kim and Lee (2006) for more information}

  \item{grp.names}{character vector of group names }

  \item{mn}{logical vector of length two. If the first element is \code{FALSE}, do use
    nominal response model and multiple-choice model parameters when estimating
    mean/mean and mean/sigma linking constants. If the second element is \code{FALSE},
    do not use NRM and MCM parameters with the Stocking-Lord calibration.}

  \item{...}{further arguments passed to or from other methods }
}

\details{
If \code{x} contains only two elements, \code{common} should be a matrix.  If \code{x} 
contains more than two elements, \code{common} should be a list.  In any of the \code{common} 
matrices the first column identifies the common items for the first group of two adjacent 
list elements in \code{x}.  The second column in \code{common} identifies the corresponding 
set of common items from the next list element in \code{x}. For example, if \code{x} 
contains only two list elements, a single set of common items links them together.  If 
item 4 in group one (row 4 in slot \code{pars}) is the same as item 6 in group two, the 
first row of \code{common} would be \code{"4,6"}.

\code{weights} can be a list or a list of lists.  The purpose of this object is to specify
the theta values to integrate over in the characteristic curve methods as well as any 
weights associated with the theta values.  See Kim and Lee (2006) or Kolen and Brennan (2004)
for more information of these weights.The function \code{\link{as.weight}} can be used
to facilitate the creation of this object.  If \code{weights} is missing, the default is
to use equally spaced theta values ranging from -4 to 4 with an increment of 0.05 and 
theta weights equal to one for all theta values.

To better understand the elements of \code{weights}, let us assume for a moment that \code{x}
has parameters for only two groups.  In this instance, \code{weights} would be a single 
list with length two. The first element should be an n x 2 matrix of theta values 
corresponding to group 1 and group 2 respectively.  The second list element should be an 
n x 2 matrix of weights corresponding the the theta values.  If \code{x} contains multiple 
groups, a single \code{weights} object can be supplied, and the same set of thetas and weights
will be used for all adjacent groups. However, if \code{Wj} is the \code{weights} object
for adjacent groups \code{j}, \code{weights} can be a list with multiple \code{W} objects
(i.e. a separate list of theta values and theta weights for each adjacent group in \code{x}).

In generally, including a different \code{W} for each adjacent group is unnecessary, but
the most likely scenario for its use is when the intent is to integrate over empirical
theta values instead of equally spaced or random values.
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
  Returns an object of class \code{\linkS4class{link}}
}

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

\author{ Jonathan Weeks \email{weeksjp@gmail.com} }

\examples{
# Create irt.pars object with two groups (all dichotomous items)
# rescale the item parameters using the Stocking-Lord linking constants
pm <- as.poly.mod(36)
x <- as.irt.pars(KB04$pars, KB04$common, cat=list(rep(2,36),rep(2,36)), 
  poly.mod=list(pm,pm))
out <- plink(x, rescale="SL")
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
x <- as.irt.pars(pars, common, cat, pm, grp.names=paste("grade",3:8,sep=""))
out <- plink(x)
summary(out)
constants <- link.con(out) # Extract linking constants

# Create an irt.pars object and a sep.pars object for two groups of
# nominal response model items. Compare symmetric and non-symmetric minimization
# Note: This example may take a minute or two to run
pm <- as.poly.mod(60, "nrm", 1:60)
pars1 <- as.irt.pars(act.nrm$yr97, cat=rep(5,60), poly.mod=pm)
pars2 <- sep.pars(act.nrm$yr98, cat=rep(5,60), poly.mod=pm)
out <- plink(list(pars1, pars2), matrix(1:60,60,2))
out1 <- plink(list(pars1, pars2), matrix(1:60,60,2), symmetric=FALSE)
summary(out, descrip=TRUE)
summary(out1, descrip=TRUE)


# Compute linking constants for two groups with multiple-choice model
# item parameters. Rescale theta values and item parameters using
# the Haebara linking constants
# Note: This example may take a minute or two to run
theta <- rnorm(100) # In practice, estimated theta values would be used
pm <- as.poly.mod(60, "mcm", 1:60)
x <- as.irt.pars(act.mcm, common=matrix(1:60,60,2), cat=list(rep(6,60),
  rep(6,60)), poly.mod=list(pm,pm))
out <- plink(x, ability=list(theta,theta), rescale="HB", symmetric=FALSE)
pars.out <- link.pars(out)
ability.out <- link.ability(out)
summary(out, descrip=TRUE)


# Compute linking constants for two groups using mixed-format items and 
# a mixed placement of common items. Compare calibrations with the
# inclusion or exclusion of NRM items. This example uses the dgn dataset.
pm1=as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group1)
pm2=as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group2)
x=as.irt.pars(dgn$pars,dgn$common,dgn$cat,list(pm1,pm2))
out <- plink(x) # Run with NRM items included in the calibration
out1 <- plink(x,mn=c(TRUE,TRUE)) # Run with NRM items excluded from the calibration
summary(out)
summary(out1)
}

\keyword{methods}
