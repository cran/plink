\name{equate-methods}
\docType{methods}
\alias{equate}
\alias{equate-methods}
\alias{equate,list-method}
\alias{equate,irt.pars-method}
\title{IRT True Score and Observed Score Equating}
\description{
This function conducts IRT true score and observed score equating for unidimensional 
single-format or mixed-format item parameters for two or more groups.  This function
supports all item response models available in \code{plink} with the exception of 
the multiple-choice model.
}

\usage{
equate(x, method="TSE", true.scores, ts.low, base.grp=1, score=1, 
  startval, weights1, weights2, syn.weights, ...) 

## S4 method for signature 'list'
equate(x, method, true.scores, ts.low, base.grp, score, startval, 
  weights1, weights2, syn.weights, ...)

## S4 method for signature 'irt.pars', 'ANY'
equate(x, method, true.scores, ts.low, base.grp, score, startval, 
  weights1, weights2, syn.weights, ...)
}

\arguments{
  \item{x}{an object of class \code{\linkS4class{irt.pars}} with two or more groups
    or the output from \code{plink} containing rescaled item parameters. }

  \item{method}{character vector identifying the equating method(s) to use. 
    Values can include \code{"TSE"} and \code{"OSE"}.}

  \item{true.scores}{numeric vector of true score values to be equated}

  \item{ts.low}{logical value. If TRUE, interpolate values for the equated 
    true scores in the range of observed scores from one to the value below 
    the lowest estimated true score (a rounded sum of guessing parameters)}

  \item{base.grp}{integer identifying the group for the base scale}

  \item{score}{if \code{score} = 1, score responses for the true-score equating
    method with zero for the lowest category and k-1 for the highest, k, category
    for each item.  If \code{score} = 2, score responses with one for the
    lowest category and k for the highest, k, category for each item.  A vector or
    list of scoring weights for each response category can be supplied, but this is 
    only recommended for advanced users.}

  \item{startval}{integer starting value for the first value of \code{true.score}}

  \item{weights1}{list containing information about the theta values and weights
    to be used in the observed score equating for population 1. See below for more 
    details. }

  \item{weights2}{list containing information about the theta values and weights
    to be used in the observed score equating for population 2. See below for more 
    details. }

  \item{syn.weights}{vector of length two or a list containing vectors of length two
    with synthetic population weights to be used for each pair of tests for populations 
    1 and 2 respectively. If missing, weights of 0.5 will be used for both populations
    for all groups.  If \code{syn.weights} is a list, there should be k-1 elements for
    k groups.}

  \item{...}{further arguments passed to or from other methods. See below for details.}
}

\details{
\code{weights1} can be a list or a list of lists.  The purpose of this object is to specify
the theta values for population 1 to integrate over in the observed score equating as well as 
any weights associated with the theta values.  The function \code{\link{as.weight}} 
can be used to facilitate the creation of this object.  If \code{weights1} is missing, the 
default is to use equally spaced theta values ranging from -4 to 4 with an increment of 0.05 
and normal density weights for all groups.

To better understand the elements of \code{weights1}, let us assume for a moment that \code{x}
has parameters for only two groups.  In this instance, \code{weights1} would be a single list 
with length two.  The first element should be a vector of theta values corresponding to points 
on the base scale.  The second list element should be a vector of weights corresponding the 
theta values.  If \code{x} contains more than two groups, a single \code{weights1} object can 
be supplied, and the same set of thetas and weights will be used for all adjacent groups.  
However, a separate list of theta values and weights for each adjacent group in \code{x} can be 
supplied.

The specification of \code{weights2} is the same as that for \code{weights1}, although the 
theta values and weights for this object correspond to theta values for population 2.
This argument is only used when the synthetic weight associated with population 2 is greater
than zero. If \code{weights2} is missing, the same theta values and weights used for 
\code{weights1} will be used for \code{weights2}.

For both equating methods, response probabilities are computed using the functions \code{drm}, 
\code{grm}, \code{gpcm}, and \code{nrm} for the associated models respectively. Various 
arguments from these functions can be passed to \code{equate}. Specifically, the argument 
\code{incorrect} can be passed to \code{drm} and \code{catprob} can be passed to \code{grm}. 
In the functions \code{drm}, \code{grm}, and \code{gpcm} there is an argument \code{D}
for the value of a scaling constant. In \code{plink}, a single argument \code{D} can be passed
that will be applied to all applicable models, or arguments \code{D.drm}, \code{D.grm}, and 
\code{D.gpcm} can be specified for each model respectively. If an argument is specified for \code{D}
and, say \code{D.drm}, the values for \code{D.grm} and \code{D.gpcm} (if applicable) will be
set equal to \code{D}. If only \code{D.drm} is specified, the values for \code{D.grm} and 
\code{D.gpcm} (if applicable) will be set to 1.
}

\value{
  Returns a matrix of equated true scores or a list of equated observed scores with associated
marginal distributions or a list combining these two objects}

\references{ 
  Kolen, M. J. (1981). 
  Comparison of traditional and item response theory methods for equating tests. 
  \emph{Journal of Educational Measurement}, \emph{18}(1), 1-11.

  Kolen, M. J. & Brennan, R. L. (2004)
  \emph{Test Equating, Scaling, and Linking} (2nd ed.).
  New York: Springer

  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\examples{
# IRT true score and observed score examples from 
# Kolen & Brennan (2004, ch. 6)
pm <- as.poly.mod(36)
x <- as.irt.pars(KB04$pars, KB04$common, exclude=list(27,NA),
  cat=list(rep(2,36),rep(2,36)), poly.mod=list(pm,pm))
out <- plink(x, rescale="MS", base.grp=2, D=1.7)

# Create the quadrature points and weights
wt <- as.weight(
  theta=c(-5.2086,-4.163,-3.1175,-2.072,-1.0269,0.0184,
    1.0635,2.109,3.1546,4.2001),
  weight=c(0.000101,0.00276,0.03021,0.142,0.3149,0.3158,
    0.1542,0.03596,0.003925,0.000186))

# Conduct the equating
equate(out,method=c("TSE","OSE"),weights1=wt, synth.weights=c(1,0),D=1.7)

# Conduct true score equating for specific true scores
equate(out, true.scores=7:15, ts.low=FALSE, D=1.7)


# Observed score equating for mixed-format tests
pm1 <- as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group1)
pm2 <- as.poly.mod(55,c("drm","gpcm","nrm"),dgn$items$group2)
x <- as.irt.pars(dgn$pars,dgn$common,dgn$cat,list(pm1,pm2))
out <- plink(x, rescale="HB") 
OSE <- equate(out, method="OSE", score=2)

# Display the equated scores
OSE[[1]]


# Multiple group equating
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
out <- plink(x, rescale="SL")


# True score equating
equate(out)

# True score equating with the base group changed to 3
equate(out, base.grp=3)

# Observed score equating (These data are for non-equivalent groups, but
# this example is included to illustrate the multigroup capabilities)
OSE <- equate(out, method="OSE", base.grp=3)

# Display the equated scores for each group
OSE[[1]]

}

\keyword{methods}
