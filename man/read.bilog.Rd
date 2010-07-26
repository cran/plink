\name{read.bilog}
\alias{read.bilog}
\alias{read.parscale}
\alias{read.multilog}
\alias{read.testfact}
\alias{read.icl}
\alias{read.bmirt}
\alias{read.erm}
\alias{read.ltm}
\title{Import Parameters from IRT Software}
\description{
This function imports item and/or ability parameters from BILOG-MG 3, PARSCALE 4, 
MULTILOG 7, TESTFACT 4, ICL, BMIRT, and ltm.
}

\usage{
read.bilog(file, ability = FALSE, pars.only = TRUE, as.irt.pars = TRUE)

read.parscale(file, ability = FALSE, loc.out = FALSE, pars.only = TRUE, 
  as.irt.pars = TRUE)

read.multilog(file, cat, poly.mod, ability = FALSE, contrast = "dev", 
  drm.3PL = TRUE, loc.out = FALSE, as.irt.pars = TRUE)

read.testfact(file, ability = FALSE, guessing = FALSE, bifactor = FALSE, 
  as.irt.pars = TRUE)

read.icl(file, poly.mod, ability = FALSE, loc.out = FALSE, 
  as.irt.pars = TRUE)

read.bmirt(file, ability = FALSE, sign.adjust = TRUE, loc.out = FALSE, 
  pars.only = TRUE, as.irt.pars = TRUE)

read.erm(x, loc.out = FALSE, as.irt.pars = TRUE)

read.ltm(x, loc.out = FALSE, as.irt.pars = TRUE)
}

\arguments{
  \item{file}{filename of file containing the item or ability parameters}

  \item{ability}{if \code{TRUE}, \code{file} contains ability parameters}

  \item{pars.only}{if \code{TRUE}, only the item/ability parameters will be imported 
     (i.e., any other information like standard errors will be dropped). }

  \item{loc.out}{if \code{TRUE}, the step/threshold parameters will be reformated to 
	be deviations from a location parameter }

  \item{as.irt.pars}{if \code{TRUE}, the parameters will be output as an 
	\code{\linkS4class{irt.pars}} object (this is only applicable to item
	parameters)}

  \item{cat}{vector with the number of response categories for each item. For 
      multiple-choice model items, \code{cat} is the number of response categories 
      plus one (the additional category is for 'do not know')}

  \item{poly.mod}{a \code{\linkS4class{poly.mod}} object.  See the documentation for
      the function \code{\link{as.irt.pars}} for more information on creating this 
      object.}

  \item{contrast}{an object identifying the type of contrast(s) used to estimate the
      various parameters for each item.  See below for more details.}

  \item{drm.3PL}{logical value indicating whether the dichotomous items (if applicable)
      were modeled using the three parameter logistic model (3PL)}

  \item{guessing}{logical value indicating whether a guessing parameter was modeled}

  \item{bifactor}{logical value indicating whether the bifactor model was used to
      estimate the item/ability parameters}

  \item{sign.adjust}{logical value indicating whether the difficulty/step parameters 
  should be multiplied by -1 to make them consistent with common formulations of 
  multidimensional response models}

  \item{x}{output object from one of the following functions in the \code{eRm} package:
      \code{LLTM}, \code{LPCM}, \code{LRSM},\code{PCM}, \code{RM}, or \code{RSM} or one
      of the following functions in the \code{ltm} package: \code{rasch}, \code{ltm}, 
      \code{tpm}, \code{grm}, or \code{gpcm} }
}

\details{
The file extensions for the item parameter and ability files respectively are as follows:
.par and .sco for BILOG-MG, PARSCALE, and MULTILOG, .par and .fsc for TESTFACT, and 
.par and .ss for BMIRT. For \code{ICL}, the file extensions are specified by the user, 
and for \code{ltm}, the name of the output object is specified by the user.

When item parameters are estimated in MULTILOG for models other than the 1PL, 2PL, and GRM,
the program estimates (and returns) contrast parameters.  MULTILOG implements three types
of contrasts: deviation, polynomial, and triangle (see Thissen & Steinberg, 1986 for more
information). A single type of contrast can be used for all parameters (a, b, and c) for all
items or different contrasts can be specified for individual parameters and individual items.
If a single type of contrast is used for all parameters for all items, a character value 
can be specified for the \code{contrast} argument: "dev", "poly", or "tri" for the three
types of contrasts respectively. When different contrasts are used, \code{contrast} should
be a list of length nine. The list elements should be ordered as follows "dev.a","poly.a",
"tri.a","dev.c","poly.c","tri.c","dev.d", "poly.d","tri.d" where the first three elements 
correspond to the various contrasts for the slope parameters, the next three elements 
correspond to the contrasts for the category parameters, and the last three elements 
correspond to the contrasts for the lower asymptote (guessing parameters). There are two
approaches that can be implemented using this list 1) character vectors with the model
names "drm", "grm", "gpcm", "nrm", and "mcm" indicating that the given parameters for
all items associated with the given model should be transformed using the specified contrast.
In instances where a model is not included for a given parameter (for any of the contrasts)
the parameters will be transformed using deviation contrasts.  2) numeric vectors 
identifying the contrasts used for given parameters for given items can be specified.  
It is only necessary to include item numbers for the various parameter/contrast combinations 
when deviation contrasts are not used. See below for examples of how to formulate this argument.
}

\value{
Returns a data.frame or an object of class \code{\linkS4class{irt.pars}} if 
\code{as.irt.pars} = \code{TRUE}.
}

\references{ 
  Hanson, B. A. (2002). 
  IRT command language [Computer Program].
  URL \url{http://www.b-a-h.com/software/irt/icl/}

  Mair, P & Hatzinger, R. (2007).
  Extended Rasch modeling: The eRm package for the application of IRT models in R.
  URL \url{http://www.jstatsoft.org/v20/i09}

  Muraki, E. & Bock, R. D. (2003). 
  PARSCALE 4: IRT item analysis and test scoring for rating scale data [Computer Program].
  Chicago, IL: Scientific Software International. URL \url{http://www.ssicentral.com}

  Rizopoulos, D. (2006). 
  ltm: An R package for latent variable modeling and item response theory analyses. 
  \emph{Journal of Statistical Software}, \emph{17}(5), 1-25. 
  URL \url{http://www.jstatsoft.org/v17/i05/}

  Thissen, D. (2003). 
  MULTILOG 7: Multiple, categorical item analysis and test scoring using item response theory [Computer Program].
  Chicago, IL: Scientific Software International. URL \url{http://www.ssicentral.com}

  Thissen, D. & Steinberg, L. (1986).
  A taxonomy of item response models.
  \emph{Psychometrika}, \emph{51}(4), 567-577.

  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}

  Wood, R., Wilson, D. T., Muraki, E., Schilling, S. G., Gibbons, R., & Bock, R. D. (2003).
  TESTFACT 4: Test scoring, item statistics, and item factor analysis [Computer Program].
  Chicago, IL: Scientific Software International. URL \url{http://www.ssicentral.com}

  Yao, L. (2008). 
  BMIRT: Bayesian multivariate item response theory [Computer Program]. 
  Monterey, CA: CTB/McGraw-Hill.

  Zimowski, M. F., Muraki, E., Mislevy, R. J., & Bock, R. D. (2003). 
  BILOG-MG 3: Multiple-group IRT analysis and test maintenance for binary items [Computer Program].
  Chicago, IL: Scientific Software International. URL \url{http://www.ssicentral.com}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com} }

\note{These functions are currently unable to handle output generated when subtests are used.}

\examples{
# Illustration of how to formulate the contrast argument. Say that we 
# have 20 items where the first 15 are modeled using the 3PL and the 
# last five are modeled using the GPCM.  For the 3PL items, deviation 
# contrasts are commonly used for all of the parameters, whereas with 
# the GPCM items, polynomial contrasts are typically used for the slope 
# parameters and triangle contrasts are used for the category parameters. 
# The contrast argument could be specified as follows

contrast <- vector("list",9)
# Note: the list elements do not need to be named for read.multilog
names(contrast) <- c("dev.a","poly.a","tri.a","dev.c","poly.c","tri.c",
  "dev.d", "poly.d","tri.d") 
contrast$poly.a <- 16:20
contrast$tri.c <- 16:20

# The object could alternatively be formatted as follows
contrast <- vector("list",9)
names(contrast) <- c("dev.a","poly.a","tri.a","dev.c","poly.c","tri.c",
  "dev.d","poly.d","tri.d") 
contrast$dev.a <- 1:15
contrast$poly.a <- 16:20
contrast$dev.c <- 1:15
contrast$tri.c <- 16:20
contrast$dev.d <- 1:15
}

\keyword{ utilities }
