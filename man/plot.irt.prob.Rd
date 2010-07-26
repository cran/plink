\name{plot.irt.prob}
\alias{plot.irt.prob}
\alias{plot.sep.pars}
\alias{plot.irt.pars}
\alias{plot.list}
\title{Plot Item Stuff}
\description{
This function plots item response curves/surfaces using the \code{lattice} 
package in addition to vector plots, contour plots, and level plots for 
multidimensional items, and comparison plots for examining parameter drift. 
}

\usage{
## S3 method for signature 'irt.pars'
\method{plot}{irt.pars}(x, y, ..., type, separate, combine, items, item.names, 
    item.nums, panels, drift, groups, grp.names, sep.mod, 
    drift.sd, save.hist)

## S3 method for signature 'irt.prob'
\method{plot}{irt.prob}(x, y, ..., type, separate, combine, items, item.names, 
    item.nums, panels, save.hist)

## S3 method for signature 'sep.pars'
\method{plot}{sep.pars}(x, y, ..., type, separate, combine, items, item.names, 
    item.nums, panels, save.hist) 

## S3 method for signature 'list'
\method{plot}{irt.pars}(x, y, ..., type, separate, combine, items, item.names, 
    item.nums, panels, drift, groups, grp.names, sep.mod, 
    drift.sd, save.hist)
}

\arguments{
  \item{x}{object of class \code{irt.prob}, \code{irt.pars}, or \code{sep.pars}.  For the
  later two classes, probabilities are computed using the function \code{mixed} before
  plotting the curves.}

  \item{y}{this is an argument in the generic plot function, but is not used for these methods}

  \item{\dots}{further arguments passed to or from other methods}

  \item{type}{type of plot to produce (for multidimensional models). Values can include
      "wireframe", "contourplot", "levelplot", "vectorplot1", "vectorplot2", "vectorplot3". 
      The default if \code{type} is missing is "wireframe".  The first three types of plots 
      are based on the correspondingly named lattice plots; properties can be controlled 
      using arguments associated with the respective methods. See the details below for
      more information on multidimensional plots.}

  \item{separate}{logical value.  If \code{TRUE}, plot the item category curves or 
      surfaces for polytomous items in separate panels}

  \item{combine}{vector identifying the number of response categories to plot in each
  	panel.  If \code{NULL}, the curves will be grouped by item. \code{combine} is 
	typically used to plot curves for more than one item in a panel. }

  \item{items}{numeric vector identifying the items to plot}

  \item{item.names}{vector of item names for use in labeling the panels }

  \item{item.nums}{logical value.  If \code{TRUE}, include item numbers on the vector plots. }

  \item{panels}{number of panels to display in the output window. If the number of items 
  	is greater than \code{panels}, the plots will be created on multiple pages.  If 
	\code{NULL}, all panels will be plotted on a single page.}

  \item{drift}{character vector identifying the plots to create to examine item
      parameter drift. Acceptable values are "a", "b", "c" for the various parameters
      respectively, "pars" to compare all of these parameters, "TCC" to compare test
      characteristic curves, "ICC" to compare item characteristic curves, or "all" to
      produce all of these plots.}

  \item{groups}{numeric vector identifying the groups in \code{x} that should be included
      for comparing parameter drift. The values should correspond to the group number of 
      the lowest group of each pair of adjacent groups in \code{x}.}

  \item{grp.names}{character vector of group names to use when creating the drift plots}

  \item{sep.mod}{logical value. If \code{TRUE} use different markers in the drift plots
      to identify parameters related to different item response models}

  \item{drift.sd}{numeric value identifying the number of standard deviations to use when
      creating the perpendicular confidence region for the drift comparison plots. If 
      missing, \code{drift.sd} will default to 3.}

  \item{save.hist}{logical value. If \code{FALSE}, the plot history should be deleted
      (only applies to Windows); otherwise, any previously saved plots will be retained.}
}

\details{
All of the plots, with the exception of the vector plots, are based on the \code{lattice} package.
For the unidimensional plots, any arguments associated with the \code{xyplot} function can be
included.  For the multidimensional plots, any arguments associated with the \code{wireframe},
\code{contourplot}, or \code{levelplot} functions (corresponding to the argument specified in 
\code{type}) can be passed as well. 

For the multidimensional plots, vector plots are limited to two dimensions.  The other three
types of plots can be used with multiple dimensions by plotting the results for the first
two dimensions conditional on each combination of theta values for the other dimensions.  The 
function can handle up to ten dimensions; however, the usefulness of the plots will likely
diminish if more than four dimensions are modeled.\cr

There are three types of vector plots.  The first two types create a set of vectors that
originate at the location of the multidimensional difficulty and characterize the length
of the vectors at the appropriate angles, relative to the axes, based on the multidimensional 
discrimination. Setting \code{type} equal to "vectorplot2" will also include the reference
composite based on an eigenvalue decomposition of the slope parameters. The third type of 
vector plot is similar to the first with the key distinction that the lines are drawn from 
the origin to the points of multidimensional difficulty. When there are more than two dimensions,
a panel of vectorplots is created for each pair of dimensions.
}

\note{When multiple pages are created, the PgUp and PgDn buttons can be used to view the 
plots on different pages on a Windows machine. For other operating systems, it is suggested
that the multipage plot be sent to an external graphics file.}

\references{
  Weeks, J. P. (2010) 
  \bold{plink}: An R package for linking mixed-format tests using IRT-based methods. 
  \emph{Journal of Statistical Software}, \bold{35(12)}, 1--33. 
  URL \url{http://www.jstatsoft.org/v35/i12/}
}

\author{ Jonathan P. Weeks \email{weeksjp@gmail.com}, Adam Wyse (Contributing Author)}

\seealso{
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}},
}

\examples{
###### Unidimensional Examples ######
# Compute probabilites for three dichotomous (3PL) items and two polytomous
# (gpcm) items then plot the item characteristic/category curves
dichot <- matrix(c(1.2, .8, .9, 2.3, -1.1, -.2, .24, .19, .13),3,3)
poly <- matrix(c(.64, -1.8, -.73, .45, NA, .88, .06, 1.4, 1.9, 2.6),
  2,5,byrow=TRUE)
pars <- rbind(cbind(dichot,matrix(NA,3,2)),poly)
cat <- c(2,2,2,4,5)
pm <- as.poly.mod(5, c("drm","gpcm"), list(1:3,4:5))
x <- mixed(pars, cat, pm) 
plot(x)

# Compute probabilites for six items using the 3PL and plot all of the item
# characteristic curves on a single panel using the combine argument
a <- c(.71,.96,.36,1.05,1.76,.64)
b <- c(-.16,1.18,2.03,1.09,.82,-1.56)
c <- c(.22,.17,.24,.08,.20,.13)
theta <- seq(-3,3,.2)
pars <- cbind(a,b,c)
x <- drm(pars,theta)
plot(x,combine=6,item.names="Items 1-6",auto.key=list(space="right"))

# Compute probabilties for the nominal response model using the ACT
# mathematics data. Plot the item category curves for a subset of 9 items.
x <- nrm(act.nrm[[1]], rep(5,60))
plot(x, items=c(2,7,12,20,35,41,48,57,60),auto.key=list(space="right"))

# Create irt.pars object with two groups (all dichotomous items)
# rescale the item parameters using the Stocking-Lord linking constants
# Create drift plots 
pm <- as.poly.mod(36)
x <- as.irt.pars(KB04$pars, KB04$common, cat=list(rep(2,36),rep(2,36)), 
  poly.mod=list(pm,pm))
out <- plink(x, rescale="SL", base.grp=2, D=1.7)
plot(out$pars, drift="all", grp.names=c("Form X","Form Y"), 
item.nums=TRUE)

# Compute linking constants for six groups using mixed-format items 
# Create drift plots with distinct markers for the two response models
pm1 <- as.poly.mod(41, c("drm", "gpcm"), reading$items[[1]])
pm2 <- as.poly.mod(70, c("drm", "gpcm"), reading$items[[2]])
pm3 <- as.poly.mod(70, c("drm", "gpcm"), reading$items[[3]])
pm4 <- as.poly.mod(70, c("drm", "gpcm"), reading$items[[4]])
pm5 <- as.poly.mod(72, c("drm", "gpcm"), reading$items[[5]])
pm6 <- as.poly.mod(71, c("drm", "gpcm"), reading$items[[6]])
pm <- list(pm1, pm2, pm3, pm4, pm5, pm6)
x <- as.irt.pars(reading$pars, reading$common, reading$cat, pm, base.grp=4)
out <- plink(x, rescale="HB")
plot(out$pars, drift=c("a","b"), sep.mod=TRUE)

###### Multidimensional Examples ######
# Compute response probabilities for 20 items from the 
# reckase9 data. Probabilites are modeled using the M2PL. 
# Each of the six plot types will be created.
x <- drm(reckase9[[1]][[1]][31:50,],dimensions=2)

# Wireframe plot
plot(x, drape=TRUE, item.names=paste("Item",31:50))
# Contour Plot
plot(x, type="contourplot",labels=FALSE,cuts=20)
# Level Plot
plot(x, type="levelplot",cuts=20)

# Use all the items for the vector plots
x <- drm(reckase9[[1]][[1]],dimensions=2)
# Vector Plot 1
plot(x, type="vectorplot1",item.nums=FALSE)
# Vector Plot 2
plot(x, type="vectorplot2",item.nums=FALSE)
# Vector Plot 3
plot(x, type="vectorplot3",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))


# Compute response probabilities for a single three-category item using 
# the multidimensional generalized partial credit model for three 
# dimensions. Plot the conditional item category surfaces in a single 
# panel and then the separated item category surfaces in separate panels
pars <- matrix(c(1.1999,0.5997,0.8087,2.1730,-1.4576),1,5)
x <- gpcm(pars,3,dimensions=3,theta=-4:4)

# plot combined item category surfaces. Rotate the plots using the
# screen argument
plot(x, screen=list(z=-30,x=-60))

# plot separated item category surfaces
x <- gpcm(pars,3,dimensions=3)
plot(x,separate=TRUE,drape=TRUE,panels=1)
}

\keyword{ misc }
