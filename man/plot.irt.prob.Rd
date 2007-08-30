\name{plot.irt.prob}
\alias{plot.irt.prob}
\alias{plot.sep.pars}
\alias{plot.irt.pars}
\title{Plot Item Characteristic/Category Curves}
\description{
This function plots item characteristic curves or item category curves using the
\code{lattice} package.
}

\usage{
\method{plot}{irt.prob}(x, y, ..., combine = NULL, item.names = NULL, item.lab = TRUE, 
	panels = 20)

\method{plot}{sep.pars}(x, y, ...)

\method{plot}{irt.pars}(x, y, ...)
}

\arguments{
  \item{x}{object of class \code{irt.prob}, \code{irt.pars}, or \code{sep.pars}.  For the
  later two classes, probabilities are computed using the function \code{mixed} before
  plotting the curves.}

  \item{y}{this argument is not used for these methods}

  \item{\dots}{further arguments passed to or from other methods}

  \item{combine}{vector identifying the number of response categories to plot in each
  	panel.  If \code{NULL}, the curves will be grouped by item. \code{combine} is 
	typically used to plot curves for more than one item in a panel. }

  \item{item.names}{vector of item names for use in labeling the panels. }

  \item{item.lab}{if \code{TRUE}, include the item label strip at the top of each panel}

  \item{panels}{number of panels to display in the output window. If the number of items 
  	is greater than \code{panels}, the plots will be created on multiple pages.  If 
	\code{NULL}, all panels will be plotted on a single page.}
}

\author{ Jonathan Weeks \email{weeksjp@gmail.com} }

\note{When multiple pages are created, the PgUp and PgDn buttons can be used to view the 
plots on different pages}

\seealso{
\code{\linkS4class{irt.prob}}, \code{\linkS4class{irt.pars}}, \code{\linkS4class{sep.pars}},
}

\keyword{ misc }
