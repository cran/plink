##   This function computes response probabilities for items
##   modeled using the graded response model and the 
##   multidimensional graded response model

setGeneric("grm", function(x, cat, theta, dimensions=1, catprob=FALSE, D=1, location=FALSE, items, ...) standardGeneric("grm"))



setMethod("grm", signature(x="matrix", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("grm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("grm", signature(x="list", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



##   For this method the objects, cat, dimensions, and location are contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("grm", signature(x="irt.pars", cat="ANY"), function(x, cat, theta, dimensions, catprob, D, location, items, ...) {
	
	##   Loop through all groups. In this scenario, a list of {irt.prob} objects will be returned
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], x@dimensions[i], x@location[i], loc.out=FALSE, ...)
			out[[i]] <- grm(tmp, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, x@dimensions, x@location, loc.out=FALSE, ...)
		callGeneric()
	}
	
})



##   For this method the objects, cat, dimensions, and location are contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("grm", signature(x="sep.pars", cat="ANY"), function(x, cat, theta, dimensions, catprob, D, location, items, ...) {

	##   The equation to compute probabilities is not (actually) parameterized using
	##   the location/threshold-deviation formulation. As such, in instances where a location
	##   parameter is included, it is necessary to reformat the parameters appropriately
	if (x@location==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(x@n[1], x@model, x@items)
		x <- sep.pars(pars, x@cat, pm, x@dimensions, location=TRUE, loc.out=FALSE)
	}
	
	##   Number of dimensions
	dimensions <- x@dimensions
	
	##   Identify the grm items
	if (missing(items)) items <- 1:x@n[1]
	tmp.items <- x@items$grm
	items <- tmp.items[tmp.items%in%items]
	
	##   Number of items
	n <- length(items)
	
	##   Extract the grm items and rescale the slope parameters (if necessary)
	a <- as.matrix(x@a[items,1:dimensions])
	b <- as.matrix(x@b[items,])
	
	##   If there is only a single item, the matrices specified above will have
	##   the wrong orientation. For example, the threshold parameters for this item
	##   will be in different rows of the matrix b instead of being in a matrix
	##   with a single row and multiple columns. As such, these matrices need
	##   to be transposed
	if (n==1) {
		a <- t(a)
		b <- t(b)
	}
	pars <- list(a=a, b=b, c=x@c[items,])
	
	cat <- x@cat[items]
	
	##   Check to see if the argument {D.grm} was passed via the function {mixed}
	dots <- list(...)
	if (length(dots$D.grm)) D <- dots$D.grm
	
	##   Generate theta values if {theta} is missing
	##   Different values should be generated depending on the number of dimensions
	if (missing(theta)) {
		if (dimensions==1) {
			theta <- seq(-4,4,.05) 
		} else if (dimensions %in% 2:3) {
			theta <- seq(-4,4,.5)
		} else {
			theta <- -4:4
		}
	}
	
	if (dimensions==1) {
		##   If the user (purposefully or accidentally) specifies {theta} as a matrix
		##   or a list instead of a vector for the unidimensional case, turn all of the 
		##   values into a vector
		if (is.matrix(theta)) {
			if (ncol(theta)>1) {
				theta <- as.vector(theta)
			}
		} else if (is.list(theta)) {
			theta <- unlist(theta)
		}
		theta <- as.matrix(theta)
		colnames(theta) <- "theta1"
		
	}else if (dimensions>1) {
		##   If, in the multidimensional case, only a vector of theta values is 
		##   supplied, treat this as a vector for each dimension then create all
		##   permutations of these values. If {theta} is formatted as a matrix
		##   or list from the outset, just find the permutations
		if (is.vector(theta)) {
			tmp <- vector("list", dimensions)
			for (i in 1:dimensions) {
				tmp[[i]] <- theta
			}
			theta <- as.matrix(expand.grid(tmp))
			colnames(theta) <- paste("theta",1:dimensions,sep="")
		} else if (is.list(theta)) {
			theta <- as.matrix(expand.grid(theta))
			colnames(theta) <- paste("theta",1:dimensions,sep="")
		} else if (is.matrix(theta)) {
			if (ncol(theta)>1) {
				colnames(theta) <- paste("theta",1:dimensions,sep="")
			} else {
				tmp <- vector("list", dimensions)
				for (i in 1:dimensions) {
					tmp[[i]] <- theta
				}
				theta <- as.matrix(expand.grid(tmp))
				colnames(theta) <- paste("theta",1:dimensions,sep="")
			}
		}
	}
	
	if (length(x@model[x@model!="grm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the grm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	##   Initialize object to hold the response probabilities
	p <- NULL 
	
	##   Compute category probabilities
	if (catprob==TRUE) { 
		for (i in 1:n) {
			ct <- cat[i]-1
			
			##   Compute the probabilities for the lowest category
			if (dimensions==1) {
				cp <- 1-1/(1+exp(-D*a[i]*(theta-b[i,1])))
			} else {
				cp <- 1-1/(1+exp(-(theta %*% a[i,]+b[i,1])))
			}
			p <- cbind(p, cp)
			colnames(p)[ncol(p)] <- paste("item_",items[i],".0",sep="")
			
			for (k in 1:ct) {
				if (dimensions==1) {
					if (k<ct) {
						cp <- (1/(1+exp(-D*a[i]*(theta-b[i,k]))))-(1/(1+exp(-D*a[i]*(theta-b[i,k+1]))))
					} else if (k==ct) {
						##   Compute the probabilities for the highest category
						cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
					}
				} else {
					if (k<ct) {
						cp <- (1/(1+exp(-(theta %*% a[i,]+b[i,k]))))-(1/(1+exp(-(theta %*% a[i,]+b[i,k+1]))))
					} else if (k==ct) {
						##   Compute the probabilities for the highest category
						cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
					}
				}
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",items[i],".",k,sep="")
			}
		}
		
	# Compute cumulative probabilities
	} else if (catprob==FALSE) { 
		for (i in 1:n) {
			for (k in 1:(cat[i]-1)) {
				if (dimensions==1) {
					cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
				} else {
					cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
				}
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",items[i],".",k,sep="")
			}
		}
	}
	
	p <- data.frame(cbind(theta,p))
	if (catprob==FALSE) cat <- cat-1
	
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="grm"], dimensions=dimensions, D=c(D.grm=D), pars=pars, model="grm", items=list(grm=1:n))
	return(p)
	
})