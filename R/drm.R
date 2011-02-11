##   This function computes response probabilities for items
##   modeled using the Rasch model, 1PL, 2PL, 3PL, M1PL, M2PL, or the M3PL

setGeneric("drm", function(x, theta, dimensions=1, D=1, incorrect=FALSE, print.mod=FALSE, items, ...) standardGeneric("drm"))


##   This method applies when {x} is a  vector of difficulty parameters

setMethod("drm", signature(x="numeric"), function(x, theta, dimensions, D, incorrect, print.mod, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(length(x))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
	
})



setMethod("drm", signature(x="matrix"), function(x, theta, dimensions, D, incorrect, print.mod, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
	
})



setMethod("drm", signature(x="data.frame"), function(x, theta, dimensions, D, incorrect, print.mod, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
	
})



setMethod("drm", signature(x="list"), function(x, theta, dimensions, D, incorrect, print.mod, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
	
})



setMethod("drm", signature(x="irt.pars"), function(x, theta, dimensions, D, incorrect, print.mod, items, ...) {
	
	##   Loop through all groups. In this scenario, a list of {irt.prob} objects will be returned
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], dimensions=x@dimensions[i], ...)
			out[[i]] <- drm(tmp, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, dimensions=x@dimensions, ...)
		callGeneric()
	}
	
})



setMethod("drm", signature(x="sep.pars"), function(x, theta, dimensions, D, incorrect, print.mod, items, ...) {
	
	##   Number of dimensions
	dimensions <- x@dimensions
	
	##   Identify the dichotomous items
	if (missing(items)) items <- 1:x@n[1]
	tmp.items <- x@items$drm
	items <- tmp.items[tmp.items%in%items]
	
	##   Number of items
	n <- length(items)
	
	##   Extract the dichotomous items
	a <- as.matrix(x@a[items,1:dimensions])
	if (n==1) a <- t(a)
	b <- x@b[items,1]
	c <- x@c[items,1]
	pars <- list(a=a, b=b, c=c)
	
	##   Check to see if the argument {D.drm} was passed via the function {mixed}
	dots <- list(...)
	if (length(dots$D.drm)) D <- dots$D.drm
	
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
	
	if (length(x@model[x@model!="drm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the dichotomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	##   Initialize object to hold the response probabilities
	p <- NULL
	
	##   Compute the response probabilities
	for (i in 1:length(b)) {
		if (dimensions==1) {
			##   This is the equation for the 3PL
			cp <- c[i]+(1-c[i])/(1+exp(-D*a[i]*(theta-b[i])))
		} else {
			##   In the multidimensional case D is typically set equal to 1
			a[i,] <- a[i,]*D
			
			##   This is the equation for the M3PL
			cp <- c[i]+(1-c[i])/(1+exp(-(theta %*% a[i,]+b[i])))
		}
		if (incorrect==TRUE) {
			p <- cbind(p,(1-cp),cp)
			colnames(p)[(ncol(p)-1):ncol(p)] <- paste("item_",items[i],".",c(0,1),sep="")
		} else {
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",items[i],".1",sep="")
		}
	}
	
	##   Identify the number of columns in p for each item
	if (incorrect==TRUE) cat <- rep(2,n) else cat <- rep(1,n)
	
	p <- data.frame(cbind(theta,p))
	if (print.mod==TRUE) cat(paste(x@mod.lab,"\n"))
	
	##   Create and return the irt.prob object
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="drm"], dimensions=dimensions, D=c(D.drm=D), pars=pars, model="drm", items=list(drm=1:n))
	return(p)
	
})