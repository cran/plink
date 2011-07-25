##   This function computes response probabilities for items
##   modeled using the partial credit model, generalized
##   partial credit model, multidimensional partial credit model,
##   and the multidimensional generalized partial credit model

setGeneric("gpcm", function(x, cat, theta, dimensions=1, D=1, location=FALSE, print.mod=FALSE, items, ...) standardGeneric("gpcm"))



setMethod("gpcm", signature(x="matrix", cat="numeric"), function(x, cat, theta, dimensions, D, location, print.mod, items, ...) {

	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"gpcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("gpcm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, dimensions, D, location, print.mod, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"gpcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("gpcm", signature(x="list", cat="numeric"), function(x, cat, theta, dimensions, D, location, print.mod, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])),"gpcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



##   For this method the objects, cat, dimensions, and location are contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("gpcm", signature(x="irt.pars", cat="ANY"), function(x, cat, theta, dimensions, D, location, print.mod, items, ...) {
	
	##   Loop through all groups. In this scenario, a list of {irt.prob} objects will be returned
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], x@dimensions[i], x@location[i], loc.out=FALSE, ...)
			out[[i]] <- gpcm(tmp, ...)
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
setMethod("gpcm", signature(x="sep.pars", cat="ANY"), function(x, cat, theta, dimensions, D, location, print.mod, items, ...) {
	
	##   The equation to compute probabilities is not (actually) parameterized using
	##   the location/step-deviation formulation. As such, in instances where a location
	##   parameter is included, it is necessary to reformat the parameters appropriately
	if (x@location==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(x@n[1], x@model, x@items)
		x <- sep.pars(pars, x@cat, pm, x@dimensions, location=TRUE, loc.out=FALSE)
	}
	
	##   Number of dimensions
	dimensions <- x@dimensions
	
	##   Identify the gpcm items
	if (missing(items)) items <- 1:x@n[1]
	tmp.items <- x@items$gpcm
	items <- tmp.items[tmp.items%in%items]
	
	##   Number of items
	n <- length(items)
	
	##   Extract the gpcm items 
	a <- as.matrix(x@a[items,1:dimensions])
	b <- as.matrix(x@b[items,])
	
	##   If there is only a single item, the matrices specified above will have
	##   the wrong orientation. For example, the step parameters for this item
	##   will be in different rows of the matrix b instead of being in a matrix
	##   with a single row and multiple columns. As such, these matrices need
	##   to be transposed
	if (n==1) {
		a <- t(a)
		b <- t(b)
	}
	pars <- list(a=a, b=b, c=x@c[items,])
	
	cat <- x@cat[items]
	
	##   Check to see if the argument {D.gpcm} was passed via the function {mixed}
	dots <- list(...)
	if (length(dots$D.gpcm)) D <- dots$D.gpcm
	
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
	
	if (length(x@model[x@model!="gpcm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the gpcm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	##   Initialize object to hold the response probabilities
	p <- NULL 
	
	for (i in 1:nrow(b)) {
		##   Object used to accumulate step parameters
		dif <- 0 
		
		##   Object for the denominator in the final equation
		den <- NULL 
		
		##   Compute the denominator
		for (k in 0:(cat[i]-1)) {
			if (k>=1) dif <- dif+b[i,k]
			if (dimensions==1) {
				d <- exp(D*a[i]*(k*theta-dif))
			} else {
				d <- exp(k*(theta %*% a[i,])+dif)
			}
			den <- cbind(den, d)
		}
		den <- apply(den,1,sum)
		
		##   Compute the response probabilities
		dif <- 0
		for (k in 0:(cat[i]-1)) {
			if (k>=1) dif <- dif+b[i,k]
			if (dimensions==1) {
				##   This is the equation for the generalized partial credit model
				cp <- (exp(D*a[i]*(k*theta-dif)))/den
			} else {
				##   This is the equation for the MD generalized partial credit model
				cp <- exp(k*(theta %*% a[i,])+dif)/den
			}
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",items[i],".",k,sep="")
		}
	}
		
	p <- data.frame(cbind(theta,p))
	if (print.mod==TRUE) cat(paste(x@mod.lab,"\n"))
	
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="gpcm"], dimensions=dimensions, D=c(D.gpcm=D), pars=pars, model="gpcm", items=list(gpcm=1:n))
	return(p)
	
})