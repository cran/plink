setGeneric("drm", function(x, theta, dimensions=1, D=1, incorrect=FALSE, print.mod=FALSE, ...) standardGeneric("drm"))

setMethod("drm", signature(x="numeric"), function(x, theta, dimensions, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(length(x))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
})

setMethod("drm", signature(x="matrix"), function(x, theta, dimensions, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
})

setMethod("drm", signature(x="data.frame"), function(x, theta, dimensions, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
})

setMethod("drm", signature(x="list"), function(x, theta, dimensions, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])))
	x <- sep.pars(x, poly.mod=poly.mod, dimensions=dimensions, ...)
	callGeneric()
})

setMethod("drm", signature(x="irt.pars"), function(x, theta, dimensions, D, incorrect, print.mod, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], dimensions=x@dimensions[i], ...)
			out[[i]] <- drm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, dimensions=x@dimensions, ...)
		callGeneric()
	}
})

setMethod("drm", signature(x="sep.pars"), function(x, theta, dimensions, D, incorrect, print.mod, ...) {
	dots <- list(...)
	if (length(dots$D.drm)) D <- dots$D.drm
	items <- x@items$drm
	n <- length(items)
	dimensions <- x@dimensions
	a <- as.matrix(x@a[items,1:dimensions])*D
	if (n==1) a <- t(a)
	b <- x@b[items,1]
	c <- x@c[items,1]
	pars <- list(a=a/D, b=b, c=c)
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
		b <- -b*matrix(a,n,1)
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
	
	# Compute item probabilities
	p <- NULL
	for (i in 1:length(b)) {
		cp <- c[i]+(1-c[i])/(1+exp(-(theta %*% a[i,]+b[i])))
		if (incorrect==TRUE) {
			p <- cbind(p,(1-cp),cp)
			cat <- rep(2,n)
			colnames(p)[(ncol(p)-1):ncol(p)] <- paste("item_",i,".",c(0,1),sep="")
		} else {
			p <- cbind(p,cp)
			cat <- rep(1,n)
			colnames(p)[ncol(p)] <- paste("item_",i,".1",sep="")
		}
	}
	p <- data.frame(cbind(theta,p))
	if (print.mod==TRUE) cat(paste(x@mod.lab,"\n"))
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="drm"], dimensions=dimensions, pars=pars, model="drm", items=list(drm=1:n))
	return(p)
})