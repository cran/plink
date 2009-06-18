setGeneric("grm", function(x, cat, theta, dimensions=1, catprob=FALSE, D=1, location=FALSE, ...) standardGeneric("grm"))

setMethod("grm", signature(x="matrix", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
})

setMethod("grm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
})

setMethod("grm", signature(x="list", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
})

setMethod("grm", signature(x="irt.pars"), function(x, cat, theta, dimensions, catprob, D, location, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], x@dimensions[i], x@location[i], loc.out=FALSE, ...)
			out[[i]] <- grm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, x@dimensions, x@location, loc.out=FALSE, ...)
		callGeneric()
	}
})

setMethod("grm", signature(x="sep.pars"), function(x, cat, theta, dimensions, catprob, D, location, ...) {
	if (x@loc.out==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(x@n[1], x@model, x@items)
		x <- sep.pars(pars, x@cat, pm, x@dimensions, location=TRUE, loc.out=FALSE)
	}
	dots <- list(...)
	if (length(dots$D.grm)) D <- dots$D.grm
	items <- x@items$grm
	n <- length(items)
	dimensions <- x@dimensions
	a <- as.matrix(x@a[items,1:dimensions])*D
	b <- as.matrix(x@b[items,])
	if (n==1) {
		a <- t(a)
		b <- t(b)
	}
	pars <- list(a=a/D, b=b, c=x@c[items,])
	cat <- x@cat[items]
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
		b <- -b*matrix(a,n,ncol(b))
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
	
	if (length(x@model[x@model!="grm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the grm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	p <- NULL 
	# Compute category probabilities
	if (catprob==TRUE) { 
		for (i in 1:n) {
			ct <- cat[i]-1
			cp <- 1-1/(1+exp(-(theta %*% a[i,]+b[i,1])))
			p <- cbind(p, cp)
			colnames(p)[ncol(p)] <- paste("item_",i,".0",sep="")
			for (k in 1:ct) {
				if (k<ct) {
					cp <- (1/(1+exp(-(theta %*% a[i,]+b[i,k]))))-(1/(1+exp(-(theta %*% a[i,]+b[i,k+1]))))
				} else if (k==ct) {
					cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
				}
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",i,".",k,sep="")
			}
		}
	# Compute cumulative probabilities
	} else if (catprob==FALSE) { 
		for (i in 1:n) {
			for (k in 1:(cat[i]-1)) {
				cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",i,".",k,sep="")
			}
		}
	}
	p <- data.frame(cbind(theta,p))
	if (catprob==FALSE) cat <- cat-1
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="grm"], dimensions=dimensions, pars=pars, model="grm", items=list(grm=1:n))
	return(p)
})