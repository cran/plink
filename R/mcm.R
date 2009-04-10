setGeneric("mcm", function(x, cat, theta, dimensions=1, ...) standardGeneric("mcm"))

setMethod("mcm", signature(x="matrix", cat="numeric"), function(x, cat, theta, dimensions, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"mcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
})

setMethod("mcm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, dimensions, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"mcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
})

setMethod("mcm", signature(x="list", cat="numeric"), function(x, cat, theta, dimensions, ...) {
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])),"mcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
})

setMethod("mcm", signature(x="irt.pars"), function(x, cat, theta, dimensions, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], dimensions=x@dimensions[i], ...)
			out[[i]] <- mcm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, dimensions=x@dimensions, ...)
		callGeneric()
	}
})

setMethod("mcm", signature(x="sep.pars"), function(x, cat, theta, dimensions, ...) {
	items <- x@items$mcm
	n <- length(items)
	dimensions <- x@dimensions
	a <- x@a[items,] #group by dimensions then categories. i.e. for dimension m and category k (amk)
	b <- x@b[items,]
	c <- x@c[items,]
	pars <- list(a=a, b=b, c=c)
	if (n==1) {
		a <- t(a)
		b <- t(b)
		c <- t(c)
	}
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
		if (is.numeric(theta)) {
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
	
	if (length(x@model[x@model!="mcm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the mcm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	# Compute category probabilities
	p <- NULL 
	for (i in 1:n) {
		den <- NULL # Compute the denominator
		a1 <- a[i,][!is.na(a[i,])]
		b1 <- b[i,][!is.na(b[i,])]
		c1 <- c[i,][!is.na(c[i,])]
		for (k in 1:cat[i]) {
			tmp <- (k-1)*dimensions
			tmp1 <- tmp+dimensions
			d <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])
			den <- cbind(den, d)
		}
		den <- apply(den,1,sum)
		for (k in 1:cat[i]) {
			tmp <- (k-1)*dimensions
			tmp1 <- tmp+dimensions
			if (k==1) {
				cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k]))/den
			} else {
				cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k])+c1[k-1]*(exp((theta %*% a1[1:dimensions])+b1[1])))/den
			}
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",i,".",k-1,sep="")
		}
	}
	p <- data.frame(cbind(theta,p))
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="mcm"], dimensions=dimensions, pars=pars, model="mcm", items=list(mcm=1:n))
	return(p)
})