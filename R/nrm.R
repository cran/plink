setGeneric("nrm", function(x, cat, theta=seq(-4,4,.05), ...) standardGeneric("nrm"))

setMethod("nrm", signature(x="matrix", cat="numeric"), function(x, cat, theta, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"nrm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, ...)
	callGeneric()
})

setMethod("nrm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"nrm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, ...)
	callGeneric()
})

setMethod("nrm", signature(x="list", cat="numeric"), function(x, cat, theta, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(as.matrix(x[[1]])),"nrm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, ...)
	callGeneric()
})

setMethod("nrm", signature(x="irt.pars"), function(x, cat, theta, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], ...)
			out[[i]] <- nrm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, ...)
		callGeneric()
	}
})

setMethod("nrm", signature(x="sep.pars"), function(x, cat, theta, ...) {
	items <- x@items$nrm
	n <- nrow(as.matrix(x@a[items,]))
	a <- x@a[items,]
	if (is.vector(a)) a <- t(a)
	b <- x@b[items,]
	if (is.vector(b)) b <- t(b)
	
	if (length(x@model[x@model!="nrm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the nrm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	cat <- x@cat[items]
	
	# Compute category probabilities
	p <- NULL 
	for (i in 1:n) {
		den <- NULL # Compute the denominator
		for (k in 1:cat[i]) {
			d <- exp(a[i,k]*theta+b[i,k])
			den <- cbind(den, d)
		}
		den <- apply(den,1,sum)
		for (k in 1:cat[i]) {
			cp <- (exp(a[i,k]*theta+b[i,k]))/den
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",i,".",k,sep="")
		}
	}
	p <- data.frame(cbind(theta,p))
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab="Nominal Response Model", model="nrm", items=list(nrm=1:n))
	return(p)
})
