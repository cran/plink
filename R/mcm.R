setGeneric("mcm", function(x, cat, theta=seq(-4,4,.05), ...) standardGeneric("mcm"))

setMethod("mcm", signature(x="matrix", cat="numeric"), function(x, cat, theta, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"mcm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, ...)
	callGeneric()
})

setMethod("mcm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"mcm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, ...)
	callGeneric()
})

setMethod("mcm", signature(x="list", cat="numeric"), function(x, cat, theta, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(as.matrix(x[[1]])),"mcm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, ...)
	callGeneric()
})

setMethod("mcm", signature(x="irt.pars"), function(x, cat, theta, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], ...)
			out[[i]] <- mcm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, ...)
		callGeneric()
	}
})

setMethod("mcm", signature(x="sep.pars"), function(x, cat, theta, ...) {
	items <- x@items$mcm
	n <- nrow(as.matrix(x@a[items,]))
	a <- x@a[items,]
	if (is.vector(a)) a <- t(a)
	b <- x@b[items,]
	if (is.vector(b)) b <- t(b)
	c <- x@c[items,]
	if (is.vector(c)) c <- t(c)
	
	if (length(x@model[x@model!="mcm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the mcm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
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
			if (k==1) {
				cp <- (exp(a[i,k]*theta+b[i,k]))/den
			} else {
				cp <- (exp(a[i,k]*theta+b[i,k])+c[i,(k-1)]*(exp(a[i,1]*theta+b[i,1])))/den
			}
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",i,".",k-1,sep="")
		}
	}
	p <- data.frame(cbind(theta,p))
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab="Multiple-Choice Model", model="mcm", items=list(mcm=1:n))
	return(p)
})
