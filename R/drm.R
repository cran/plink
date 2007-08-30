setGeneric("drm", function(x, theta=seq(-4,4,.05), D=1.7, incorrect=FALSE, print.mod=FALSE, ...) standardGeneric("drm"))

setMethod("drm", signature(x="numeric"), function(x, theta, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) pm <- as.poly.mod(length(x))
	x <- sep.pars(x, poly.mod=pm, ...)
	callGeneric()
})

setMethod("drm", signature(x="matrix"), function(x, theta, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) pm <- as.poly.mod(nrow(x))
	x <- sep.pars(x, poly.mod=pm, ...)
	callGeneric()
})

setMethod("drm", signature(x="data.frame"), function(x, theta, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) pm <- as.poly.mod(nrow(x))
	x <- sep.pars(x, poly.mod=pm, ...)
	callGeneric()
})

setMethod("drm", signature(x="list"), function(x, theta, D, incorrect, print.mod, ...) {
	if(!hasArg(poly.mod)) pm <- as.poly.mod(nrow(as.matrix(x[[1]])))
	x <- sep.pars(x, poly.mod=pm, ...)
	callGeneric()
})

setMethod("drm", signature(x="irt.pars"), function(x, theta, D, incorrect, print.mod, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], ...)
			out[[i]] <- drm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, ...)
		callGeneric()
	}
})

setMethod("drm", signature(x="sep.pars"), function(x, theta, D, incorrect, print.mod, ...) {
	items <- x@items$drm
	n <- nrow(as.matrix(x@a[items,]))
	a <- x@a[items,1]
	b <- x@b[items,1]
	c <- x@c[items,1]
	
	if (length(x@model[x@model!="drm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the dichotomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	# Compute item probabilities
	p <- NULL 
	for (i in 1:length(b)) {
		cp <- c[i]+(1-c[i])/(1+exp(-D*a[i]*(theta-b[i])))
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
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="drm"], model="drm", items=list(drm=1:n))
	return(p)
})
