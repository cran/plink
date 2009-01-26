setGeneric("grm", function(x, cat, theta=seq(-4,4,.05), catprob=FALSE, D=1.7, location=FALSE, ...) standardGeneric("grm"))

setMethod("grm", signature(x="matrix", cat="numeric"), function(x, cat, theta, catprob, D, location, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"grm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, location, ...)
	callGeneric()
})

setMethod("grm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, catprob, D, location, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"grm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, location, ...)
	callGeneric()
})

setMethod("grm", signature(x="list", cat="numeric"), function(x, cat, theta, catprob, D, location, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(as.matrix(x[[1]])),"grm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, location, ...)
	callGeneric()
})

setMethod("grm", signature(x="irt.pars"), function(x, cat, theta, catprob, D, location, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], ...)
			out[[i]] <- grm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, ...)
		callGeneric()
	}
})

setMethod("grm", signature(x="sep.pars"), function(x, cat, theta, catprob, D, location, ...) {
	if (x@loc.out==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(list(x@model,x@items))
		x <- sep.pars(pars, x@cat, pm, location=TRUE)
	}
	items <- x@items$grm
	n <- nrow(as.matrix(x@a[items,]))
	a <- x@a[items,1]
	b <- x@b[items,]
	if (is.vector(b)) b <- t(b)
	
	if (length(x@model[x@model!="grm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the grm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	cat <- x@cat[items]
	
	p <- NULL 
	# Compute category probabilities
	if (catprob==TRUE) { 
		for (i in 1:n) {
			ct <- cat[i]-1
			cp <- 1-1/(1+exp(-D*a[i]*(theta-b[i,1])))
			p <- cbind(p, cp)
			colnames(p)[ncol(p)] <- paste("p",i,".0",sep="")
			for (k in 1:ct) {
				if (k<ct) {
					cp <- (1/(1+exp(-D*a[i]*(theta-b[i,k]))))-(1/(1+exp(-D*a[i]*(theta-b[i,k+1]))))
				} else if (k==ct) {
					cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
				}
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",i,".",k-1,sep="")
			}
		}
	# Compute cumulative probabilities
	} else if (catprob==FALSE) { 
		for (i in 1:n) {
			for (k in 1:(cat[i]-1)) {
				cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",i,".",k,sep="")
			}
		}
	}
	p <- data.frame(cbind(theta,p))
	if (catprob==FALSE) cat <- cat-1
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab="Graded Response Model", model="grm", items=list(grm=1:n))
	return(p)
})
