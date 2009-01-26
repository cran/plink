setGeneric("gpcm", function(x, cat, theta=seq(-4,4,.05), D=1.7, location=FALSE, print.mod=FALSE, ...) standardGeneric("gpcm"))

setMethod("gpcm", signature(x="matrix", cat="numeric"), function(x, cat, theta, D, location, print.mod, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"gpcm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, location, ...)
	callGeneric()
})

setMethod("gpcm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, D, location, print.mod, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(x),"gpcm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, location, ...)
	callGeneric()
})

setMethod("gpcm", signature(x="list", cat="numeric"), function(x, cat, theta, D, location, print.mod, ...) {
	dots <- list(...)
	if(is.null(dots$poly.mod)) pm <- as.poly.mod(nrow(as.matrix(x[[1]])),"gpcm") else pm <- dots$poly.mod
	x <- sep.pars(x, cat, pm, location, ...)
	callGeneric()
})

setMethod("gpcm", signature(x="irt.pars"), function(x, cat, theta, D, location, print.mod, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], ...)
			out[[i]] <- gpcm(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, ...)
		callGeneric()
	}
})

setMethod("gpcm", signature(x="sep.pars"), function(x, cat, theta, D, location, print.mod, ...) {
	if (x@loc.out==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(list(x@model,x@items))
		x <- sep.pars(pars, x@cat, pm, location=TRUE)
	}
	items <- x@items$gpcm
	n <- nrow(as.matrix(x@a[items,]))
	a <- x@a[items,1]
	b <- x@b[items,]
	if (is.vector(b)) b <- t(b)
	
	if (length(x@model[x@model!="gpcm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the gpcm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	cat <- x@cat[items]
	
	# Compute category probabilities
	p <- NULL 
	for (i in 1:nrow(b)) {
		dif <- 0 # Difference between subsequent step parameters
		den <- NULL # Compute the denominator
		for (k in 0:(cat[i]-1)) {
			if (k==1) dif <- b[i,k] else if (k>1) dif <- dif+b[i,k]
			d <- exp(D*a[i]*(k*theta-dif))
			den <- cbind(den, d)
		}
		den <- apply(den,1,sum)
		dif <- 0
		for (k in 0:(cat[i]-1)) {
			if (k==1) dif <- b[i,k] else if (k>1) dif <- dif+b[i,k]
			cp <- (exp(D*a[i]*(k*theta-dif)))/den
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",i,".",k,sep="")
		}
	}
	p <- data.frame(cbind(theta,p))
	if (print.mod==TRUE) cat(paste(x@mod.lab,"\n"))
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="gpcm"], model="gpcm", items=list(gpcm=1:n))
	return(p)
})
