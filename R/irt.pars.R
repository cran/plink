setGeneric("as.irt.pars", function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL, ...) standardGeneric("as.irt.pars"))

setMethod("as.irt.pars", signature(x="numeric", common="missing"), function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL,...) {
	x <- sep.pars(x, ...)
	callGeneric()
})

setMethod("as.irt.pars", signature(x="matrix", common="missing"), function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL,...) {
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, location=TRUE, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, ...)
	}
	callGeneric()
})

setMethod("as.irt.pars", signature(x="data.frame", common="missing"), function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL,...) {
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, location=TRUE, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, ...)
	}
	callGeneric()
})

setMethod("as.irt.pars", signature(x="list", common="missing"), function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL,...) {
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, location=TRUE, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, ...) 
	}
	callGeneric()
})

setMethod("as.irt.pars", signature(x="sep.pars", common="missing"), function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL,...) {
	ni <- x@n[1]
	if (length(x@model)==1) {
		if (x@model=="drm") {
			pars <- cbind(x@a,x@b,x@c)
		} else {
			pars <- cbind(x@a,x@b)
		}
	} else {
		n.a <- ncol(x@a)
		n.b <- ncol(x@b)
		n.c <- ncol(x@c)
		pars <- matrix(NA,ni,n.a+n.b+n.c)
		pars[,1:n.a] <- x@a
		for (j in 1:length(x@model)) {
			mod <- x@model[j]
			items <- x@items[[j]]
			if (mod=="nrm"|mod=="mcm") {
				pars[items,(n.a+1):(n.a+n.b)] <- x@b[items,]
				pars[items,(n.a+n.b+1):(n.a+n.b+n.c)] <- x@c[items,]
			} else {
				pars[items,2:(n.b+1)] <- x@b[items,]
				if (mod=="drm") pars[items,3] <- x@c[items,1]
			}
		}
		tmp <- pars[,ncol(pars)]
		if (length(tmp[is.na(tmp)])==ni) pars <- pars[,-ncol(pars)]
	}
	pm <- as.poly.mod(ni,x@model,x@items)
	cat <- x@cat
	location <- x@loc.out
	
	out <- new("irt.pars",pars=pars,cat=cat,poly.mod=pm,common=NULL,location=location,groups=1)
	return(out)
})

setMethod("as.irt.pars", signature(x="list", common="matrix"), function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL,...) {
	tmp <- vector("list",length(x))
	for (i in 1:length(x)) {
		if (is.sep.pars(x[[i]]) | is.irt.pars(x[[i]])) {
			tmp[[i]] <- x[[i]]
			next
		} 
		if (location==TRUE) {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], location=TRUE, loc.out=TRUE, ...)
		} else {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], ...)
		}
	}
	return(combine.pars(tmp,common,grp.names))
})

setMethod("as.irt.pars", signature(x="list", common="list"), function(x, common, cat, poly.mod, location=FALSE, grp.names=NULL,...) {
	tmp <- vector("list",length(x))
	for (i in 1:length(x)) {
		if (is.sep.pars(x[[i]]) | is.irt.pars(x[[i]])) {
			tmp[[i]] <- x[[i]]
			next
		} 
		if (location==TRUE) {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], location=TRUE, loc.out=TRUE, ...)
		} else {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], ...)
		}
	}
	return(combine.pars(tmp,common,grp.names))
})

