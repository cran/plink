setGeneric("as.irt.pars", function(x, common, cat, poly.mod, dimensions=1, location=FALSE, grp.names, ...) standardGeneric("as.irt.pars"))

setMethod("as.irt.pars", signature(x="numeric", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	x <- sep.pars(x, cat, poly.mod, ...)
	callGeneric()
})

setMethod("as.irt.pars", signature(x="matrix", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, ...)
	}
	callGeneric()
})

setMethod("as.irt.pars", signature(x="data.frame", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, ...)
	}
	callGeneric()
})

setMethod("as.irt.pars", signature(x="list", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, ...)
	}
	callGeneric()
})

setMethod("as.irt.pars", signature(x="sep.pars", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	ni <- x@n[1]
	dimensions <- x@dimensions
	# Combine the separated item parameters into a matrix
	if (length(x@model)==1) {
		if (x@model=="drm"|x@model=="mcm") {
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
				if (mod=="mcm") pars[items,(n.a+n.b+1):(n.a+n.b+n.c)] <- x@c[items,]
			} else {
				pars[items,(dimensions+1):(dimensions+n.b)] <- x@b[items,]
				if (mod=="drm") pars[items,dimensions+2] <- x@c[items,1]
			}
		}
		tmp <- pars[,ncol(pars)]
		if (length(tmp[is.na(tmp)])==ni) pars <- pars[,-ncol(pars)]
	}
	pm <- as.poly.mod(ni,x@model,x@items)
	cat <- x@cat
	location <- x@loc.out
	
	out <- new("irt.pars",pars=pars,cat=cat,poly.mod=pm,common=NULL,location=location,groups=1,dimensions=dimensions)
	return(out)
})

setMethod("as.irt.pars", signature(x="list", common="matrix"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	tmp <- vector("list",length(x))
	if (length(dimensions)!=length(tmp)) dimensions <- rep(dimensions,length(tmp))
	if (length(location)!=length(tmp)) location <- rep(location,length(tmp))
	if (missing(grp.names)) grp.names <- paste("group",1:length(x),sep="")
	for (i in 1:length(tmp)) {
		if (is.sep.pars(x[[i]]) | is.irt.pars(x[[i]])) {
			tmp[[i]] <- x[[i]]
			next
		} 
		if (location[i]==TRUE) {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], dimensions[i], location[i], loc.out=TRUE, ...)
		} else {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], dimensions[i], location[i])
		}
	}
	return(combine.pars(tmp,common,grp.names))
})

setMethod("as.irt.pars", signature(x="list", common="list"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	tmp <- vector("list",length(x))
	if (length(dimensions)!=length(tmp)) dimensions <- rep(dimensions,length(tmp))
	if (length(location)!=length(tmp)) location <- rep(location,length(tmp))
	if (missing(grp.names)) grp.names <- paste("group",1:length(x),sep="")
	for (i in 1:length(x)) {
		if (is.sep.pars(x[[i]]) | is.irt.pars(x[[i]])) {
			tmp[[i]] <- x[[i]]
			next
		} 
		if (location[i]==TRUE) {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], dimensions[i], location[i], loc.out=TRUE, ...)
		} else {
			tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], dimensions[i], location[i], ...)
		}
	}
	return(combine.pars(tmp,common,grp.names))
})