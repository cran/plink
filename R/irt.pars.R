##   The as.irt.pars function is used to create an {irt.pars} object for 
##   one or more groups. It includes seven methods.

setGeneric("as.irt.pars", function(x, common, cat, poly.mod, dimensions=1, location=FALSE, grp.names, ...) standardGeneric("as.irt.pars"))


##   This method is used for dichotomous items only where {x} is a 
##   vector of difficulty parameters (this applies for both unidimensional
##   and multidimensional models.

setMethod("as.irt.pars", signature(x="numeric", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {

	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, ...)
	}
	callGeneric()
	
})




##   This method applies when there is only a single group

setMethod("as.irt.pars", signature(x="matrix", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, ...)
	}
	callGeneric()
	
})



##   This method applies when there is only a single group

setMethod("as.irt.pars", signature(x="data.frame", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, ...)
	}
	callGeneric()
	
})



##   This method applies when there is only a single group. It
##   assumes that the item parameters are formatted as a list

setMethod("as.irt.pars", signature(x="list", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	
	if (location==TRUE) {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=TRUE, ...)
	} else {
		x <- sep.pars(x, cat, poly.mod, dimensions, location, ...)
	}
	callGeneric()
	
})



##   This method applies when there is only a single group

setMethod("as.irt.pars", signature(x="sep.pars", common="missing"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	
	##   Identify the total number of items
	ni <- x@n[1]
	
	##   Identify the number of dimensions
	dimensions <- x@dimensions
	
	##   In the {sep.pars} object the item parameters are separated into three elements
	##   Combine the item parameters into a single matrix
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
		
		##   Add colums of NAs (if necessary) so that there will be an adequate 
		##   number of columns for all of the item parameters 
		pars <- matrix(NA,ni,n.a+n.b+n.c)
		
		##   Insert the a-parameters into the first n.a columns
		pars[,1:n.a] <- x@a
		
		##   Loop through the models and insert the b-parameters and c-parameters
		##   into the appropriate rows/columns
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
		
		##   In situations where there are no c-parameters (e.g., when all items are modeled
		##   using the GRM, GPCM, or NRM), the above code will still insert a column
		##   of NAs. Eliminate this column.
		tmp <- pars[,ncol(pars)]
		if (length(tmp[is.na(tmp)])==ni) pars <- pars[,-ncol(pars)]
	}
	pm <- as.poly.mod(ni,x@model,x@items)
	
	##   Create the irt.pars object
	out <- new("irt.pars",pars=pars,cat=x@cat,poly.mod=pm,common=NULL,location=x@location,groups=1,dimensions=dimensions)
	return(out)
	
})



##   This method applies when there are two groups, although
##   the list elements in {x} can be any combination of vectors
##   matrices, data.frames, sep.pars objects, or irt.pars objects

setMethod("as.irt.pars", signature(x="list", common="matrix"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	
	##   Given that common is a matrix, there are necessarily only two groups
	##   Create a temporary list to hold a set of sep.pars objects (or irt.pars objects
	##   if already specified as such) to be combined using the {combine.pars} function
	n <- 2
	tmp <- vector("list",n)
	dots <- list(...)
	
	##   If the number of dimensions is the same across both tests, a single
	##   value for {dimensions} can be specified. The same is true for {location}
	##   repeat these values n times
	if (length(dimensions)!=n) dimensions <- rep(dimensions,n)
	if (length(location)!=n) location <- rep(location,n)
	
	if (missing(grp.names)) grp.names <- paste("group",1:n,sep="")
	
	for (i in 1:n) {
		if (is.sep.pars(x[[i]]) | is.irt.pars(x[[i]])) {
			tmp[[i]] <- x[[i]]
			next
		} 
		
		##  Determine if a location parameter should be output
		if (length(dots$loc.out)) loc.out <- dots$loc.out  else loc.out <- location[i]
		
		##   If {x[[i]]} is not an irt.pars or sep.pars object, make it a sep.pars object
		tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], dimensions[i], location[i], loc.out=loc.out)
	}
	
	##   Return the combined irt.pars object
	return(combine.pars(tmp,common,grp.names))

})


##   This method applies when there are more than two groups.
##   The list elements in {x} can be any combination of vectors
##   matrices, data.frames, sep.pars objects, or irt.pars objects

setMethod("as.irt.pars", signature(x="list", common="list"), function(x, common, cat, poly.mod, dimensions, location, grp.names, ...) {
	
	##   Identify the number of list elements
	n <- length(x)
	
	tmp <- vector("list",n)
	dots <- list(...)
	
	##   If the number of dimensions is the same across all tests, a single
	##   value for {dimensions} can be specified. The same is true for {location}
	##   repeat these values n times
	if (length(dimensions)!=n) dimensions <- rep(dimensions,n)
	if (length(location)!=n) location <- rep(location,n)
	
	if (missing(grp.names)) grp.names <- paste("group",1:n,sep="")
	
	for (i in 1:length(x)) {
		if (is.sep.pars(x[[i]]) | is.irt.pars(x[[i]])) {
			tmp[[i]] <- x[[i]]
			next
		} 
		
		##  Determine if a location parameter should be output
		if (length(dots$loc.out)) loc.out <- dots$loc.out  else loc.out <- location[i]
		
		##   If {x[[i]]} is not an irt.pars or sep.pars object, make it a sep.pars object
		tmp[[i]] <- sep.pars(x[[i]], cat[[i]], poly.mod[[i]], dimensions[i], location[i], loc.out=loc.out)
	}
	
	##   Return the combined irt.pars object
	return(combine.pars(tmp,common,grp.names))
})