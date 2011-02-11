##   This function computes response probabilities for items
##   modeled using the multiple-choice model and the 
##   multidimensional multiple-choice model

setGeneric("mcm", function(x, cat, theta, dimensions=1, items, ...) standardGeneric("mcm"))



setMethod("mcm", signature(x="matrix", cat="numeric"), function(x, cat, theta, dimensions, items, ...) {

	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"mcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
	
})



setMethod("mcm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, dimensions, items, ...) {

	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"mcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
	
})



setMethod("mcm", signature(x="list", cat="numeric"), function(x, cat, theta, dimensions, items, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])),"mcm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
	
})



##   For this method the objects, cat and dimensionsare contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("mcm", signature(x="irt.pars", cat="ANY"), function(x, cat, theta, dimensions, items, ...) {

	##   Loop through all groups. In this scenario, a list of {irt.prob} objects will be returned
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], dimensions=x@dimensions[i], ...)
			out[[i]] <- mcm(tmp, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, dimensions=x@dimensions, ...)
		callGeneric()
	}
	
})



##   For this method the objects, cat and dimensionsare contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("mcm", signature(x="sep.pars", cat="ANY"), function(x, cat, theta, dimensions, items, ...) {

	##   Identify the mcm items
	if (missing(items)) items <- 1:x@n[1]
	tmp.items <- x@items$mcm
	items <- tmp.items[tmp.items%in%items]
	
	##   Number of items
	n <- length(items)
	
	dimensions <- x@dimensions
	
	##   Extract the mcm items
	##   When creating the object with the item parameters they should be 
	##   grouped first by dimensions then by categories.  For example, for 
	##   an item with 2 dimensions and 4 "actual" categories we would have 
	##   (a11,a12,a13,a14,a15,a21,a22,a23,a24,a25) where a11 and a21
	##   correspond to the "do not know" category for each dimension.
	##   This grouping applies for both the slope parameters (named a)
	##   and the category parameters (named b)
	a <- x@a[items,] 
	b <- x@b[items,]
	c <- x@c[items,]
	
	##   If there is only a single item, the matrices specified above will have
	##   the wrong orientation. For example, the slope parameters for this item
	##   will be in different rows of the matrix a instead of being in a matrix
	##   with a single row and multiple columns. As such, these matrices need
	##   to be transposed
	if (n==1) {
		a <- t(a)
		b <- t(b)
		c <- t(c)
	}
	
	pars <- list(a=a, b=b, c=c)
	cat <- x@cat[items]
	
	##   Generate theta values if {theta} is missing
	##   Different values should be generated depending on the number of dimensions
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
		##   If the user (purposefully or accidentally) specifies {theta} as a matrix
		##   or a list instead of a vector for the unidimensional case, turn all of the 
		##   values into a vector
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
			##   If, in the multidimensional case, only a vector of theta values is 
			##   supplied, treat this as a vector for each dimension then create all
			##   permutations of these values. If {theta} is formatted as a matrix
			##   or list from the outset, just find the permutations
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
	
	##   Initialize object to hold the response probabilities
	p <- NULL 
	
	for (i in 1:n) {
		##   Object for the denominator in the final MMCM equation
		den <- NULL
		
		##   Because of how the parameters are organized in {x}
		##   there may be NAs in various rows. Remove these NAs
		##   before computing the response probabilites
		a1 <- a[i,][!is.na(a[i,])]
		b1 <- b[i,][!is.na(b[i,])]
		c1 <- c[i,][!is.na(c[i,])]
		
		##   Compute the denominator
		for (k in 1:cat[i]) {
			tmp <- (k-1)*dimensions
			tmp1 <- tmp+dimensions
			d <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])
			den <- cbind(den, d)
		}
		den <- apply(den,1,sum)
		
		##   Compute the response probabilities
		for (k in 1:cat[i]) {
			tmp <- (k-1)*dimensions
			tmp1 <- tmp+dimensions
			if (k==1) {
				cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k]))/den
			} else {
				cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k])+c1[k-1]*(exp((theta %*% a1[1:dimensions])+b1[1])))/den
			}
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",items[i],".",k-1,sep="")
		}
	}
	p <- data.frame(cbind(theta,p))
	
	p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="mcm"], dimensions=dimensions, D=c(D=1), pars=pars, model="mcm", items=list(mcm=1:n))
	return(p)
	
})