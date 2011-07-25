##   This function computes response probabilities for items
##   modeled using any combination of the models available 
##   in the functions {drm}, {gpcm}, {grm}, {mcm}, and {nrm}

setGeneric("mixed", function(x, cat, poly.mod, theta, dimensions=1, items, ...) standardGeneric("mixed"))



setMethod("mixed", signature(x="numeric", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, items, ...) {
	
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("mixed", signature(x="matrix", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, items, ...) {
	
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("mixed", signature(x="data.frame", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, items, ...) {
	
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("mixed", signature(x="list", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, items, ...) {
	
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
	
})



##   For this method the objects, cat and dimensionsare contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("mixed", signature(x="irt.pars", cat="ANY"), function(x, cat, poly.mod, theta, dimensions, items, ...) {
	
	##   Loop through all groups. In this scenario, a list of {irt.prob} objects will be returned
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], loc.out=FALSE, location=x@location[i], dimensions=x@dimensions[i], ...)
			out[[i]] <- mixed(tmp, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, loc.out=FALSE, location=x@location, dimensions=x@dimensions, items, ...)
		callGeneric()
	}
})



##   For this method the objects, cat and dimensionsare contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("mixed", signature(x="sep.pars", cat="ANY"), function(x, cat, poly.mod, theta, dimensions, items, ...) {

	##   The equation to compute probabilities is not (actually) parameterized using
	##   the location/(step or threshold)-deviation formulation. As such, in instances where a location
	##   parameter is included, it is necessary to reformat the parameters appropriately
	if (x@location==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(x@n[1], x@model, x@items)
		x <- sep.pars(pars, x@cat, pm, location=TRUE, loc.out=FALSE, x@dimensions, ...)
	}
	
	##   Total number of items
	if (missing(items)) {
		items <- seq(1,length(x@cat))
		item.nms <- items
	} else {
		## Recreate the set of item parameters
		tmp <- as.irt.pars(x)
		tmp@pars <- tmp@pars[items,]
		if (length(items)==1) tmp@pars <- t(tmp@pars)
		tmp@pars <- tmp@pars[,apply(!is.na(tmp@pars),2,sum)>0]
		if (length(items)==1) tmp@pars <- t(tmp@pars)
		tmp@cat <- tmp@cat[items]
		
		n <- length(items)
		pm <- tmp@poly.mod
		mod <- NULL
		it <- vector("list",5)
		for (i in 1:length(pm@model)) {
			if (length(pm@items[[i]][pm@items[[i]]%in%items])) {
				mod <- c(mod, pm@model[i])
				for (j in 1:n) {
					if (items[j]%in%pm@items[[i]]) it[[i]] <- c(it[[i]], j)
				}
			}
		}
		flag <- NULL
		for (i in 1:5) {
			if (length(it[[i]])==0) flag <- c(flag,i)
		}
		it <- it[-flag]
		tmp@poly.mod <- as.poly.mod(n,mod,it)
		x <- sep.pars(tmp)
		item.nms <- items
		items <- seq(1,length(x@cat))
	}
	
	
	pars <- list(a=x@a, b=x@b, c=x@c)
	dots <- list(...)
	cat <- x@cat
	
	mod <- x@model
	dimensions <- x@dimensions
	
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
	
	##   Initialize object to hold the response probabilities
	p <- NULL
	
	##   Loop through all of the item response models and compute response probabilities
	for (i in 1:length(mod)) {
		if (mod[i]=="drm") tmp <- suppressWarnings(plink::drm(x, theta=theta, dimensions=dimensions, items=items, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="gpcm") tmp <- suppressWarnings(plink::gpcm(x, theta=theta, dimensions=dimensions, items=items, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="grm") tmp <- suppressWarnings(plink::grm(x, theta=theta, dimensions=dimensions, items=items, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="mcm") tmp <- suppressWarnings(plink::mcm(x, theta=theta, dimensions=dimensions, items=items, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="nrm") tmp <- suppressWarnings(plink::nrm(x, theta=theta, dimensions=dimensions, items=items, ...)@prob[,-c(1:dimensions)])
		p <- cbind(p,as.matrix(tmp))
	}
	
	
	##   When the argument incorrect is included, two columns of
	##   probabilities are created for each dichotomous items. The 
	##   number of categories (for p.cat in the output) needs to be
	##   adjusted accordingly
	if ("drm"%in%mod) {
		if (length(dots$incorrect)) {
			if (dots$incorrect==FALSE) cat[x@items$drm] <- 1
		} else {
			cat[x@items$drm] <- 1
		}
	}
	
	##   When the argument catprob is FALSE (i.e., for cumulative
	##   probabilities for the graded response model), the number
	##   of categories (for p.cat in the output) needs to be
	##   adjusted accordingly
	if ("grm"%in%mod) {
		if (length(dots$catprob)) {
			if (dots$catprob==FALSE)  cat[x@items$grm] <- cat[x@items$grm]-1
		} else {
			cat[x@items$grm] <- cat[x@items$grm]-1
		}
	}
	
	##   Based on the included models, determine the values that need to be returned for D
	if (length(dots$D)==0) {
		D <- NULL
		if (length(dots$D.drm)) D <- c(D, D.drm=dots$D.drm)
		if (length(dots$D.gpcm)) D <- c(D, D.gpcm=dots$D.gpcm)
		if (length(dots$D.grm)) D <- c(D, D.grm=dots$D.grm)
		if (is.null(D)) D <- 1
	} else {
		D <- dots$D
	}
	
	##   Resort all of the items and categories using the original sort order
	sort <- unlist(x@items)
	cat1 <- cat[sort]
	
	##   Sort the columns in p so they correspond to the original ordering of the items
	sort <- rep(sort,cat1)
	p <- p[,order(sort)]
	
	##   Label the items and categories
	lab=NULL
	for(i in 1:length(cat1)) {
		lab=c(lab,paste("item_",item.nms[i],".",seq(1,cat[i]),sep=""))
	}
	
	##   If there is only a single polytomous item, p needs to be transposed
	if (is.vector(p) & length(cat)>1) p <- t(p)
	
	p <- data.frame(p)
	names(p) <- lab
	p <- cbind(theta,p)
	
	p <- new("irt.prob", p.cat=cat, prob=p, mod.lab=x@mod.lab, dimensions=dimensions, D=D, pars=pars, model=x@model, items=x@items)
	return(p)
	
})