setGeneric("mixed", function(x, cat, poly.mod, theta, dimensions=1, ...) standardGeneric("mixed"))

setMethod("mixed", signature(x="numeric", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, ...) {
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
})

setMethod("mixed", signature(x="matrix", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, ...) {
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
})

setMethod("mixed", signature(x="data.frame", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, ...) {
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
})

setMethod("mixed", signature(x="list", cat="numeric"), function(x, cat, poly.mod, theta, dimensions, ...) {
	x <- sep.pars(x, cat, poly.mod, dimensions, loc.out=FALSE, ...)
	callGeneric()
})

setMethod("mixed", signature(x="irt.pars"), function(x, cat, poly.mod, theta, dimensions, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], loc.out=FALSE, dimensions=x@dimensions[i], ...)
			out[[i]] <- mixed(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, loc.out=FALSE, dimensions=x@dimensions, ...)
		callGeneric()
	}
})

setMethod("mixed", signature(x="sep.pars"), function(x, cat, poly.mod, theta, dimensions, ...) {
	if (x@loc.out==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(x@n[1], x@model, x@items)
		x <- sep.pars(pars, x@cat, pm, location=TRUE, loc.out=FALSE, x@dimensions, ...)
	}
	pars <- list(a=x@a, b=x@b, c=x@c)
	dots <- list(...)
	cat <- x@cat
	items <- seq(1,length(cat))
	mod <- x@model
	dimensions <- x@dimensions
	p <- NULL
	for (i in 1:length(mod)) {
		if (mod[i]=="drm") tmp <- suppressWarnings(drm(x, theta=theta, dimensions=dimensions, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="gpcm") tmp <- suppressWarnings(gpcm(x, theta=theta, dimensions=dimensions, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="grm") tmp <- suppressWarnings(grm(x, theta=theta, dimensions=dimensions, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="mcm") tmp <- suppressWarnings(mcm(x, theta=theta, dimensions=dimensions, ...)@prob[,-c(1:dimensions)])
		if (mod[i]=="nrm") tmp <- suppressWarnings(nrm(x, theta=theta, dimensions=dimensions, ...)@prob[,-c(1:dimensions)])
		p <- cbind(p,as.matrix(tmp))
	}
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
	
	if ("drm"%in%mod) {
		if (length(dots$incorrect)) {
			if (dots$incorrect==FALSE) cat[x@items$drm] <- 1
		} else {
			cat[x@items$drm] <- 1
		}
	}
	
	if ("grm"%in%mod) {
		if (length(dots$catprob)) {
			if (dots$catprob==FALSE)  cat[x@items$grm] <- cat[x@items$grm]-1
		} else {
			cat[x@items$grm] <- cat[x@items$grm]-1
		}
	}
	sort <- unlist(x@items)
	cat1 <- cat[sort]
	sort <- rep(sort,cat1)
	p <- p[,order(sort)]
	lab=NULL
	for(i in 1:length(cat1)) {
		lab=c(lab,paste("item_",i,".",seq(1,cat[i]),sep=""))
	}
	if (is.vector(p) & length(cat)>1) p <- t(p)
	p <- data.frame(p)
	names(p) <- lab
	p <- cbind(theta,p)
	p <- new("irt.prob", p.cat=cat, prob=p, mod.lab=x@mod.lab, dimensions=dimensions, pars=pars, model=x@model, items=x@items)
	return(p)
})