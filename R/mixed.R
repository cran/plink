setGeneric("mixed", function(x, cat, poly.mod, theta=seq(-4,4,.05), ...) standardGeneric("mixed"))

setMethod("mixed", signature(x="numeric", cat="numeric"), function(x, cat, poly.mod, theta, ...) {
	x <- sep.pars(x, cat, poly.mod, ...)
	callGeneric()
})

setMethod("mixed", signature(x="matrix", cat="numeric"), function(x, cat, poly.mod, theta, ...) {
	x <- sep.pars(x, cat, poly.mod, ...)
	callGeneric()
})

setMethod("mixed", signature(x="data.frame", cat="numeric"), function(x, cat, poly.mod, theta, ...) {
	x <- sep.pars(x, cat, poly.mod, ...)
	callGeneric()
})

setMethod("mixed", signature(x="list", cat="numeric"), function(x, cat, poly.mod, theta, ...) {
	x <- sep.pars(x, cat, poly.mod, ...)
	callGeneric()
})

setMethod("mixed", signature(x="irt.pars"), function(x, cat, poly.mod, theta, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], ...)
			out[[i]] <- mixed(tmp, ...)
		}
		names(out) <- paste("Group",1:x@groups,sep="")
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, ...)
		callGeneric()
	}
})

setMethod("mixed", signature(x="sep.pars"), function(x, cat, poly.mod, theta, ...) {
	if (x@loc.out==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(list(x@model,x@items))
		x <- sep.pars(pars, x@cat, pm, location=TRUE)
	}
	dots <- list(...)
	cat <- x@cat
	items <- seq(1,length(cat))
	mod <- x@model
	p <- NULL
	for (i in 1:length(mod)) {
		if (i==1) {
			if (mod[i]=="drm") tmp <- suppressWarnings(drm(x, ...)@prob)
			if (mod[i]=="gpcm") tmp <- suppressWarnings(gpcm(x, ...)@prob)
			if (mod[i]=="grm") tmp <- suppressWarnings(grm(x, ...)@prob)
			if (mod[i]=="mcm") tmp <- suppressWarnings(mcm(x, ...)@prob)
			if (mod[i]=="nrm") tmp <- suppressWarnings(nrm(x, ...)@prob)
		} else {
			if (mod[i]=="drm") tmp <- suppressWarnings(drm(x, ...)@prob[,-1])
			if (mod[i]=="gpcm") tmp <- suppressWarnings(gpcm(x, ...)@prob[,-1])
			if (mod[i]=="grm") tmp <- suppressWarnings(grm(x, ...)@prob[,-1])
			if (mod[i]=="mcm") tmp <- suppressWarnings(mcm(x, ...)@prob[,-1])
			if (mod[i]=="nrm") tmp <- suppressWarnings(nrm(x, ...)@prob[,-1])
		}
		p <- cbind(p,as.matrix(tmp))
	}
	if ("drm"%in%mod) {
		if (!is.null(dots$incorrect)) {
			if (dots$incorrect==FALSE) cat[x@items$drm] <- 1
		} else {
			cat[x@items$drm] <- 1
		}
	}
	
	if ("grm"%in%mod) {
		if (!is.null(dots$catprob)) {
			if (dots$catprob==FALSE)  cat[x@items$grm] <- cat[x@items$grm]-1
		} else {
			cat[x@items$grm] <- cat[x@items$grm]-1
		}
	}
	
	sort <- unlist(x@items)
	cat1 <- cat[sort]
	sort <- c(1,rep(sort,cat1)+1)
	p <- p[,order(sort)]
	lab=NULL
	for(i in 1:length(cat1)) {
		lab=c(lab,paste("item_",i,".",seq(1,cat[i]),sep=""))
	}
	p <- data.frame(p)
	colnames(p) <- c("theta",lab)
	p <- new("irt.prob", p.cat=cat, prob=p, mod.lab=x@mod.lab, model=x@model, items=x@items)
	return(p)
})
