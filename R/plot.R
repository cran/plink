plot.irt.prob <- function(x, y, ..., combine=NULL, item.names=NULL, item.lab=TRUE, panels=20) {
	require(lattice)
	options(graphics.record=TRUE)
	if (exists(".SavedPlots")) rm(.SavedPlots,envir=.GlobalEnv)
	if (item.lab==TRUE) strip <- strip.custom(bg="lightblue") else strip <- FALSE
	
	if (length(combine)) {
		cat <- combine 
		if (sum(combine)<sum(x@p.cat)) {
			warning("{combine} did not identify all items. The specified subset will be plotted.")
			x@prob <- x@prob[,1:(sum(cat)+1)]  
		} else if (sum(combine)>sum(x@p.cat)) {
			warning("{combine} identified too many items. The original item categories will be plotted.")
			cat <- x@p.cat
		} 
	} else {
		cat <- x@p.cat
	}
	theta <- x@prob$theta
	nt <- length(theta)
	ni <- length(cat)
	if (ncol(x@prob)>2) {
		sx <- stack(x@prob[,-1]) 
	} else {
		sx <- data.frame(values=x@prob[,-1],ind=factor(rep(colnames(x@prob)[2],nt)))
	}
	id <- NULL
	cid <- NULL
	
	for (i in 1:ni) {
		tmp <- rep(i,nt*cat[i])
		id <- c(id,tmp)
		for (j in 1:cat[i]) {
			tmpc <- rep(j,nt)
			cid <- c(cid, tmpc)
		}
	}
	if (missing(item.names)) {
		id <- factor(id,seq(1:ni),paste("Item",1:ni))
	} else { 
		id <- factor(id,seq(1:ni),item.names)
	}
	out <- cbind(rep(theta,sum(cat)),id,cid,sx)
	colnames(out)[1] <- "theta"
	
	if (!is.null(panels)) {
		if (ni>panels) {
			cat("Use PgUp and PgDn to view different plot pages\n")
			xyplot(values~theta|id,out,type="l",as.table=TRUE,ylab="Probability",xlab="Theta",groups=cid,par.strip.text=list(cex=0.7),strip=strip,layout=c(0,panels),...)
		} else {
			xyplot(values~theta|id,out,type="l",as.table=TRUE,ylab="Probability",xlab="Theta",groups=cid,par.strip.text=list(cex=0.7),strip=strip,...)
		}
	} else {
		xyplot(values~theta|id,out,type="l",as.table=TRUE,ylab="Probability",xlab="Theta",groups=cid,par.strip.text=list(cex=0.7),strip=strip,...)
	}
}

plot.sep.pars <- function(x, y, ...) {
	x <- mixed(x, ...)
	plot(x, ...)
}

plot.irt.pars <- function(x, y, ...) {
	tmp <- mixed(x, ...)
	if (x@groups>1) {
		warning("There is more than one group in {x}. No plots were produced.")
	} else {
		plot(tmp, ...)
	}
}

