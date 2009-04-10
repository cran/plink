setGeneric("sep.pars", function(x, cat, poly.mod, dimensions=1, location=FALSE, loc.out=FALSE, ...) standardGeneric("sep.pars"))

setMethod("sep.pars", signature(x="numeric"), function(x, cat, poly.mod, dimensions, location, loc.out, ...) {
	x <- as.matrix(x)
	callGeneric()
})

setMethod("sep.pars", signature(x="data.frame"), function(x, cat, poly.mod, dimensions, location, loc.out, ...) {
	x <- as.matrix(x)
	callGeneric()
})

setMethod("sep.pars", signature(x="irt.pars"), function(x, cat, poly.mod, dimensions, location, loc.out, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			out[[i]] <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], x@dimensions[i], x@location[[i]], loc.out, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		out <- sep.pars(x@pars, x@cat, x@poly.mod, x@dimensions, x@location, loc.out, ...)
		return(out)
	}
})

setMethod("sep.pars", signature(x="matrix"), function(x, cat, poly.mod, dimensions, location, loc.out, ...) {
	if (dimensions<1) stop("You must specify at least one dimension")
	x <- as.matrix(x)
	colnames(x) <- NULL
	rownames(x) <- NULL
	ni <- nrow(x) # Number of items
	if (missing(cat)) cat <- rep(2,ni) # Default to dichotomous items
	if (missing(poly.mod)) poly.mod <- as.poly.mod(ni)
	sort <- unlist(poly.mod@items) # Sorting vector
	
	 # Extract dichotomous items
	dichot <- as.matrix(x[poly.mod@items$drm,])
	ndi <- length(poly.mod@items$drm)
	if (ndi>0) {
		ncd <- length(dichot[1,!is.na(dichot[1,])]) # Number of dichotomous parameters
		if (ncd==1) {
			da <- matrix(1,ndi,dimensions) 
			db <- as.matrix(dichot[,1])
			dc <- matrix(0,ndi,1) 
			if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
		} else if (ncd==dimensions+1) {
			da <- as.matrix(dichot[,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(dichot[,dimensions+1])
			dc <- matrix(0,ndi,1) 
			if (length(da[da==da[1]])==ndi*dimensions) { 
				if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
			} else {
				if (dimensions==1) dmod <- "2PL" else dmod <- "M2PL"
			}
		} else if (ncd==dimensions+2) {
			da <- as.matrix(dichot[,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(dichot[,dimensions+1])
			dc <- as.matrix(dichot[,dimensions+2])
			if (length(da[da==da[1]])==ndi*dimensions) { 
				if (length(dc[dc==0])==ndi) { 
					if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
				} else {
					if (dimensions==1) dmod <- "1PL with Guessing" else dmod <- "M1PL with Guessing"
				}
			} else { 
				if (length(dc[dc==0])==ndi) { 
					if (dimensions==1) dmod <- "2PL" else dmod <- "M2PL"
				} else {
					if (dimensions==1) dmod <- "3PL" else dmod <- "M3PL"
				}
			}
		}
	} else {
		dmod <- NULL
	}

	# Extract polytomous items
	npi <- ni-ndi # Number of polytomous items
	if (npi>0) {
		mod <- poly.mod@model
		mod.it <- poly.mod@items
		ppars <- vector("list",length(mod[mod!="drm"]))
		max.a <- max.b <- max.c <- 0 
		pnum <- NULL
		mods <- NULL
		pnames <- NULL
		for (i in 1:length(mod)) {
			if (mod[i]=="drm") next
			pit <- mod.it[[i]]
			pcat <- cat[pit]
			poly <- as.matrix(x[pit,])
			if (npi==1) poly <- t(poly)
			nc <- ncol(poly)
			np <- nrow(poly)
			mc <- max(pcat)
			if (mod[i]=="nrm") {
				pa <-poly[,1:(mc*dimensions)]
				pb <- poly[,(mc*dimensions+1):(mc*(dimensions+1))]
				if (npi==1) {
					pa <- t(pa)
					pb <- t(pb)
				}
				if (dimensions==1) pmod <- "Nominal Response Model" else pmod <- "Multidimensional Nominal Response Model"
			} else if (mod[i]=="mcm") {
				pa <- poly[,1:(mc*dimensions)]
				pb <- poly[,(mc*dimensions+1):(mc*(dimensions+1))]
				pc <- poly[,(mc*(dimensions+1)+1):(mc*(dimensions+1)+mc-1)]
				if (npi==1) {
					pa <- t(pa)
					pb <- t(pb)
					pc <- t(pc)
				}
				if (dimensions==1) pmod <- "Multiple-Choice Model" else pmod <- "Multidimensional Multiple-Choice Model"
			} else if (mod[i]=="grm") {
				pa <- as.matrix(poly[,1:dimensions])
				if (npi==1) pa <- t(pa)
				if (location==FALSE) {
					pb <- poly[,(dimensions+1):(dimensions+mc-1)]
					if (npi==1) pb <- t(pb)
				} else {
					pb <- matrix(poly[,dimensions+1],np,mc-1)+poly[,(dimensions+2):(dimensions+mc)]
				}
				if (dimensions==1) pmod <- "Graded Response Model" else pmod <- "Multidimensional Graded Response Model"
			} else if (mod[i]=="gpcm") {
				len.p <- length(poly[1,][!is.na(poly[1,])])
				if (location==FALSE) {
					if (len.p==pcat[1]-1) {
						pa <- matrix(1,np,dimensions)
						pb <- poly[,1:(mc-1)]
						if (npi==1) pb <- t(pb)
						if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
					} else if (len.p==pcat[1]-1+dimensions) {
						pa <- as.matrix(poly[,1:dimensions])
						pb <- as.matrix(poly[,(dimensions+1):(dimensions+mc-1)])
						if (npi==1) {
							pa <- t(pa)
							pb <- t(pb)
						}
						if (length(pa[pa==pa[1]])==np*dimensions) {         
							if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
						} else {
							if (dimensions==1) pmod <- "Generalized Partial Credit Model" else pmod <- "Multidimensional Generalized Partial Credit Model"
						}
					}
				} else {
					if (len.p==pcat[1]) { 
						pa <- matrix(1,np,dimensions)
						pb <- matrix(poly[,1],np,mc-1)+poly[,2:mc]
						if (npi==1) pb <- t(pb)
						if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
					} else if (len.p==pcat[1]+dimensions) {
						pa <- as.matrix(poly[,1:dimensions])
						pb <- matrix(poly[,dimensions+1],np,mc-1)+poly[,(dimensions+2):(dimensions+mc)]
						if (npi==1) {
							pa <- t(pa)
							pb <- t(pb)
						}
						if (length(pa[pa==pa[1]])==np*dimensions) {         
							if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
						} else {
							if (dimensions==1) pmod <- "Generalized Partial Credit Model" else pmod <- "Multidimensional Generalized Partial Credit Model"
						}
					} 
				}
			}
			if (mod[i]!="mcm") pc <- matrix(NA,np,1)
			if (mod[i]=="gpcm"|mod[i]=="grm") {
				if (loc.out==TRUE) {
					pbm <- apply(pb,1,mean,na.rm=TRUE)
					pb <- cbind(pbm,pb-pbm)
				}
			}
			if (ncol(pa)>max.a) max.a <- ncol(pa)
			if (ncol(pb)>max.b) max.b <- ncol(pb)
			if (ncol(pc)>max.c) max.c <- ncol(pc)
			ppars[[i]] <- list(a=pa,b=pb,c=pc)
			pnum <- c(pnum,np)
			mods <- c(mods,pmod)
			pnames <- c(pnames,paste("Poly.",mod[i],sep=""))
		} 
		pa <- pb <- pc <- NULL
		pmod <- mods
		for(i in 1:length(mod)) {
			if (mod[i]=="drm") next
			a <- ppars[[i]]$a
			b <- ppars[[i]]$b
			c <- ppars[[i]]$c
			tmp.a <- cbind(a,matrix(NA,nrow(a),max.a-ncol(a)))
			tmp.b <- cbind(b,matrix(NA,nrow(b),max.b-ncol(b)))
			tmp.c <- cbind(c,matrix(NA,nrow(c),max.c-ncol(c)))
			pa <- rbind(pa,tmp.a)
			pb <- rbind(pb,tmp.b)
			pc <- rbind(pc,tmp.c)
		}
	} else {
		pmod <- NULL
		ppars <- NULL
	}
	
	# Compile dichotomous and polytomous item parameters
	if (ndi>0) {
		if (npi>0) {
			# Compile a
			da <- cbind(da,matrix(NA,ndi,ncol(pa)-ncol(da)))
			colnames(da) <- colnames(pa) <- paste("a",1:ncol(da),sep="")
			a <- rbind(da,pa)
			
			# Compile b
			db <- cbind(db,matrix(NA,ndi,ncol(pb)-ncol(db)))
			colnames(db) <- colnames(pb) <- paste("b",1:ncol(db),sep="")
			b <- rbind(db,pb)
			
			# Compile c
			dc <- cbind(dc,matrix(NA,length(dc),ncol(pc)-ncol(dc)))
			colnames(dc) <- colnames(pc) <- paste("c",1:ncol(dc),sep="")
			c <- rbind(dc,pc)
		} else { # Only dichotomous items
			a <- da
			b <- db
			c <- dc
			colnames(a) <- paste("a",1:ncol(da),sep="")
			colnames(b) <- paste("b",1:ncol(db),sep="")
			colnames(c) <- paste("c",1:ncol(dc),sep="")
		}
	} else { # Only polytomous items
		a <- pa
		b <- pb
		c <- pc
		colnames(a) <- paste("a",1:ncol(pa),sep="")
		colnames(b) <- paste("b",1:ncol(pb),sep="")
		colnames(c) <- paste("c",1:ncol(pc),sep="")
	}
	
	# Sort the parameters according to their original order
	a <- as.matrix(a[order(sort),])
	b <- as.matrix(b[order(sort),])
	c <- as.matrix(c[order(sort),])
	if (ni==1) {
		a <- t(a)
		b <- t(b)
		c <- t(c)
	}
	if (length(ppars)>1) {
		n <- c(ni,ndi,npi,pnum)
		names(n) <- c("Total","Dichot","Poly.all",pnames)
	} else {
		n <- c(ni,ndi,npi)
		names(n) <- c("Total","Dichot","Poly")
	}
	sep <- new("sep.pars", a=a, b=b, c=c, cat=cat, mod.lab=c(dmod,pmod), model=poly.mod@model, items=poly.mod@items, loc.out=loc.out, n=n, dimensions=dimensions)
	return(sep)
})

setMethod("sep.pars", signature(x="list"), function(x, cat, poly.mod, dimensions, location, loc.out, ...) {
	if (dimensions<1) stop("You must specify at least one dimension")
	for (i in 1:length(x)) {
		x[[i]] <- as.matrix(x[[i]])
	}
	ni <- nrow(x[[1]])
	if (missing(cat)) cat <- rep(2,ni) # Default to dichotomous items
	if (missing(poly.mod)) poly.mod <- as.poly.mod(ni)
	sort <- unlist(poly.mod@items) # Sorting vector
	
	 # Extract dichotomous items
	ndi <- length(poly.mod@items$drm)
	if (ndi>0) {
		tmp <- poly.mod@items$drm
		if (length(x)==1) {
			da <- matrix(1,ndi,dimensions)
			db <- as.matrix(x[[1]][tmp,1])
			dc <- matrix(0,ndi,1)
			if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
		} else if (length(x)==2) {
			da <- as.matrix(x[[1]][tmp,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(x[[2]][tmp,1])
			dc <- matrix(0,ndi,1)
			if (length(da[da==da[1]])==ndi*dimensions) {
				if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
			} else {
				if (dimensions==1) dmod <- "2PL" else dmod <- "M2PL"
			}
		} else {
			da <- as.matrix(x[[1]][tmp,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(x[[2]][tmp,1])
			dc <- as.matrix(x[[3]][tmp,1])
			if (length(da[da==da[1]])==ndi*dimensions) {
				if (length(dc[dc==0])==ndi) {
					if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
				} else {
					if (dimensions==1) dmod <- "1PL with Guessing" else dmod <- "M1PL with Guessing"
				}
			} else { 
				if (length(dc[dc==0])==ndi) {
					if (dimensions==1) dmod <- "2PL" else dmod <- "M2PL"
				} else {
					if (dimensions==1) dmod <- "3PL" else dmod <- "M3PL"
				}
			}
		}
	} else {
		dmod <- NULL
	}
	
	# Extract polytomous items
	npi <- ni-ndi
	if (npi>0) {
		mod <- poly.mod@model
		mod.it <- poly.mod@items
		ppars <- vector("list",length(mod[mod!="drm"]))
		max.a <- max.b <- max.c <- 0
		pnum <- NULL
		mods <- NULL
		pnames <- NULL
		for (i in 1:length(mod)) {
			if (mod[i]=="drm") next
			if (mod[i]=="nrm") {
				pa <- x[[1]][poly.mod@items$nrm,]
				pb <- x[[2]][poly.mod@items$nrm,]
				if (npi==1) {
					pa <- t(pa)
					pb <- t(pb)
				}
				np <- nrow(pb)
				if (dimensions==1) pmod <- "Nominal Response Model" else pmod <- "Multidimensional Nominal Response Model"
			} else if (mod[i]=="mcm") {
				pa <- x[[1]][poly.mod@items$mcm,]
				pb <- x[[2]][poly.mod@items$mcm,]
				pc <- x[[3]][poly.mod@items$mcm,]
				if (npi==1) {
					pa <- t(pa)
					pb <- t(pb)
					pc <- t(pc)
				}
				np <- nrow(pb)
				if (dimensions==1) pmod <- "Multiple-Choice Model" else pmod <- "Multidimensional Multiple-Choice Model"
			} else if(mod[i]=="grm"){
				pa <- as.matrix(x[[1]][poly.mod@items$grm,1:dimensions])
				pb <- as.matrix(x[[2]][poly.mod@items$grm,])
				if (npi==1) {
					pa <- t(pa)
					pb <- t(pb)
				}
				if (location==TRUE) {
					pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
				}
				np <- nrow(pb)
				if (dimensions==1) pmod <- "Graded Response Model" else pmod <- "Multidimensional Graded Response Model"
			} else if (mod[i]=="gpcm") {
				if (length(x)==1) {
					pb <- as.matrix(x[[1]][poly.mod@items$gpcm,])
					if (npi==1) pb <- t(pb)
					if (location==TRUE) {
						pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					np <- nrow(pb)
					pa <- matrix(1,np,dimensions)
					if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
				} else {
					pa <- as.matrix(x[[1]][poly.mod@items$gpcm,1:dimensions])
					pb <- as.matrix(x[[2]][poly.mod@items$gpcm,])
					if (npi==1) {
						pa <- t(pa)
						pb <- t(pb)
					}
					if (location==TRUE) {
						pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					np <- nrow(pb)
					if (length(pa[pa==pa[1]])==np*dimensions) {         
						if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
					} else {
						if (dimensions==1) pmod <- "Generalized Partial Credit Model" else pmod <- "Multidimensional Generalized Partial Credit Model"
					}
				}
			}
			if (mod[i]!="mcm") pc <- matrix(NA,np,1)
			if (mod[i]=="gpcm"|mod[i]=="grm") {
				if (loc.out==TRUE) {
					pbm <- apply(pb,1,mean,na.rm=TRUE)
					pb <- cbind(pbm,pb-pbm)
				}
			}
			if (ncol(pa)>max.a) max.a <- ncol(pa)
			if (ncol(pb)>max.b) max.b <- ncol(pb)
			if (ncol(pc)>max.c) max.c <- ncol(pc)
			ppars[[i]] <- list(a=pa,b=pb,c=pc)
			pnum <- c(pnum,np)
			mods <- c(mods,pmod)
			pnames <- c(pnames,paste("Poly.",mod[i],sep=""))
		} 
		pa <- pb <- pc <- NULL
		pmod <- mods
		for(i in 1:length(mod)) {
			if (mod[i]=="drm") next
			a <- ppars[[i]]$a
			b <- ppars[[i]]$b
			c <- ppars[[i]]$c
			tmp.a <- cbind(a,matrix(NA,nrow(a),max.a-ncol(a)))
			tmp.b <- cbind(b,matrix(NA,nrow(b),max.b-ncol(b)))
			tmp.c <- cbind(c,matrix(NA,nrow(c),max.c-ncol(c)))
			pa <- rbind(pa,tmp.a)
			pb <- rbind(pb,tmp.b)
			pc <- rbind(pc,tmp.c)
		}
	} else {
		pmod <- NULL
		ppars <- NULL
	}

	# Compile dichotomous and polytomous item parameters
	if (ndi>0) {
		if (npi>0) {
			# Compile a
			da <- cbind(da,matrix(NA,ndi,ncol(pa)-ncol(da)))
			colnames(da) <- colnames(pa) <- paste("a",1:ncol(da),sep="")
			a <- rbind(da,pa)
			
			# Compile b
			db <- cbind(db,matrix(NA,ndi,ncol(pb)-ncol(db)))
			colnames(db) <- colnames(pb) <- paste("b",1:ncol(db),sep="")
			b <- rbind(db,pb)
			
			# Compile c
			dc <- cbind(dc,matrix(NA,length(dc),ncol(pc)-ncol(dc)))
			colnames(dc) <- colnames(pc) <- paste("c",1:ncol(dc),sep="")
			c <- rbind(dc,pc)
		} else { # Only dichotomous items
			a <- da
			b <- db
			c <- dc
			colnames(a) <- paste("a",1:ncol(da),sep="")
			colnames(b) <- paste("b",1:ncol(db),sep="")
			colnames(c) <- paste("c",1:ncol(dc),sep="")
		}
	} else { # Only polytomous items
		a <- pa
		b <- pb
		c <- pc
		colnames(a) <- paste("a",1:ncol(pa),sep="")
		colnames(b) <- paste("b",1:ncol(pb),sep="")
		colnames(c) <- paste("c",1:ncol(pc),sep="")
	}
	
	# Sort the parameters according to their original order
	a <- as.matrix(a[order(sort),])
	b <- as.matrix(b[order(sort),])
	c <- as.matrix(c[order(sort),])
	if (ni==1) {
		a <- t(a)
		b <- t(b)
		c <- t(c)
	}
	if (length(ppars)>1) {
		n <- c(ni,ndi,npi,pnum)
		names(n) <- c("Total","Dichot","Poly.all",pnames)
	} else {
		n <- c(ni,ndi,npi)
		names(n) <- c("Total","Dichot","Poly")
	}
	sep <- new("sep.pars", a=a, b=b, c=c, cat=cat, mod.lab=c(dmod,pmod), model=poly.mod@model, items=poly.mod@items, loc.out=loc.out, n=n, dimensions=dimensions)
	return(sep)
})

combine.pars <- function(x, common, grp.names) {
	if (missing(grp.names)) grp.names <- paste("group",1:length(x),sep="")
	for (i in 1:length(x)) {
		if (!is.sep.pars(x[[i]]) & !is.irt.pars(x[[i]])) stop(paste("list element",i,"is not an object of class {irt.pars} or class {sep.pars}"))
	}
	if (is.matrix(common) | is.data.frame(common)) common <- list(common)
	for (i in 1:length(common)) {
		if (is.data.frame(common[[i]])) common[[i]] <- as.matrix(common[[i]])
	} 
	pars <- cat <- pm <- com <- list(NULL)
	location <- dimensions <- NULL
	n <- length(x) # Number of objects (sets of parameters) to combine
	for (i in 1:n) {
		if (is.irt.pars(x[[i]])) {
			if (is.list(x[[i]]@pars)) {
				pars <- c(pars,x[[i]]@pars)
				cat <- c(cat,x[[i]]@cat)
				pm <- c(pm,x[[i]]@poly.mod)
			} else {
				pars[[length(pars)+1]] <- x[[i]]@pars
				cat[[length(cat)+1]] <- x[[i]]@cat
				pm[[length(pm)+1]] <- x[[i]]@poly.mod
			}
			location <- c(location,x[[i]]@location)
			dimensions <- c(dimensions,x[[i]]@dimensions)
			
		} else if (is.sep.pars(x[[i]])) {
			ni <- x[[i]]@n[1]
			if (length(x[[i]]@model)==1) {
				if (x[[i]]@model=="drm"|x[[i]]@model=="mcm") {
					t.pars <- cbind(x[[i]]@a,x[[i]]@b,x[[i]]@c)
				} else {
					t.pars <- cbind(x[[i]]@a,x[[i]]@b)
				}
			} else {
				n.a <- ncol(x[[i]]@a)
				n.b <- ncol(x[[i]]@b)
				n.c <- ncol(x[[i]]@c)
				t.pars <- matrix(NA,ni,n.a+n.b+n.c)
				t.pars[,1:n.a] <- x[[i]]@a
				for (j in 1:length(x[[i]]@model)) {
					mod <- x[[i]]@model[j]
					items <- x[[i]]@items[[j]]
					if (mod=="nrm"|mod=="mcm") {
						t.pars[items,(n.a+1):(n.a+n.b)] <- x[[i]]@b[items,]
						if (mod=="mcm") t.pars[items,(n.a+n.b+1):(n.a+n.b+n.c)] <- x[[i]]@c[items,]
					} else {
						t.pars[items,(x[[i]]@dimensions+1):(x[[i]]@dimensions+n.b)] <- x[[i]]@b[items,]
						if (mod=="drm") t.pars[items,x[[i]]@dimensions+2] <- x[[i]]@c[items,1]
					}
				}
				tmp <- t.pars[,ncol(t.pars)]
				if (length(tmp[is.na(tmp)])==ni) t.pars <- t.pars[,-ncol(t.pars)]
			}
			t.pm <- as.poly.mod(ni,x[[i]]@model,x[[i]]@items)
			pars[[length(pars)+1]] <- t.pars
			cat[[length(cat)+1]] <- x[[i]]@cat
			pm[[length(pm)+1]] <- t.pm
			location <- c(location,x[[i]]@loc.out)
			dimensions <- c(dimensions,x[[i]]@dimensions)
		}
	} 
	# Eliminate the NULL elements 
	pars <- pars[-1]
	cat <- cat[-1]
	pm <- pm[-1]
	
	nc <- length(pars)-1 # Number of common item elements
	if (length(common)==nc) {
		com <- common
	} else {
		com <- list(NULL)
		for (i in 1:n) {
			if (i<n) {
				if (is.irt.pars(x[[i]])) {
					if (!is.null(x[[i]]@common)) {
						if (is.list(x[[i]]@common)) {
							com <- c(com,x[[i]]@common) 
						} else {
							com[[length(com)+1]] <- x[[i]]@common
						}
						com[[length(com)+1]] <- common[[i]]
					} else {
						com[[length(com)+1]] <- common[[i]]
					}
				} else if (is.sep.pars(x[[i]])) {
					com[[length(com)+1]] <- common[[i]]
				}
			} else {
				if (is.irt.pars(x[[i]])) {
					if (!is.null(x[[i]]@common)) {
						if (is.list(x[[i]]@common)) {
							com <- c(com,x[[i]]@common) 
						} else {
							com[[length(com)+1]] <- x[[i]]@common
						}
					}
				} 
			}
		}
	}
	if (is.null(com[[1]])) com <- com[-1]
	
	out <- new("irt.pars",pars=pars,cat=cat,poly.mod=pm,common=com,location=location,groups=n,dimensions=dimensions)
	if (is.null(grp.names)) {
		grp.names <- paste("group",1:n,sep="")
	} else if (length(grp.names)!=n) {
		warning("The number of group names does not match the number of groups specified in {pars}")
		grp.names <- paste("group",1:n,sep="")
	}
		
	names(out@pars) <- names(out@cat) <- names(out@poly.mod) <- grp.names
	if (n==2) {
		if (is.list(out@common)) out@common <- out@common[[1]]
		colnames(out@common) <- c(grp.names[1],grp.names[2])
	} else if (n>2) {
		cn <- NULL
		for (i in 1:(n-1)) {
			colnames(out@common[[i]]) <- c(grp.names[i],grp.names[i+1])
			cn <- c(cn, paste(grp.names[i],".",grp.names[i+1],sep=""))
		}
		names(out@common) <- cn
	}
	return(out)
}