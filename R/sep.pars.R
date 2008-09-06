setGeneric("sep.pars", function(x, cat, poly.mod, location=FALSE, loc.out=FALSE, ...) standardGeneric("sep.pars"))

setMethod("sep.pars", signature(x="numeric"), function(x, cat, poly.mod, location=FALSE, loc.out=FALSE, ...) {
	x <- as.matrix(x)
	callGeneric()
})

setMethod("sep.pars", signature(x="data.frame"), function(x, cat, poly.mod, location=FALSE, loc.out=FALSE, ...) {
	x <- as.matrix(x)
	callGeneric()
})
	
setMethod("sep.pars", signature(x="irt.pars"), function(x, cat, poly.mod, location=FALSE, loc.out=FALSE, ...) {
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			out[[i]] <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], x@location[[i]], loc.out, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		out <- sep.pars(x@pars, x@cat, x@poly.mod, x@location, loc.out, ...)
		return(out)
	}
})

setMethod("sep.pars", signature(x="matrix"), function(x, cat, poly.mod, location=FALSE, loc.out=FALSE, ...) {
	colnames(x) <- NULL
	rownames(x) <- NULL
	ni <- nrow(x) # Number of items
	if (missing(cat)) cat <- rep(2,ni) # Default to dichotomous items
	if (missing(poly.mod)) poly.mod <- as.poly.mod(ni)
	sort <- unlist(poly.mod@items) # Sorting vector
	
	 # Extract dichotomous items
	dichot <- x[poly.mod@items$drm,]
	ndi <- length(poly.mod@items$drm)
	if (is.vector(dichot)) {
		if (ndi==1) dichot <- t(dichot) else dichot <- as.matrix(dichot)
	}
	if (ndi>0) {
		ncd <- length(dichot[1,!is.na(dichot[1,])]) # Number of dichotomous parameters
		if (ncd==1) {
			da <- rep(1,ndi)
			db <- dichot
			dc <- rep(0,ndi)
			dmod <- "1PL"
		} else if (ncd==2) {
			da <- dichot[,1]
			db <- dichot[,2]
			dc <- rep(0,ndi)
			if (length(da[da==da[1]])==ndi) dmod <- "1PL" else dmod <- "2PL"
		} else {
			da <- dichot[,1]
			db <- dichot[,2]
			dc <- dichot[,3]
			if (length(da[da==da[1]])==ndi) {
				if (length(dc[dc==0])==ndi) dmod <- "1PL" else dmod <- "1PL with Guessing"
			} else { 
				if (length(dc[dc==0])==ndi) dmod <- "2PL" else dmod <- "3PL"
			}
		}
		if (ndi==1) {
			if (da==1) {
				if (dc==0) dmod <- "1PL" else dmod <- "1PL with Guessing"
			} else {
				if (dc==0) dmod <- "2PL" else dmod <- "3PL"
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
			pit <- mod.it[[i]]
			pcat <- cat[pit]
			poly <- x[pit,]
			if (is.vector(poly)) poly <- t(poly)
			nc <- ncol(poly)
			np <- nrow(poly)
			mc <- max(pcat)
			if (mod[i]=="nrm") {
				pa <- poly[,1:mc]
				if (is.vector(pa)) pa <- t(pa) 
				pb <- poly[,(mc+1):(2*mc)]
				if (is.vector(pb)) pb <- t(pb) 
				pmod <- "Nominal Response Model"
			} else if (mod[i]=="mcm") {
				pa <- poly[,1:mc]
				if (is.vector(pa)) pa <- t(pa) 
				pb <- poly[,(mc+1):(2*mc)]
				if (is.vector(pb)) pb <- t(pb) 
				pc <- poly[,(2*mc+1):(2*mc+mc-1)]
				if (is.vector(pc)) pc <- t(pc) 
				pmod <- "Multiple-Choice Model"
			} else {
				len.p <- length(poly[1,][!is.na(poly[1,])])
				if (location==FALSE) {
					if (len.p==pcat[1]-1) { 
						pa <- rep(1,np)
						pb <- poly[,1:(mc-1)]
						if (is.vector(pb)) pb <- t(pb) 
						if (mod[i]=="grm") {
							pmod <- "Graded Response Model with Constant Discrimination"
						} else {
							pmod <- "Partial Credit Model"
						}
					} else if (len.p==pcat[1]) { 
						pa <- poly[,1]
						pb <- poly[,2:mc]
						if (is.vector(pb)) pb <- t(pb) 
						if (mod[i]=="grm") {
							pmod <- "Graded Response Model"
						} else {
							if (length(pa[pa==pa[1]])==np) {         
								pmod <- "Partial Credit Model"
							} else {
								pmod <- "Generalized Partial Credit Model"
							}
						}
					}
				} else {
					if (len.p==pcat[1]) { 
						pa <- rep(1,np)
						pb <- matrix(poly[,1],np,mc-1)+poly[,2:mc]
						if (is.vector(pb)) pb <- t(pb) 
						if (mod[i]=="grm") {
							pmod <- "Graded Response Model with Constant Discrimination"
						} else {
							pmod <- "Partial Credit Model"
						}
					} else if (len.p==pcat[1]+1) {
						pa <- poly[,1]
						pb <- matrix(poly[,2],np,mc-1)+poly[,3:(mc+1)]
						if (is.vector(pb)) pb <- t(pb) 
						if (mod[i]=="grm") {
							pmod <- "Graded Response Model"
						} else {
							if (length(pa[pa==pa[1]])==np) {         
								pmod <- "Partial Credit Model"
							} else {
								pmod <- "Generalized Partial Credit Model"
							}
						}
					} 
				}
				if (npi==1) {
					if (mod[i]=="gpcm") {
						if (pa==1) pmod <- "Partial Credit Model" else pmod <- "Generalized Partial Credit Model"
					}
				}
			}
			if (mod[i]!="mcm") pc <- rep(NA,np)
			if (mod[i]=="gpcm"|mod[i]=="grm") {
				if (loc.out==TRUE) {
					pbm <- apply(pb,1,mean,na.rm=TRUE)
					pb <- cbind(pbm,pbm-pb)
				}
			}
			pa <- as.matrix(pa)
			pb <- as.matrix(pb)
			pc <- as.matrix(pc)
			if (ncol(pa)>max.a) max.a <- ncol(pa)
			if (ncol(pb)>max.b) max.b <- ncol(pb)
			if (ncol(pc)>max.c) max.c <- ncol(pc)
			ppars[[i]] <- list(a=pa,b=pb,c=pc)
			pnum <- c(pnum,np)
			mods <- c(mods,pmod)
			pnames <- c(pnames,paste("Poly.",mod[i],sep=""))
		} 
		pa <- pb <- pc <- NULL
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
			if (ncol(pa)>1) {
				da <- cbind(da,matrix(NA,length(da),ncol(pa)-1))
				colnames(da) <- colnames(pa) <- paste("a",1:ncol(da),sep="")
				a <- rbind(da,pa)
			} else {
				a <- c(da,pa)
			}
			
			# Compile b
			db <- cbind(db,matrix(NA,length(db),ncol(pb)-1))
			colnames(db) <- colnames(pb) <- paste("b",1:ncol(db),sep="")
			b <- rbind(db,pb)
			
			# Compile c
			if (ncol(pc)>1) {
				dc <- cbind(dc,matrix(NA,length(dc),ncol(pc)-1))
				colnames(dc) <- colnames(pc) <- paste("c",1:ncol(dc),sep="")
				c <- rbind(dc,pc)
			} else {
				c <- c(dc,pc)
			}
		} else { # Only dichotomous items
			a <- da
			b <- db
			c <- dc
		}
	} else { # Only polytomous items
		a <- pa
		b <- pb
		c <- pc
	}
	
	# Sort the parameters according to their original order
	if (is.vector(a)) a <- as.matrix(a[order(sort)]) else a <- as.matrix(a[order(sort),])
	if (nrow(a)!=ni) a <- t(a)
	if (is.vector(b))  b <- as.matrix(b[order(sort)]) else b <- as.matrix(b[order(sort),])
	if (nrow(b)!=ni) b <- t(b)
	if (is.vector(c)) c <- as.matrix(c[order(sort)]) else c <- as.matrix(c[order(sort),])
	if (nrow(c)!=ni) c <- t(c)
	if (length(ppars)>1) {
		n <- c(ni,ndi,npi,pnum)
		names(n) <- c("Total","Dichot","Poly.all",pnames)
		pmod <- mods
	} else {
		n <- c(ni,ndi,npi)
		names(n) <- c("Total","Dichot","Poly")
	}
	sep <- new("sep.pars", a=a, b=b, c=c, cat=cat, mod.lab=c(dmod,pmod), model=poly.mod@model, items=poly.mod@items, loc.out=loc.out, n=n)
	return(sep)
})

setMethod("sep.pars", signature(x="list"), function(x, cat, poly.mod, location=FALSE, loc.out=FALSE, ...) {
	if (is.vector(x[[1]])) ni <- length(x[[1]]) else ni <- nrow(x[[1]])
	if (missing(cat)) cat <- rep(2,ni) # Default to dichotomous items
	if (missing(poly.mod)) poly.mod <- as.poly.mod(ni)
	sort <- unlist(poly.mod@items) # Sorting vector
	
	 # Extract dichotomous items
	ndi <- length(poly.mod@items$drm)
	if (ndi>0) {
		if (length(x)==1) {
			if (is.vector(x[[1]])) db <- x[[1]][poly.mod@items$drm] else db <- x[[1]][poly.mod@items$drm,1]
			da <- rep(1,length(db))
			dc <- rep(0,length(db))
			dmod <- "1PL"
		} else if (length(x)==2) {
			if (is.vector(x[[1]])) da <- x[[1]][poly.mod@items$drm] else da <- x[[1]][poly.mod@items$drm,1]
			if (is.vector(x[[2]])) db <- x[[2]][poly.mod@items$drm] else db <- x[[2]][poly.mod@items$drm,1]
			dc <- rep(0,length(db))
			if (length(da[da==da[1]])==ndi) dmod <- "1PL" else dmod <- "2PL"
		} else {
			if (is.vector(x[[1]])) da <- x[[1]][poly.mod@items$drm] else da <- x[[1]][poly.mod@items$drm,1]
			if (is.vector(x[[2]])) db <- x[[2]][poly.mod@items$drm] else db <- x[[2]][poly.mod@items$drm,1]
			if (is.vector(x[[3]])) dc <- x[[3]][poly.mod@items$drm] else dc <- x[[3]][poly.mod@items$drm,1]
			if (length(da[da==da[1]])==ndi) {
				if (length(dc[dc==0])==ndi) dmod <- "1PL" else dmod <- "1PL with Guessing"
			} else { 
				if (length(dc[dc==0])==ndi) dmod <- "2PL" else dmod <- "3PL"
			}
		}
		if (ndi==1) {
			if (da==1) {
				if (dc==0) dmod <- "1PL" else dmod <- "1PL with Guessing"
			} else {
				if (dc==0) dmod <- "2PL" else dmod <- "3PL"
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
				if (is.vector(pa)) pa <- t(pa) 
				pb <- x[[2]][poly.mod@items$nrm,]
				if (is.vector(pb)) pb <- t(pb) 
				np <- nrow(pb)
				pc <- rep(NA,np)
				pmod <- "Nominal Response Model"
			} else if (mod[i]=="mcm") {
				pa <- x[[1]][poly.mod@items$mcm,]
				if (is.vector(pa)) pa <- t(pa) 
				pb <- x[[2]][poly.mod@items$mcm,]
				if (is.vector(pb)) pb <- t(pb) 
				pc <- x[[3]][poly.mod@items$mcm,]
				if (is.vector(pc)) pc <- t(pc) 
				np <- nrow(pb)
				pmod <- "Multiple-Choice Model"
			} else if(mod[i]=="grm"){
				if (length(x)==1) {
					if (is.vector(x[[1]])) pb <- x[[1]][poly.mod@items$grm] else pb <- x[[1]][poly.mod@items$grm,]
					if (location==TRUE) {
						if (ncol(pb)>1) pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					if (is.vector(pb)) pb <- t(pb) 
					np <- nrow(pb)
					pa <- rep(1,np)
					pc <- rep(NA,np)
				} else {
					if (is.vector(x[[1]])) pa <- x[[1]][poly.mod@items$grm] else pa <- x[[1]][poly.mod@items$grm,1]
					if (is.vector(x[[2]])) pb <- x[[2]][poly.mod@items$grm] else pb <- x[[2]][poly.mod@items$grm,]
					if (location==TRUE) {
						browser()
						if (ncol(pb)>1) pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					if (is.vector(pb)) pb <- t(pb) 
					np <- nrow(pb)
					pc <- rep(NA,np)
				}
				if (loc.out==TRUE) {
					pbm <- apply(pb,1,mean,na.rm=TRUE)
					pb <- cbind(pbm,pbm-pb)
				}
				pmod <- "Graded Response Model"
			} else if (mod[i]=="gpcm") {
				if (length(x)==1) {
					if (is.vector(x[[1]])) pb <- x[[1]][poly.mod@items$gpcm] else pb <- x[[1]][poly.mod@items$gpcm,]
					if (location==TRUE) {
						if (ncol(pb)>1) pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					if (is.vector(pb)) pb <- t(pb) 
					np <- nrow(pb)
					pa <- rep(1,np)
					pc <- rep(NA,np)
					pmod <- "Partial Credit Model"
				} else {
					if (is.vector(x[[1]])) pa <- x[[1]][poly.mod@items$gpcm] else pa <- x[[1]][poly.mod@items$gpcm,1]
					if (is.vector(x[[2]])) pb <- x[[2]][poly.mod@items$gpcm] else pb <- x[[2]][poly.mod@items$gpcm,]
					if (location==TRUE) {
						if (ncol(pb)>1) pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					if (is.vector(pb)) pb <- t(pb) 
					np <- nrow(pb)
					pc <- rep(NA,np)
					if (length(pa[pa==pa[1]])==np) {         
						pmod <- "Partial Credit Model"
					} else {
						pmod <- "Generalized Partial Credit Model"
					}
					if (npi==1) {
						if (pa==1) pmod <- "Partial Credit Model" else pmod <- "Generalized Partial Credit Model"
					}
				}
				if (loc.out==TRUE) {
					pbm <- apply(pb,1,mean,na.rm=TRUE)
					pb <- cbind(pbm,pbm-pb)
				}
			}
			pa <- as.matrix(pa)
			pb <- as.matrix(pb)
			pc <- as.matrix(pc)
			if (ncol(pa)>max.a) max.a <- ncol(pa)
			if (ncol(pb)>max.b) max.b <- ncol(pb)
			if (ncol(pc)>max.c) max.c <- ncol(pc)
			ppars[[i]] <- list(a=pa,b=pb,c=pc)
			pnum <- c(pnum,np)
			mods <- c(mods,pmod)
			pnames <- c(pnames,paste("Poly.",mod[i],sep=""))
		} 
		pa <- pb <- pc <- NULL
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
			if (ncol(pa)>1) {
				da <- cbind(da,matrix(NA,length(da),ncol(pa)-1))
				colnames(da) <- colnames(pa) <- paste("a",1:ncol(da),sep="")
				a <- rbind(da,pa)
			} else {
				a <- c(da,pa)
			}
			
			# Compile b
			db <- cbind(db,matrix(NA,length(db),ncol(pb)-1))
			colnames(db) <- colnames(pb) <- paste("b",1:ncol(db),sep="")
			b <- rbind(db,pb)
			
			# Compile c
			if (ncol(pc)>1) {
				dc <- cbind(dc,matrix(NA,length(dc),ncol(pc)-1))
				colnames(dc) <- colnames(pc) <- paste("c",1:ncol(dc),sep="")
				c <- rbind(dc,pc)
			} else {
				c <- c(dc,pc)
			}
		} else { # Only dichotomous items
			a <- da
			b <- db
			c <- dc
		}
	} else { # Only polytomous items
		a <- pa
		b <- pb
		c <- pc
	}
	
	# Sort the parameters according to their original order
	if (is.vector(a)) a <- as.matrix(a[order(sort)]) else a <- as.matrix(a[order(sort),])
	if (nrow(a)!=ni) a <- t(a)
	if (is.vector(b))  b <- as.matrix(b[order(sort)]) else b <- as.matrix(b[order(sort),])
	if (nrow(b)!=ni) b <- t(b)
	if (is.vector(c)) c <- as.matrix(c[order(sort)]) else c <- as.matrix(c[order(sort),])
	if (nrow(c)!=ni) c <- t(c)
	if (length(ppars)>1) {
		n <- c(ni,ndi,npi,pnum)
		names(n) <- c("Total","Dichot","Poly.all",pnames)
		pmod <- mods
	} else {
		n <- c(ni,ndi,npi)
		names(n) <- c("Total","Dichot","Poly")
	}
	sep <- new("sep.pars", a=a, b=b, c=c, cat=cat, mod.lab=c(dmod,pmod), model=poly.mod@model, items=poly.mod@items, loc.out=loc.out, n=n)
	return(sep)
})


combine.pars <- function(x, common, grp.names=NULL) {
	for (i in 1:length(x)) {
		if (!is.sep.pars(x[[i]]) & !is.irt.pars(x[[i]])) stop(paste("list element",i,"is not an object of class {irt.pars} or class {sep.pars}"))
	}
	if (is.matrix(common) | is.data.frame(common)) common <- list(common)
	for (i in 1:length(common)) {
		if (is.data.frame(common[[i]])) common[[i]] <- as.matrix(common[[i]])
	} 
	pars <- cat <- pm <- com <- list(NULL)
	location <- NULL
	n <- length(x)
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
						t.pars[items,(n.a+n.b+1):(n.a+n.b+n.c)] <- x[[i]]@c[items,]
					} else {
						t.pars[items,2:(n.b+1)] <- x[[i]]@b[items,]
						if (mod=="drm") t.pars[items,3] <- x[[i]]@c[items,1]
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
		}
	} 
	pars <- pars[-1]
	cat <- cat[-1]
	pm <- pm[-1]
	
	nc <- length(pars)-1
	if (length(common)==nc) {
		com <- common
	} else {
		c.flag <- 1
		for (i in 1:(n-1)) {
			if (is.irt.pars(x[[i]])) {
				if (!is.null(x[[i]]@common)) {
					if (is.list(x[[i]]@common)) com <- c(com,x[[i]]@common) else com[[length(com)+1]] <- x[[i]]@common
					com[[length(com)+1]] <- common[[c.flag]]
					c.flag <- c.flag+1
				} else {
					com[[length(com)+1]] <- common[[c.flag]]
					c.flag <- c.flag+1
				}
			} else if (is.sep.pars(x[[i]])) {
				com[[length(com)+1]] <- common[[c.flag]]
				c.flag <- c.flag+1
			}
		}
	}
	if (is.null(com[[1]])) com <- com[-1]
	
	out <- new("irt.pars",pars=pars,cat=cat,poly.mod=pm,common=com,location=location,groups=n)
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
