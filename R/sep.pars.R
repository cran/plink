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
	
	##   Loop through all groups. In this scenario, a list of {sep.pars} objects will be returned
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
	
	colnames(x) <- NULL
	rownames(x) <- NULL
	
	##   Number of items
	ni <- nrow(x)
	
	##   If missing {cat}, assume that all items are dichotomous
	if (missing(cat)) cat <- rep(2,ni)
	
	##   If missing {poly.mod}, assume that all items are dichotomous
	if (missing(poly.mod)) poly.mod <- as.poly.mod(ni)
	
	##   Create a sorting vector for ordering the items after separating out the parameters
	sort <- NULL
	
	##  Check to see if any dichotomous items were improperly specified as polytomous items
	if ("grm" %in% poly.mod@model) {
		tmp.cat <- cat[poly.mod@items$grm]
		if (length(tmp.cat[tmp.cat==2])) {
			if ("drm" %in% poly.mod@model) {
				poly.mod@items$drm <- c(poly.mod@items$drm, poly.mod@items$grm[tmp.cat==2])
			} else {
				poly.mod@model <- c(poly.mod@model, "drm")
				poly.mod@items$drm <- poly.mod@items$grm[tmp.cat==2]
			}
			poly.mod@items$drm <- poly.mod@items$drm[order(poly.mod@items$drm)]
			poly.mod@items$grm <- poly.mod@items$grm[-poly.mod@items$grm[tmp.cat==2]]
			if (length(poly.mod@items$grm)==0) {
				poly.mod@model <- poly.mod@model[poly.mod@model!="grm"]
				poly.mod@items$grm <- NULL
			}
			cat("One or more {grm} items was specified with 2 response categories.\n")
			cat("These items were re-specified as dichotomous items\n")
		}
	}
	
	if ("gpcm" %in% poly.mod@model) {
		tmp.cat <- cat[poly.mod@items$gpcm]
		if (length(tmp.cat[tmp.cat==2])) {
			if ("drm" %in% poly.mod@model) {
				poly.mod@items$drm <- c(poly.mod@items$drm, poly.mod@items$gpcm[tmp.cat==2])
			} else {
				poly.mod@model <- c(poly.mod@model, "drm")
				poly.mod@items$drm <- poly.mod@items$gpcm[tmp.cat==2]
			}
			poly.mod@items$drm <- poly.mod@items$drm[order(poly.mod@items$drm)]
			poly.mod@items$gpcm <- poly.mod@items$gpcm[-poly.mod@items$gpcm[tmp.cat==2]]
			if (length(poly.mod@items$gpcm)==0) {
				poly.mod@model <- poly.mod@model[poly.mod@model!="gpcm"]
				poly.mod@items$gpcm <- NULL
			}
			cat("One or more {gpcm} items was specified with 2 response categories.\n")
			cat("These items were re-specified as dichotomous items\n")
		}
	}
	
	##   Number of dichotomous items
	ndi <- length(poly.mod@items$drm)
	
	##   Extract the dichotomous items
	dichot <- as.matrix(x[poly.mod@items$drm,])
	if (ndi==1) dichot <- t(dichot)
	
	if (ndi>0) {
		##   Identify the dichotomous items for re-ordering all of the items later
		sort <- c(sort,poly.mod@items$drm)
		
		##   Number of dichotomous parameters. This is used
		##   to determine the specific item response model
		ncd <- length(dichot[1,!is.na(dichot[1,])])
		
		##   When there is only a difficulty parameter
		if (ncd==1) {
			da <- matrix(1,ndi,dimensions) 
			db <- as.matrix(dichot[,1])
			dc <- matrix(0,ndi,1) 
			if (dimensions==1) dmod <- "Rasch" else dmod <- "Multidimensional Rasch"
		
		##   When there are two parameters this could be a 1PL or 2PL model
		} else if (ncd==dimensions+1) {
			da <- as.matrix(dichot[,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(dichot[,dimensions+1])
			dc <- matrix(0,ndi,1) 
			if (length(da[da==da[1]])==ndi*dimensions) { 
				if (da[1]==1) {
					if (dimensions==1) dmod <- "Rasch" else dmod <- "Multidimensional Rasch"
				} else {
					if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
				}
			} else {
				if (dimensions==1) dmod <- "2PL" else dmod <- "M2PL"
			}
		
		##   When there are three parameters this could be a 1PL, 2PL, or 3PL model
		} else if (ncd==dimensions+2) {
			da <- as.matrix(dichot[,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(dichot[,dimensions+1])
			dc <- as.matrix(dichot[,dimensions+2])
			if (length(da[da==da[1]])==ndi*dimensions) { 
				if (length(dc[dc==0])==ndi) { 
					if (da[1]==1) {
						if (dimensions==1) dmod <- "Rasch" else dmod <- "Multidimensional Rasch"
					} else {
						if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
					}
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
	
	##   If there are no dichotomous items set the following object equal to NULL
	##   This is necessary for compiling the outpur
	} else {
		dmod <- NULL
	}

	
	##   Number of polytomous items
	npi <- ni-ndi
	
	##   Extract the polytomous items
	if (npi>0) {
		##   Vector identifying all of the item response models used
		mod <- poly.mod@model
		
		##   List identifying the items associated with each model
		mod.it <- poly.mod@items
		
		##   Object to hold the separated item parameters
		##   for each polytomous model
		ppars <- vector("list",length(mod[mod!="drm"]))
		
		##   Objects used to identify the maximum number of
		##   columns in each parameter matrix across models.
		##   This is necessary when combining the parameters 
		##   over multiple models
		max.a <- max.b <- max.c <- 0 
		
		##   Initialize a vector identified the number of items
		##   associated with each polytomous model
		pnum <- NULL
		
		##   Initialize a vector identifying the polytomous models
		mods <- NULL
		
		##   Initialize a vector for the descriptive names of the polytomous models
		pnames <- NULL
		
		##   Loop through all of the polytomous models
		for (i in 1:length(mod)) {
		
			##   Skip any dichotomous items
			if (mod[i]=="drm") next
			
			##   Items associated with the given model
			pit <- mod.it[[i]]
			
			##   Number of items associated with the given model
			np <- length(pit)
			
			##   Response categories associated with the given set of items
			pcat <- cat[pit]
			
			##   Extract the item parameters for the given items
			##   associated with the given model
			poly <- as.matrix(x[pit,])
			
			##   Transpose this matrix of parameters if there is only 
			##   one item associated with the given model
			if (np==1) poly <- t(poly)
			
			##   Total number of columns in the extracted matrix of parameters
			nc <- ncol(poly)
			
			##   Maximum number of response categories across items for this given model
			mc <- max(pcat)
			
			##   Separate the parameters for the nominal response model
			if (mod[i]=="nrm") {
			
				##   Identify the NRM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				pa <-poly[,1:(mc*dimensions)]
				pb <- poly[,(mc*dimensions+1):(mc*(dimensions+1))]
				
				##   If there is only one NRM item, transpose the parameter matrices
				if (np==1) {
					pa <- t(pa)
					pb <- t(pb)
				}
				if (dimensions==1) pmod <- "Nominal Response Model" else pmod <- "Multidimensional Nominal Response Model"
				
			##   Separate the parameters for the multiple-choice model
			} else if (mod[i]=="mcm") {
				
				##   Identify the MCM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				pa <- poly[,1:(mc*dimensions)]
				pb <- poly[,(mc*dimensions+1):(mc*(dimensions+1))]
				pc <- poly[,(mc*(dimensions+1)+1):(mc*(dimensions+1)+mc-1)]
				
				##   If there is only one MCM item, transpose the parameter matrices
				if (np==1) {
					pa <- t(pa)
					pb <- t(pb)
					pc <- t(pc)
				}
				if (dimensions==1) pmod <- "Multiple-Choice Model" else pmod <- "Multidimensional Multiple-Choice Model"
				
			##   Separate the parameters for the graded response model
			} else if (mod[i]=="grm") {
				
				##   Identify the GRM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				pa <- as.matrix(poly[,1:dimensions])
				if (np==1) pa <- t(pa)
				
				if (location==FALSE) {
					pb <- poly[,(dimensions+1):(dimensions+mc-1)]
					if (np==1) pb <- t(pb)
					
				} else {
					pb <- matrix(poly[,dimensions+1],np,mc-1)+poly[,(dimensions+2):(dimensions+mc)]
				}
				
				if (dimensions==1) pmod <- "Graded Response Model" else pmod <- "Multidimensional Graded Response Model"
				
			##   Separate the parameters for the generalized partial credit model
			} else if (mod[i]=="gpcm") {
			
				##   Identify the GPCM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				##   Identify the number of parameters in the first row of the 
				##   item parameter matrix. This information is used to distinguish
				##   between the PCM and GPCM
				len.p <- length(poly[1,][!is.na(poly[1,])])
				
				if (location==FALSE) {
					
					if (len.p==pcat[1]-1) {
						pa <- matrix(1,np,dimensions)
						pb <- poly[,1:(mc-1)]
						if (np==1) pb <- t(pb)
						if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
						
					} else if (len.p==pcat[1]-1+dimensions) {
						pa <- as.matrix(poly[,1:dimensions])
						pb <- as.matrix(poly[,(dimensions+1):(dimensions+mc-1)])
						
						##   If there is only one GPCM item, transpose the parameter matrices
						if (np==1) {
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
						if (np==1) pb <- t(pb)
						if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
						
					} else if (len.p==pcat[1]+dimensions) {
						pa <- as.matrix(poly[,1:dimensions])
						if (np==1) pa <- t(pa)
						pb <- matrix(poly[,dimensions+1],np,mc-1)+poly[,(dimensions+2):(dimensions+mc)]
						
						if (length(pa[pa==pa[1]])==np*dimensions) {         
							if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
						} else {
							if (dimensions==1) pmod <- "Generalized Partial Credit Model" else pmod <- "Multidimensional Generalized Partial Credit Model"
						}
					} 
				}
				
			}
			
			##   Add NA values for the c parameters for all polytomous
			##   models with the exception of the MCM
			if (mod[i]!="mcm") pc <- matrix(NA,np,1)
			
			##   For the GRM and GPCM, reformat the threshold/step parameters 
			##   to include a location parameter (if specified using loc.out=TRUE)
			if (mod[i]=="gpcm"|mod[i]=="grm") {
				if (loc.out==TRUE) {
					pbm <- apply(pb,1,mean,na.rm=TRUE)
					pb <- cbind(pbm,pb-pbm)
				}
			}
			
			##   Update the max.a, max.b, and max.c columns (if necessary)
			if (ncol(pa)>max.a) max.a <- ncol(pa)
			if (ncol(pb)>max.b) max.b <- ncol(pb)
			if (ncol(pc)>max.c) max.c <- ncol(pc)
			
			##   Compile the separated item parameters for the given model
			ppars[[i]] <- list(a=pa,b=pb,c=pc)
			
			##   Update these objects for the given polytomous model
			pnum <- c(pnum,np)
			mods <- c(mods,pmod)
			pnames <- c(pnames,paste("Poly.",mod[i],sep=""))
			
		} 
		
		
		pa <- pb <- pc <- NULL
		pmod <- mods
		
		##   Compile all of the item parameters from the polytomous
		##   models into parameter-specific matrices. Loop through
		##   all of the polytomous models
		for(i in 1:length(mod)) {
		
			if (mod[i]=="drm") next
			
			##   Create temporary objects with the separated item
			##   parameters for the given model
			a <- ppars[[i]]$a
			b <- ppars[[i]]$b
			c <- ppars[[i]]$c
			
			##   Add columns of NAs to the right of the parameter blocks
			##   for the given model if necessary
			tmp.a <- cbind(a,matrix(NA,nrow(a),max.a-ncol(a)))
			tmp.b <- cbind(b,matrix(NA,nrow(b),max.b-ncol(b)))
			tmp.c <- cbind(c,matrix(NA,nrow(c),max.c-ncol(c)))
			
			##   Stack the blocks (for each parameter) for each model
			##   on top of one another
			pa <- rbind(pa,tmp.a)
			pb <- rbind(pb,tmp.b)
			pc <- rbind(pc,tmp.c)
		}
	
	##   If there are no polytomous items set the following objects equal to NULL
	##   This is necessary for compiling the outpur
	} else {
		pmod <- NULL
		ppars <- NULL
	}
	
	
	##   If there are a combination of  dichotomous and polytomous items
	if (ndi>0) {
		if (npi>0) {
			##   Combine the discrimination/slope parameters
			da <- cbind(da,matrix(NA,ndi,ncol(pa)-ncol(da)))
			colnames(da) <- colnames(pa) <- paste("a",1:ncol(da),sep="")
			a <- rbind(da,pa)
			
			##   Combine the difficulty/threshold/step/category  parameters
			db <- cbind(db,matrix(NA,ndi,ncol(pb)-ncol(db)))
			colnames(db) <- colnames(pb) <- paste("b",1:ncol(db),sep="")
			b <- rbind(db,pb)
			
			##   Combine the lower asymptote  parameters
			dc <- cbind(dc,matrix(NA,length(dc),ncol(pc)-ncol(dc)))
			colnames(dc) <- colnames(pc) <- paste("c",1:ncol(dc),sep="")
			c <- rbind(dc,pc)
			
		##   If there are only dichotomous items
		} else {
			a <- da
			b <- db
			c <- dc
			colnames(a) <- paste("a",1:ncol(da),sep="")
			colnames(b) <- paste("b",1:ncol(db),sep="")
			colnames(c) <- paste("c",1:ncol(dc),sep="")
		}
		
	##   If there are only polytomous items
	} else {
		a <- pa
		b <- pb
		c <- pc
		colnames(a) <- paste("a",1:ncol(pa),sep="")
		colnames(b) <- paste("b",1:ncol(pb),sep="")
		colnames(c) <- paste("c",1:ncol(pc),sep="")
	}
	
	##   Sort the parameters according to their original order
	a <- as.matrix(a[order(sort),])
	b <- as.matrix(b[order(sort),])
	c <- as.matrix(c[order(sort),])
	
	##   If there is only one item, transpose these matrices
	if (ni==1) {
		a <- t(a)
		b <- t(b)
		c <- t(c)
	}
	
	##   Compile the numbers of items associated with the various
	##   models for use in the output
	if (length(ppars)>1) {
		n <- c(ni,ndi,npi,pnum)
		names(n) <- c("Total","Dichot","Poly.all",pnames)
	} else {
		n <- c(ni,ndi,npi)
		names(n) <- c("Total","Dichot","Poly")
	}
	
	sep <- new("sep.pars", a=a, b=b, c=c, cat=cat, mod.lab=c(dmod,pmod), model=poly.mod@model, items=poly.mod@items, location=loc.out, n=n, dimensions=dimensions)
	return(sep)
})



setMethod("sep.pars", signature(x="list"), function(x, cat, poly.mod, dimensions, location, loc.out, ...) {
	
	if (dimensions<1) stop("You must specify at least one dimension")
	
	##   Format the item parameters in each list element as a matrix
	for (i in 1:length(x)) {
		x[[i]] <- as.matrix(x[[i]])
	}
	
	##   Total number of items
	ni <- nrow(x[[1]])
	
	##   If missing {cat}, assume that all items are dichotomous
	if (missing(cat)) cat <- rep(2,ni)
	
	##   If missing {poly.mod}, assume that all items are dichotomous
	if (missing(poly.mod)) poly.mod <- as.poly.mod(ni)
	
	##   Create a sorting vector for ordering the items after separating out the parameters
	sort <- NULL
	
	##  Check to see if any dichotomous items were improperly specified as polytomous items
	if ("grm" %in% poly.mod@model) {
		tmp.cat <- cat[poly.mod@items$grm]
		if (length(tmp.cat[tmp.cat==2])) {
			if ("drm" %in% poly.mod@model) {
				poly.mod@items$drm <- c(poly.mod@items$drm, poly.mod@items$grm[tmp.cat==2])
			} else {
				poly.mod@model <- c(poly.mod@model, "drm")
				poly.mod@items$drm <- poly.mod@items$grm[tmp.cat==2]
			}
			poly.mod@items$drm <- poly.mod@items$drm[order(poly.mod@items$drm)]
			poly.mod@items$grm <- poly.mod@items$grm[-poly.mod@items$grm[tmp.cat==2]]
			if (length(poly.mod@items$grm)==0) {
				poly.mod@model <- poly.mod@model[poly.mod@model!="grm"]
				poly.mod@items$grm <- NULL
			}
			cat("One or more {grm} items was specified with 2 response categories.\n")
			cat("These items were re-specified as dichotomous items\n")
		}
	}
	
	if ("gpcm" %in% poly.mod@model) {
		tmp.cat <- cat[poly.mod@items$gpcm]
		if (length(tmp.cat[tmp.cat==2])) {
			if ("drm" %in% poly.mod@model) {
				poly.mod@items$drm <- c(poly.mod@items$drm, poly.mod@items$gpcm[tmp.cat==2])
			} else {
				poly.mod@model <- c(poly.mod@model, "drm")
				poly.mod@items$drm <- poly.mod@items$gpcm[tmp.cat==2]
			}
			poly.mod@items$drm <- poly.mod@items$drm[order(poly.mod@items$drm)]
			poly.mod@items$gpcm <- poly.mod@items$gpcm[-poly.mod@items$gpcm[tmp.cat==2]]
			if (length(poly.mod@items$gpcm)==0) {
				poly.mod@model <- poly.mod@model[poly.mod@model!="gpcm"]
				poly.mod@items$gpcm <- NULL
			}
			cat("One or more {gpcm} items was specified with 2 response categories.\n")
			cat("These items were re-specified as dichotomous items\n")
		}
	}
	
	##   Number of dichotomous items
	ndi <- length(poly.mod@items$drm)
	
	if (ndi>0) {
		##   Identify the GPCM items for re-ordering all of the items later
		sort <- c(sort,poly.mod@items$drm)
		
		##   Extract dichotomous items
		tmp <- poly.mod@items$drm
		
		##   When there is only a difficulty parameter
		if (length(x)==1) {
			da <- matrix(1,ndi,dimensions)
			db <- as.matrix(x[[1]][tmp,1])
			dc <- matrix(0,ndi,1)
			if (dimensions==1) dmod <- "Rasch" else dmod <- "Multidimensional Rasch"
		
		##   When there are two parameters this could be a 1PL or 2PL model
		} else if (length(x)==2) {
			da <- as.matrix(x[[1]][tmp,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(x[[2]][tmp,1])
			dc <- matrix(0,ndi,1)
			if (length(da[da==da[1]])==ndi*dimensions) { 
				if (da[1]==1) {
					if (dimensions==1) dmod <- "Rasch" else dmod <- "Multidimensional Rasch"
				} else {
					if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
				}
			} else {
				if (dimensions==1) dmod <- "2PL" else dmod <- "M2PL"
			}
			
		##   When there are three parameters this could be a 1PL, 2PL, or 3PL model
		} else {
			da <- as.matrix(x[[1]][tmp,1:dimensions])
			if (ndi==1) da <- t(da)
			db <- as.matrix(x[[2]][tmp,1])
			dc <- as.matrix(x[[3]][tmp,1])
			if (length(da[da==da[1]])==ndi*dimensions) { 
				if (length(dc[dc==0])==ndi) { 
					if (da[1]==1) {
						if (dimensions==1) dmod <- "Rasch" else dmod <- "Multidimensional Rasch"
					} else {
						if (dimensions==1) dmod <- "1PL" else dmod <- "M1PL"
					}
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
		
	##   If there are no dichotomous items set the following object equal to NULL
	##   This is necessary for compiling the outpur
	} else {
		dmod <- NULL
	}
	
	##   Number of polytomous items
	npi <- ni-ndi
	
	##   Extract the polytomous items
	if (npi>0) {
		##   Vector identifying all of the item response models used
		mod <- poly.mod@model
		
		##   List identifying the items associated with each model
		mod.it <- poly.mod@items
		
		##   Object to hold the separated item parameters
		##   for each polytomous model
		ppars <- vector("list",length(mod[mod!="drm"]))
		
		##   Objects used to identify the maximum number of
		##   columns in each parameter matrix across models.
		##   This is necessary when combining the parameters 
		##   over multiple models
		max.a <- max.b <- max.c <- 0 
		
		##   Initialize a vector identified the number of items
		##   associated with each polytomous model
		pnum <- NULL
		
		##   Initialize a vector identifying the polytomous models
		mods <- NULL
		
		##   Initialize a vector for the descriptive names of the polytomous models
		pnames <- NULL
		
		
		##   Loop through all of the polytomous models
		for (i in 1:length(mod)) {
		
			##   Skip any dichotomous items
			if (mod[i]=="drm") next
			
			##   Items associated with the given model
			pit <- mod.it[[i]]
			
			##   Number of items associated with the given model
			np <- length(pit)
			
			##   Separate the parameters for the nominal response model
			if (mod[i]=="nrm") {
				
				##   Identify the NRM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				pa <- x[[1]][pit,]
				pb <- x[[2]][pit,]
				
				##   If there is only one NRM item, transpose the parameter matrices
				if (np==1) {
					pa <- t(pa)
					pb <- t(pb)
				}
				if (dimensions==1) pmod <- "Nominal Response Model" else pmod <- "Multidimensional Nominal Response Model"
				
				
			##   Separate the parameters for the multiple-choice model
			} else if (mod[i]=="mcm") {
				
				##   Identify the MCM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				pa <- x[[1]][pit,]
				pb <- x[[2]][pit,]
				pc <- x[[3]][pit,]
				
				##   If there is only one MCM item, transpose the parameter matrices
				if (np==1) {
					pa <- t(pa)
					pb <- t(pb)
					pc <- t(pc)
				}
				if (dimensions==1) pmod <- "Multiple-Choice Model" else pmod <- "Multidimensional Multiple-Choice Model"
				
				
			##   Separate the parameters for the graded response model
			} else if (mod[i]=="grm"){
			
				##   Identify the GRM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				pa <- as.matrix(x[[1]][pit,1:dimensions])
				pb <- as.matrix(x[[2]][pit,])
				
				##   If there is only one GRM item, transpose the parameter matrices
				if (np==1) {
					pa <- t(pa)
					pb <- t(pb)
				}
				
				if (location==TRUE) {
					pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
				}
				
				if (dimensions==1) pmod <- "Graded Response Model" else pmod <- "Multidimensional Graded Response Model"
				
			
			##   Separate the parameters for the generalized partial credit model
			} else if (mod[i]=="gpcm") {
			
				##   Identify the GPCM items for re-ordering all of the items later
				sort <- c(sort,pit)
				
				if (length(x)==1) {
					
					pa <- matrix(1,np,dimensions)
					pb <- as.matrix(x[[1]][pit,])
					if (np==1) pb <- t(pb)
					
					if (location==TRUE) {
						pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					
					if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
				} else {
				
					pa <- as.matrix(x[[1]][pit,1:dimensions])
					pb <- as.matrix(x[[2]][pit,])
					
					##   If there is only one GPCM item, transpose the parameter matrices
					if (np==1) {
						pa <- t(pa)
						pb <- t(pb)
					}
					
					if (location==TRUE) {
						pb <- matrix(pb[,1],nrow(pb),ncol(pb)-1)+pb[,-1]
					}
					
					if (length(pa[pa==pa[1]])==np*dimensions) {         
						if (dimensions==1) pmod <- "Partial Credit Model" else pmod <- "Multidimensional Partial Credit Model"
					} else {
						if (dimensions==1) pmod <- "Generalized Partial Credit Model" else pmod <- "Multidimensional Generalized Partial Credit Model"
					}
				}
			}
			##   Add NA values for the c parameters for all polytomous
			##   models with the exception of the MCM
			if (mod[i]!="mcm") pc <- matrix(NA,np,1)
			
			##   For the GRM and GPCM, reformat the threshold/step parameters 
			##   to include a location parameter (if specified using loc.out=TRUE)
			if (mod[i]=="gpcm"|mod[i]=="grm") {
				if (loc.out==TRUE) {
					pbm <- apply(pb,1,mean,na.rm=TRUE)
					pb <- cbind(pbm,pb-pbm)
				}
			}
			
			##   Update the max.a, max.b, and max.c columns (if necessary)
			if (ncol(pa)>max.a) max.a <- ncol(pa)
			if (ncol(pb)>max.b) max.b <- ncol(pb)
			if (ncol(pc)>max.c) max.c <- ncol(pc)
			
			##   Compile the separated item parameters for the given model
			ppars[[i]] <- list(a=pa,b=pb,c=pc)
			
			##   Update these objects for the given polytomous model
			pnum <- c(pnum,np)
			mods <- c(mods,pmod)
			pnames <- c(pnames,paste("Poly.",mod[i],sep=""))
			
		} 
		
		
		pa <- pb <- pc <- NULL
		pmod <- mods
		
		##   Compile all of the item parameters from the polytomous
		##   models into parameter-specific matrices. Loop through
		##   all of the polytomous models
		for(i in 1:length(mod)) {
		
			if (mod[i]=="drm") next
			
			##   Create temporary objects with the separated item
			##   parameters for the given model
			a <- ppars[[i]]$a
			b <- ppars[[i]]$b
			c <- ppars[[i]]$c
			
			##   Add columns of NAs to the right of the parameter blocks
			##   for the given model if necessary
			tmp.a <- cbind(a,matrix(NA,nrow(a),max.a-ncol(a)))
			tmp.b <- cbind(b,matrix(NA,nrow(b),max.b-ncol(b)))
			tmp.c <- cbind(c,matrix(NA,nrow(c),max.c-ncol(c)))
			
			##   Stack the blocks (for each parameter) for each model
			##   on top of one another
			pa <- rbind(pa,tmp.a)
			pb <- rbind(pb,tmp.b)
			pc <- rbind(pc,tmp.c)
		}
	
	##   If there are no polytomous items set the following objects equal to NULL
	##   This is necessary for compiling the outpur
	} else {
		pmod <- NULL
		ppars <- NULL
	}
	
	
	##   If there are a combination of  dichotomous and polytomous items
	if (ndi>0) {
		if (npi>0) {
			##   Combine the discrimination/slope parameters
			da <- cbind(da,matrix(NA,ndi,ncol(pa)-ncol(da)))
			colnames(da) <- colnames(pa) <- paste("a",1:ncol(da),sep="")
			a <- rbind(da,pa)
			
			##   Combine the difficulty/threshold/step/category  parameters
			db <- cbind(db,matrix(NA,ndi,ncol(pb)-ncol(db)))
			colnames(db) <- colnames(pb) <- paste("b",1:ncol(db),sep="")
			b <- rbind(db,pb)
			
			##   Combine the lower asymptote  parameters
			dc <- cbind(dc,matrix(NA,length(dc),ncol(pc)-ncol(dc)))
			colnames(dc) <- colnames(pc) <- paste("c",1:ncol(dc),sep="")
			c <- rbind(dc,pc)
			
		##   If there are only dichotomous items
		} else {
			a <- da
			b <- db
			c <- dc
			colnames(a) <- paste("a",1:ncol(da),sep="")
			colnames(b) <- paste("b",1:ncol(db),sep="")
			colnames(c) <- paste("c",1:ncol(dc),sep="")
		}
		
	##   If there are only polytomous items
	} else {
		a <- pa
		b <- pb
		c <- pc
		colnames(a) <- paste("a",1:ncol(pa),sep="")
		colnames(b) <- paste("b",1:ncol(pb),sep="")
		colnames(c) <- paste("c",1:ncol(pc),sep="")
	}
	
	##   Sort the parameters according to their original order
	a <- as.matrix(a[order(sort),])
	b <- as.matrix(b[order(sort),])
	c <- as.matrix(c[order(sort),])
	
	##   If there is only one item, transpose these matrices
	if (ni==1) {
		a <- t(a)
		b <- t(b)
		c <- t(c)
	}
	
	##   Compile the numbers of items associated with the various
	##   models for use in the output
	if (length(ppars)>1) {
		n <- c(ni,ndi,npi,pnum)
		names(n) <- c("Total","Dichot","Poly.all",pnames)
	} else {
		n <- c(ni,ndi,npi)
		names(n) <- c("Total","Dichot","Poly")
	}
	
	sep <- new("sep.pars", a=a, b=b, c=c, cat=cat, mod.lab=c(dmod,pmod), model=poly.mod@model, items=poly.mod@items, location=loc.out, n=n, dimensions=dimensions)
	return(sep)
})



combine.pars <- function(x, common, grp.names) {
	##   Check to see if the objects in {x} are of class {irt.pars} or class {sep.pars}
	for (i in 1:length(x)) {
		if (!is.sep.pars(x[[i]]) & !is.irt.pars(x[[i]])) stop(paste("list element",i,"is not an object of class {irt.pars} or class {sep.pars}"))
	}
	
	##   Make sure {common} is a list and that all list elements are matrices
	if (is.matrix(common) | is.data.frame(common)) common <- list(common)
	for (i in 1:length(common)) {
		if (is.data.frame(common[[i]])) common[[i]] <- as.matrix(common[[i]])
	} 
	
	##   Initialize objects that will store all of the components necessary
	##   to create an irt.pars object
	pars <- cat <- pm <- com <- list(NULL)
	location <- dimensions <- NULL
	
	##   Number of objects (sets of parameters) to combine
	n <- length(x)
	
	##   Loop through all the elements in {x}
	for (i in 1:n) {
		if (is.irt.pars(x[[i]])) {
			##   When there are two or more groups in {x}
			if (x[[i]]@groups>1) {
				pars <- c(pars,x[[i]]@pars)
				cat <- c(cat,x[[i]]@cat)
				pm <- c(pm,x[[i]]@poly.mod)
			
			##   When there is only one group in {x}
			} else {
				pars[[length(pars)+1]] <- x[[i]]@pars
				cat[[length(cat)+1]] <- x[[i]]@cat
				pm[[length(pm)+1]] <- x[[i]]@poly.mod
			}
			location <- c(location,x[[i]]@location)
			dimensions <- c(dimensions,x[[i]]@dimensions)
			
		} else if (is.sep.pars(x[[i]])) {
			
			##  Total number of items
			ni <- x[[i]]@n[1]
			
			##   Combine the item parameters into a single matrix for the given group
			if (length(x[[i]]@model)==1) {
				if (x[[i]]@model=="drm"|x[[i]]@model=="mcm") {
					t.pars <- cbind(x[[i]]@a,x[[i]]@b,x[[i]]@c)
				} else {
					t.pars <- cbind(x[[i]]@a,x[[i]]@b)
				}
			} else {
				##   Numbers of columns for each block of parameters
				n.a <- ncol(x[[i]]@a)
				n.b <- ncol(x[[i]]@b)
				n.c <- ncol(x[[i]]@c)
				
				##   Initialize a matrix for all of the combined item parameters
				t.pars <- matrix(NA,ni,n.a+n.b+n.c)
				
				##   Insert the block of slope/discrimination parameters
				t.pars[,1:n.a] <- x[[i]]@a
				
				##   Loop through all of the models
				for (j in 1:length(x[[i]]@model)) {
					
					##   Identify a given model
					mod <- x[[i]]@model[j]
					
					##   Items associated with the given model
					items <- x[[i]]@items[[j]]
					
					##   Insert the difficulty/threshold/step/category and lower asymptote parameters
					if (mod=="nrm"|mod=="mcm") {
						t.pars[items,(n.a+1):(n.a+n.b)] <- x[[i]]@b[items,]
						if (mod=="mcm") t.pars[items,(n.a+n.b+1):(n.a+n.b+n.c)] <- x[[i]]@c[items,]
					} else {
						t.pars[items,(x[[i]]@dimensions+1):(x[[i]]@dimensions+n.b)] <- x[[i]]@b[items,]
						if (mod=="drm") t.pars[items,x[[i]]@dimensions+2] <- x[[i]]@c[items,1]
					}
				}
				
				##   In situations where there are no c-parameters (e.g., when all items are modeled
				##   using the GRM, GPCM, or NRM), the above code will still insert a column
				##   of NAs via the . Eliminate this column.
				tmp <- t.pars[,ncol(t.pars)]
				if (length(tmp[is.na(tmp)])==ni) t.pars <- t.pars[,-ncol(t.pars)]
			}
			
			##   Create the poly.mod object for the given set of parameters
			t.pm <- as.poly.mod(ni,x[[i]]@model,x[[i]]@items)
			
			pars[[length(pars)+1]] <- t.pars
			cat[[length(cat)+1]] <- x[[i]]@cat
			pm[[length(pm)+1]] <- t.pm
			location <- c(location,x[[i]]@location)
			dimensions <- c(dimensions,x[[i]]@dimensions)
		}
	} 
	
	##   Eliminate the NULL elements 
	if (is.null(pars[[1]])) pars <- pars[-1]
	if (is.null(cat[[1]])) cat <- cat[-1]
	if (is.null(pm[[1]])) pm <- pm[-1]
	
	
	##   Number of groups
	ng <- length(pars)
	
	##   Identify the number of expected common item elements
	nc <- ng-1
	
	if (length(common)==nc) {
		##   A single matrix of common items
		com <- common
	} else {
		##   Initialize a list for the common item matrices
		com <- list(NULL)
		
		##   Loop through the number of list elements in {x}
		for (i in 1:n) {
			if (i<n) {
				if (is.irt.pars(x[[i]])) {
					##   If there is more than one group (which necessitates
					##   that there will be one or more common item matrices
					if (x[[i]]@groups>1) {
					
						##   Add the common item matices from the irt.pars
						##   object to the full set of common item matrices
						if (is.list(x[[i]]@common)) {
							com <- c(com,x[[i]]@common) 
						} else {
							com[[length(com)+1]] <- x[[i]]@common
						}
					}
				} 
				##   Add the common item matix from the {common} 
				##   argument to the full set of common item matrices
				com[[length(com)+1]] <- common[[i]]
			
			##   Check to see if the last element in {x} is an irt.pars object 
			##   (it may contain common item matrices)
			} else {
				if (is.irt.pars(x[[i]])) {
					if (x[[i]]@groups>1) {
						##   Add the common item matices from the irt.pars
						##   object to the full set of common item matrices
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
	
	out <- new("irt.pars",pars=pars,cat=cat,poly.mod=pm,common=com,location=location,groups=ng,dimensions=dimensions)
	
	##   Create group names if necessary
	if (missing(grp.names)) {
		grp.names <- paste("group",1:ng,sep="")
	} else if (length(grp.names)!=ng) {
		warning("The number of group names does not match the number of groups specified in {pars}")
		grp.names <- paste("group",1:ng,sep="")
	}
		
	names(out@pars) <- names(out@cat) <- names(out@poly.mod) <- grp.names
	
	if (ng==2) {
		if (is.list(out@common)) out@common <- out@common[[1]]
		colnames(out@common) <- c(grp.names[1],grp.names[2])
	} else if (ng>2) {
		##   Initialize a vector for names of the common item matrices
		cn <- NULL
		for (i in 1:(ng-1)) {
			colnames(out@common[[i]]) <- c(grp.names[i],grp.names[i+1])
			cn <- c(cn, paste(grp.names[i],".",grp.names[i+1],sep=""))
		}
		names(out@common) <- cn
	}
	
	return(out)
}
