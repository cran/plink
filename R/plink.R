setGeneric("plink", function(x, common, rescale, ability, weights, startvals, score=1, base.grp=1, symmetric=TRUE, grp.names=NULL, mn=c(FALSE,FALSE), ...) standardGeneric("plink"))

setMethod("plink", signature(x="list", common="list"), function(x, common, rescale, ability, weights, startvals, score, base.grp, symmetric, grp.names, mn, ...) {
	x <- combine.pars(x, common, ...)
	callGeneric()
})

setMethod("plink", signature(x="list", common="matrix"), function(x, common, rescale, ability, weights, startvals, score, base.grp, symmetric, grp.names, mn, ...) {
	x <- combine.pars(x, common, ...)
	callGeneric()
})

setMethod("plink", signature(x="list", common="data.frame"), function(x, common, rescale, ability, weights, startvals, score, base.grp, symmetric, grp.names, mn, ...) {
	x <- combine.pars(x, common, ...)
	callGeneric()
})

setMethod("plink", signature(x="irt.pars", common="ANY"), function(x, common, rescale, ability, weights, startvals, score, base.grp, symmetric, grp.names, mn, ...) {

	.CC <- function(startvals, to, from, weights, score, transform, symmetric, D, incorrect, catprob, mn) {
		theta1 <- weights$points[,1]
		theta2 <- weights$points[,2]
		to.f <- to
		from.t <- from
		pm <- as.poly.mod(length(to@cat),to@model,to@items)
		
		alpha <- startvals[1]
		beta <- startvals[2]
		mcnr <- NULL
		for (i in 1:length(to@model)) {
			if (to@model[i]=="mcm"|to@model[i]=="nrm") mcnr <- c(mcnr, to@items[[i]])
		}
		# Transform the item parameters for the FROM scale
		from.t@a <- as.matrix(from@a/alpha)
		from.t@b <- as.matrix(alpha * from@b + beta)
		from.t@b[mcnr,] <- from@b[mcnr,]-(beta/alpha)*from@a[mcnr,]
		# Transform the item parameters for the TO scale
		if (symmetric==TRUE) { 
			to.f@a <- as.matrix(alpha*to@a)
			to.f@b <- as.matrix((to@b-beta)/alpha)
			to.f@b[mcnr,] <- to@b[mcnr,]+(beta*to@a[mcnr,])
			prob.to.f <- .Mixed(to.f, theta1, D, incorrect, catprob)$p
			prob.from <- .Mixed(from, theta2, D, incorrect, catprob)$p
		}
			
		tmp.to <- .Mixed(to, theta1, D, incorrect, catprob)
		prob.to <- tmp.to$p
		p.cat<- tmp.to$p.cat
		prob.from.t <- .Mixed(from.t, theta2, D, incorrect, catprob)$p
		
		cat <- rep(p.cat,p.cat)
		p.mod <- rep(tmp.to$p.mod,p.cat)
		scr <- NULL
		for (h in 1:length(p.cat)) {
			scr <- c(scr,seq(1,p.cat[h]))
		}
		if (length(score)==1) {
			if (score==2) {
				scr <- scr-1
				scr[cat==1] <- 1
			}
			if (mn==FALSE) {
				scr[p.mod=="nrm"|p.mod=="mcm"] <- 0
			}
		} else {
			if (length(score)==length(scr)) {
				scr <- score
			} else {
				warning("The length of {score} does not match the number of response probabilities. Score was set to 1")
			}
		}
		scr <- matrix(scr,nrow(prob.to),length(scr),byrow=TRUE)
		# Compute sum of squares difference
		if (transform=="HB") {
			L1 <- ncol(prob.to)*sum(weights$weights[,1])
			Q1 <- weights$weights[,1]*(prob.to-prob.from.t)^2
			Q1 <- sum(apply(Q1,1,sum))
			if (symmetric==TRUE) {
				L2 <- ncol(prob.to)*sum(weights$weights[,2])
				Q2 <- weights$weights[,2]*(prob.from-prob.to.f)^2
				Q2 <- sum(apply(Q2,1,sum))
				SS <- (Q1/L1)+(Q2/L2)
			} else {
				SS <- Q1/L1
			}
		} else if (transform=="SL") {
			TCC.to <- apply(scr*prob.to, 1, sum)
			TCC.from <- apply(scr*prob.from.t, 1, sum)
			Q1 <-sum(weights$weights[,1]*(TCC.to - TCC.from)^2)
			if (symmetric==TRUE) {
				TCC.to <- apply(scr*prob.to.f, 1, sum)
				TCC.from <- apply(scr*prob.from, 1, sum)
				Q2 <-sum(weights$weights[,2]*(TCC.from- TCC.to)^2)
				SS <- Q1+Q2
			} else {
				SS <- Q1
			}
		}
	return(SS)
	}
	
	.Mixed <- function(x, theta, D, incorrect, catprob) {
		mod <- x@model
		p <- NULL
		p.cat <- NULL
		p.mod <- NULL
		for (i in 1:length(mod)) {
			if (mod[i]=="drm") tmp <- .Drm(x, theta, D, incorrect)
			if (mod[i]=="gpcm") tmp <- .Gpcm(x, theta, D)
			if (mod[i]=="grm") tmp <- .Grm(x, theta, catprob, D)
			if (mod[i]=="mcm") tmp <- .Mcm(x, theta)
			if (mod[i]=="nrm") tmp <- .Nrm(x, theta)
			p <- cbind(p, tmp$p)
			p.cat <- c(p.cat,tmp$cat)
			p.mod <- c(p.mod,tmp$mod)
		}
		return(list(p=p,p.cat=p.cat,p.mod=p.mod))
	}
	
	.Drm <- function(x, theta, D, incorrect) {
		items <- x@items$drm
		n <- nrow(as.matrix(x@a[items,]))
		a <- x@a[items,1]
		b <- x@b[items,1]
		c <- x@c[items,1]
		
		# Compute item probabilities
		p <- NULL 
		for (i in 1:length(b)) {
			cp <- c[i]+(1-c[i])/(1+exp(-D*a[i]*(theta-b[i])))
			if (incorrect==TRUE) p <- cbind(p,(1-cp),cp) else p <- cbind(p,cp)
		}
		if (incorrect==TRUE) cat <- rep(2,n) else cat <- rep(1,n)
		mod <- rep("drm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	.Gpcm <- function(x, theta, D) {
		items <- x@items$gpcm
		n <- nrow(as.matrix(x@a[items,]))
		a <- x@a[items,1]
		b <- x@b[items,]
		if (is.vector(b)) b <- t(b)
		cat <- x@cat[items]
		
		# Compute category probabilities
		p <- NULL 
		for (i in 1:nrow(b)) {
			dif <- 0 # Difference between subsequent step parameters
			den <- NULL # Compute the denominator
			for (k in 0:(cat[i]-1)) {
				if (k==1) dif <- b[i,k] else if (k>1) dif <- dif+b[i,k]
				d <- exp(D*a[i]*(k*theta-dif))
				den <- cbind(den, d)
			}
			den <- apply(den,1,sum)
			dif <- 0
			for (k in 0:(cat[i]-1)) {
				if (k==1) dif <- b[i,k] else if (k>1) dif <- dif+b[i,k]
				cp <- (exp(D*a[i]*(k*theta-dif)))/den
				p <- cbind(p,cp)
			}
		}
		mod <- rep("gpcm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	.Grm <- function(x, theta, catprob, D) {
		items <- x@items$grm
		n <- nrow(as.matrix(x@a[items,]))
		a <- x@a[items,1]
		b <- x@b[items,]
		if (is.vector(b)) b <- t(b)
		cat <- x@cat[items]
		
		p <- NULL 
		# Compute category probabilities
		if (catprob==TRUE) { 
			for (i in 1:n) {
				ct <- cat[i]-1
				cp <- 1-1/(1+exp(-D*a[i]*(theta-b[i,1])))
				p <- cbind(p, cp)
				for (k in 1:ct) {
					if (k<ct) {
						cp <- (1/(1+exp(-D*a[i]*(theta-b[i,k]))))-(1/(1+exp(-D*a[i]*(theta-b[i,k+1]))))
					} else if (k==ct) {
						cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
					}
					p <- cbind(p, cp)
				}
			}
		# Compute cumulative probabilities
		} else if (catprob==FALSE) { 
			for (i in 1:n) {
				for (k in 1:(cat[i]-1)) {
					cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
					p <- cbind(p, cp)
				}
			}
		}
		if (catprob==FALSE) cat <- cat-1
		mod <- rep("grm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	.Mcm <- function(x, theta) {
		items <- x@items$mcm
		n <- nrow(as.matrix(x@a[items,]))
		a <- x@a[items,]
		if (is.vector(a)) a <- t(a)
		b <- x@b[items,]
		if (is.vector(b)) b <- t(b)
		c <- x@c[items,]
		if (is.vector(c)) b <- t(c)
		cat <- x@cat[items]
		
		# Compute category probabilities
		p <- NULL 
		for (i in 1:n) {
			den <- NULL # Compute the denominator
			for (k in 1:cat[i]) {
				d <- exp(a[i,k]*theta+b[i,k])
				den <- cbind(den, d)
			}
			den <- apply(den,1,sum)
			for (k in 1:cat[i]) {
				if (k==1) {
					cp <- (exp(a[i,k]*theta+b[i,k]))/den
				} else {
					cp <- (exp(a[i,k]*theta+b[i,k])+c[i,(k-1)]*(exp(a[i,1]*theta+b[i,1])))/den
				}
				p <- cbind(p,cp)
			}
		}
		mod <- rep("mcm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	.Nrm <- function(x, theta) {
		items <- x@items$nrm
		n <- nrow(as.matrix(x@a[items,]))
		a <- x@a[items,]
		if (is.vector(a)) a <- t(a)
		b <- x@b[items,]
		if (is.vector(b)) b <- t(b)
		cat <- x@cat[items]
		
		# Compute category probabilities
		p <- NULL 
		for (i in 1:n) {
			den <- NULL # Compute the denominator
			for (k in 1:cat[i]) {
				d <- exp(a[i,k]*theta+b[i,k])
				den <- cbind(den, d)
			}
			den <- apply(den,1,sum)
			for (k in 1:cat[i]) {
				cp <- (exp(a[i,k]*theta+b[i,k]))/den
				p <- cbind(p,cp)
			}
		}
		mod <- rep("nrm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	.sd <- function(x) {
		z <- x[!is.na(x)]
		out <- sqrt(sum((z-mean(z,na.rm=T))^2)/length(z))
		return(out)
	}
	
	### START PLINK ###
	ng <- x@groups
	if (missing(rescale)) rescale <- NULL
	if (missing(ability)) ability <- NULL
	if (!hasArg(grp.names)) grp.names <- paste("group",1:ng,sep="")
	if (missing(weights)) weights <- as.weight()
	dots <- list(...)
	if (!is.null(dots$D)) D <- dots$D else D <- 1.7
	if (!is.null(dots$incorrect)) incorrect <- dots$incorrect else incorrect <- FALSE
	if (!is.null(dots$catprob)) catprob <- dots$catprob else catprob <- TRUE
	
	# Extract common items
	com <- vector("list",ng-1)
	catci <- vector("list",ng-1)
	pm <- vector("list",ng-1)
	for (i in 1:(ng-1)) {
		com[[i]] <- vector("list",2)
		if (ng==2) {
			it.com <- x@common
			com[[i]][[1]] <- x@pars[[i]][it.com[,1],] # Common items for the lower group
			com[[i]][[2]] <- x@pars[[i+1]][it.com[,2],] # Common items for the higher group
			catci[[i]] <- x@cat[[i]][it.com[,1]] # Category vector for the common items
		} else if (ng>2) {
			it.com <- x@common[[i]]
			com[[i]][[1]] <- x@pars[[i]][it.com[,1],] # Common items for the lower group
			com[[i]][[2]] <- x@pars[[i+1]][it.com[,2],] # Common items for the higher group
			catci[[i]] <- x@cat[[i]][it.com[,1]] # Category vector for the common items
		}
		if (i<base.grp) com[[i]] <- com[[i]][c(2,1)] # The first element is the "To" set and the second is the "From" set
		mod <-  x@poly.mod[[i]]@model
		items <- x@poly.mod[[i]]@items
		if (ng==2) it.com <- x@common[,1] else it.com <- x@common[[i]][,1]
		poly <- mod1 <- NULL
		step <- 1
	
		for (k in 1:length(mod)) {
			tmp <- seq(1,length(it.com))[it.com %in% items[[k]]]
			if (length(tmp)==0) {
				next 
			}else {
				mod1 <- c(mod1,mod[k])
				poly[[step]] <- tmp
				step <- step+1
			}
		}
		pm[[i]] <- as.poly.mod(length(it.com),mod1,poly) # Create poly.mod object for the common items
	}
	
	# Perform the calibration
	link.out <- vector("list",ng-1)
	mm <- ms <- NULL
	for (i in 1:(ng-1)) {
		tmp1 <- sep.pars(com[[i]][[1]],catci[[i]],pm[[i]],x@location[i],...)
		tmp2 <- sep.pars(com[[i]][[2]],catci[[i]],pm[[i]],x@location[i],...)
		
		a1 <- as.matrix(tmp1@a)
		a2 <- as.matrix(tmp2@a)
		b1 <- as.matrix(tmp1@b)
		b2 <- as.matrix(tmp2@b)
		c1 <- as.matrix(tmp1@c)
		c2 <- as.matrix(tmp2@c)
		
		# Transform b for mcm and nrm
		if (ncol(a1)!=ncol(b1)) {
			b1r <- as.matrix(-b1/matrix(a1,nrow(a1),ncol(b1)))
		} else {
			b1r <- as.matrix(-b1/a1)
		}
		if (ncol(a2)!=ncol(b2)) {
			b2r <- as.matrix(-b2/matrix(a2,nrow(a2),ncol(b2)))
		} else {
			b2r <- as.matrix(-b2/a2)
		}
		b1r[b1r==Inf] <- 0
		b2r[b2r==Inf] <- 0
		b1r[b1r==-Inf] <- 0
		b2r[b2r==-Inf] <- 0
		
		descrip <- NULL
		mcnr <- NULL
		pm.mod <- pm[[i]]@model
		pm.it <- pm[[i]]@items
		for (j in 1:length(pm.mod)) {
			a1k <- a1[pm.it[[j]],]
			a2k <- a2[pm.it[[j]],]
			b1k <- b1[pm.it[[j]],]
			b2k <- b2[pm.it[[j]],]
			c1k <- c1[pm.it[[j]],]
			c2k <- c2[pm.it[[j]],]
			b1kr <- b1r[pm.it[[j]],]
			b2kr <- b2r[pm.it[[j]],]
			
			if (pm.mod[j]=="mcm") {
				a <- c(length(a1k[,-1][!is.na(a1k[,-1])]),mean(a1k[,-1],na.rm=T),mean(a2k[,-1],na.rm=T),.sd(a1k[,-1]),.sd(a2k[,-1]))
				b <- c(length(b1k[,-1][!is.na(b1k[,-1])]),mean(b1k[,-1],na.rm=T),mean(b2k[,-1],na.rm=T),.sd(b1k[,-1]),.sd(b2k[,-1]))
				c <- c(length(c1k[,-1][!is.na(c1k[,-1])]),mean(c1k[,-1],na.rm=T),mean(c2k[,-1],na.rm=T),.sd(c1k[,-1]),.sd(c2k[,-1]))
				bd <- c(length(b1kr[,-1][!is.na(b1kr[,-1])]),mean(b1kr[,-1],na.rm=T),mean(b2kr[,-1],na.rm=T),.sd(b1kr[,-1]),.sd(b2kr[,-1]))
			} else {
				a <- c(length(a1k[!is.na(a1k)]),mean(a1k,na.rm=T),mean(a2k,na.rm=T),.sd(a1k),.sd(a2k))
				b <- c(length(b1k[!is.na(b1k)]),mean(b1k,na.rm=T),mean(b2k,na.rm=T),.sd(b1k),.sd(b2k))
				c <- c(length(c1k[!is.na(c1k)]),mean(c1k,na.rm=T),mean(c2k,na.rm=T),.sd(c1k),.sd(c2k))
				bd <- c(length(b1kr[!is.na(b1kr)]),mean(b1kr,na.rm=T),mean(b2kr,na.rm=T),.sd(b1kr),.sd(b2kr))
			}
			
			if (pm.mod[j]=="drm" ) {
				des <- data.frame(cbind(a,b,c))
				colnames(des) <- c("a","b","c")
			} else if (pm.mod[j]=="mcm") {
				des <- data.frame(cbind(a,b,bd,c))
				colnames(des) <- c("a","c","b (-c/a)","d")
				mcnr <- c(mcnr, pm.it[[j]])
			} else if (pm.mod[j]=="nrm") {
				des <- data.frame(cbind(a,b,bd))
				colnames(des) <- c("a","c","b (-c/a)")
				mcnr <- c(mcnr, pm.it[[j]])
			} else {
				des <- data.frame(cbind(a,b))
				colnames(des) <- c("a","b")
			}
			rownames(des) <- c("N Pars:","Mean: To","Mean: From","SD: To","SD: From")
			descrip[[j]] <- round(des,6)
		}
		
		chk.sl <- TRUE
		if (sum((mod=="nrm"|mod=="mcm")+0)==length(mod)) {
			chk.sl <- FALSE
			mn[1] <- TRUE
			cat("All items were used to compute the Mean/Mean and Mean/Sigma linking constants\n")
		}
		
		if (mn[1]==TRUE) {
			b1[mcnr,] <- b1r[mcnr,] 
			b2[mcnr,] <- b2r[mcnr,] 
			a1[mcnr,1] <- NA
			a2[mcnr,1] <- NA
			b1[mcnr,1] <- NA
			b2[mcnr,1] <- NA
		} else {
			if (!is.null(mcnr)) {
				a1 <- a1[-mcnr,]
				a2 <- a2[-mcnr,]
				b1 <- b1[-mcnr,]
				b2 <- b2[-mcnr,]
			}
		}
		a <- c(length(a1[!is.na(a1)]),mean(a1,na.rm=T),mean(a2,na.rm=T),.sd(a1),.sd(a2))
		b <- c(length(b1[!is.na(b1)]),mean(b1,na.rm=T),mean(b2,na.rm=T),.sd(b1),.sd(b2))
		desall <- round(cbind(a,b),6)
		colnames(desall) <- c("a","b")
		rownames(desall) <- c("N Pars:","Mean: To","Mean: From","SD: To","SD: From")
		descrip[[length(descrip)+1]] <- desall
		names(descrip) <- c(pm.mod,"all")
		
		A1 <- mean(a2,na.rm=TRUE)/mean(a1,na.rm=TRUE)
		A2 <- .sd(b1)/.sd(b2)
		B1 <- mean(b1,na.rm=TRUE)-A1*mean(b2,na.rm=TRUE)
		B2 <- mean(b1,na.rm=TRUE)-A2*mean(b2,na.rm=TRUE)
		mm <- round(c(A1,B1),6)
		ms <- round(c(A2,B2),6)
		names(mm) <- names(ms) <- c("Slope","Intercept")
		
		if (missing(startvals)) startvals <- c(1,0)
		if (is.list(weights[[1]][[1]])) wgt <- weights[[i]] else wgt <- weights  
		hb <- nlminb(startvals, .CC, to=tmp1, from=tmp2, weights=wgt, score=score, transform="HB", symmetric=symmetric, D=D, incorrect=incorrect, catprob=catprob, mn=mn[2])
		names(hb$par) <- c("Slope","Intercept")
		
		if (chk.sl==TRUE) {
			sl <- nlminb(startvals, .CC, to=tmp1, from=tmp2, weights=wgt, score=score, transform="SL", symmetric=symmetric, D=D, incorrect=incorrect, catprob=catprob, mn=mn[2])
			names(sl$par) <- c("Slope","Intercept")
		} else {
			sl <- list(par=NULL,iterations=NA,message=NA)
		}
		
		mm <-mm[2:1]
		ms <- ms[2:1]
		hb$par <- hb$par[2:1]
		sl$par <- sl$par[2:1]
		names(mn) <- c("Moment Methods","Characteristic Curve Methods")
		it <- c(HB=hb$iterations,SL=sl$iterations)
		con <- c(HB=hb$message,SL=sl$message)
		
		link.out[[i]] <- new("link", MM=mm, MS=ms, HB=hb$par, SL=sl$par, descriptives=descrip, iterations=it[!is.na(it)], convergence=con[!is.na(con)], base.grp=base.grp, include.mcm.nrm=mn)
		
		if (i==base.grp) {
			names(link.out)[[i]] <- paste(grp.names[i],"*/",grp.names[i+1],sep="")
		} else if (i < base.grp) {
			if ((i+1)==base.grp) {
				names(link.out)[[i]] <- paste(grp.names[i],"/",grp.names[i+1],"*",sep="")
			} else {
				names(link.out)[[i]] <- paste(grp.names[i],"/",grp.names[i+1],sep="")
			}
		} else if (i > base.grp) {
			if ((i-1)==base.grp) {
				names(link.out)[[i]] <- paste(grp.names[i],"/",grp.names[i-1],"*",sep="")
			} else {
				names(link.out)[[i]] <- paste(grp.names[i],"/",grp.names[i-1],sep="")
				
			}
		}
		
	}
	
	if (!is.null(rescale)) {
		tmp <- NULL
		for (i in 1:(ng-1)) {
			if (toupper(rescale)=="MM") tmp <- rbind(tmp, link.out[[i]]@MM)
			if (toupper(rescale)=="MS") tmp <-rbind(tmp,  link.out[[i]]@MS)
			if (toupper(rescale)=="HB") tmp <- rbind(tmp, link.out[[i]]@HB)
			if (toupper(rescale)=="SL") {
				if (is.null(link.out[[i]]@SL)) {
					tmp <- rbind(tmp, link.out[[i]]@HB)
					if (i==1) {
					warning("No linking constants were computed for the Stocking-Lord method. The linking constants from the Haebara method were used instead")
					}
				} else {
					tmp <- rbind(tmp, link.out[[i]]@SL)
				}
			}
		}
		res.meth <- matrix(NA,ng,2)
		j <- 1
		for (i in 1:ng) {
			if (i==base.grp) {
				res.meth[i,] <- c(0,1)
			} else {
				res.meth[i,] <- tmp[j,]
				j <- j+1
			}
		}
		
		out.pars <- vector("list",ng)
		out.ability <- vector("list",ng)
		for (i in 1:ng) {
			mcnr <- NULL
			pm <- x@poly.mod[[i]]@model
			for (k in 1:length(pm)) {
				if (pm[k]=="mcm"|pm[k]=="nrm") mcnr <- c(mcnr, x@poly.mod[[i]]@items[[k]])
			}
			tmp <- sep.pars(x@pars[[i]],x@cat[[i]],x@poly.mod[[i]])
			if (!is.null(ability)) tmpa <- ability[[i]]
			tmp1 <- tmp
			if (i < base.grp) {
				for (j in i:(base.grp-1)) {
					tmp1@a <- as.matrix(tmp1@a/res.meth[j,2])
					tmp1@b <-  as.matrix(res.meth[j,2]*tmp1@b + res.meth[j,1])
					tmp1@b[mcnr,] <- tmp@b[mcnr,]-(res.meth[j,1]/res.meth[j,2])*tmp@a[mcnr,]
					tmp@b[mcnr,] <- tmp1@b[mcnr,]
					if (!is.null(ability)) tmpa <- res.meth[j,2]*tmpa + res.meth[j,1]
				}
			}else if (i > base.grp) {
				for (j in i:(base.grp+1)) {
					tmp1@a <- as.matrix(tmp1@a/res.meth[j,2])
					tmp1@b <-  as.matrix(res.meth[j,2]*tmp1@b + res.meth[j,1])
					tmp1@b[mcnr,] <- tmp@b[mcnr,]-(res.meth[j,1]/res.meth[j,2])*tmp@a[mcnr,]
					tmp@b[mcnr,] <- tmp1@b[mcnr,]
					if (!is.null(ability)) tmpa <- res.meth[j,2]*tmpa + res.meth[j,1]
				}
			}
			out.pars[[i]] <- tmp1
			names(out.pars)[[i]] <- grp.names[i]
			if (!is.null(ability)) {
				out.ability[[i]] <- tmpa
				names(out.ability)[[i]] <- grp.names[i]
			}
		}
	} else {
		if (!is.null(ability)) {
			cat("Ability values will be rescaled using the Haebara linking constants.\n")
			tmp <- NULL
			for (i in 1:(ng-1)) {
				tmp <- rbind(tmp, link.out[[i]]@HB)
			}
			res.meth <- matrix(NA,ng,2)
			j <- 1
			for (i in 1:ng) {
				if (i==base.grp) {
					res.meth[i,] <- c(0,1)
				} else {
					res.meth[i,] <- tmp[j,]
					j <- j+1
				}
			}
			out.ability <- vector("list",ng)
			for (i in 1:ng) {
				tmpa <- ability[[i]]
				if (i < base.grp) {
					for (j in i:(base.grp-1)) {
						tmpa <- res.meth[j,2]*tmpa + res.meth[j,1]
					}
				}else if (i > base.grp) {
					for (j in i:(base.grp+1)) {
						tmpa <- res.meth[j,2]*tmpa + res.meth[j,1]
					}
				}
				out.ability[[i]] <- tmpa
				names(out.ability)[[i]] <- grp.names[i]
			}
		}
	}
	
	if (ng==2) {
		if (!is.null(ability)) {
			if (!is.null(rescale)) {
				return(list(link=link.out[[1]],pars=combine.pars(out.pars,x@common,grp.names),ability=out.ability))
			} else {
				return(list(link=link.out[[1]],ability=out.ability))
			}
		} else {
			if (!is.null(rescale)) {
				return(list(link=link.out[[1]],pars=combine.pars(out.pars,x@common,grp.names)))
			} else {
				return(link.out[[1]])
			}
		}
	} else {
		if (!is.null(ability)) {
			if (!is.null(rescale)) {
				return(list(link=link.out,pars=combine.pars(out.pars,x@common,grp.names),ability=out.ability))
			} else {
				return(list(link=link.out,ability=out.ability))
			}
		} else {
			if (!is.null(rescale)) {
				return(list(link=link.out,pars=combine.pars(out.pars,x@common,grp.names)))
			} else {
				return(link.out)
			}
		}
	}
})
