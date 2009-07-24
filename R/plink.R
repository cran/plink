setGeneric("plink", function(x, common, rescale, ability, method, weights.t, weights.f, startvals, score=1, base.grp=1, symmetric=FALSE, rescale.com=TRUE, grp.names=NULL, mn.exclude=0, dilation="ODL",  dim.order=NULL, ...) standardGeneric("plink"))



setMethod("plink", signature(x="list", common="list"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, score, base.grp, symmetric, rescale.com, grp.names, mn.exclude, dilation, dim.order, ...) {

	x <- combine.pars(x, common, ...)
	callGeneric()
	
})



setMethod("plink", signature(x="list", common="matrix"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, score, base.grp, symmetric, rescale.com, grp.names, mn.exclude, dilation, dim.order, ...) {

	x <- combine.pars(x, common, ...)
	callGeneric()
	
})



setMethod("plink", signature(x="list", common="data.frame"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, score, base.grp, symmetric, rescale.com, grp.names, mn.exclude, dilation, dim.order, ...) {
	
	x <- combine.pars(x, common, ...)
	callGeneric()
	
})



setMethod("plink", signature(x="irt.pars", common="ANY"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, score, base.grp, symmetric, rescale.com, grp.names, mn.exclude, dilation, dim.order, ...) {

	##   This is the function that will be minimized for the characteristic curve methods 
	##   (both unidimensional and multidimensional)
	.CC <- function(startvals, to, from, dimensions, weights.t, weights.f, score, transform, symmetric, mn.exclude, dilation, T, ...) {
		
		##   In general, the criterion that will minimized is Q = Q1+Q2
		##   where Q1 corresponds to the differences in probabilities after
		##   transforming the parameters on the FROM scale to the TO scale
		##   and Q2 corresponds to the differences in probabilities after
		##   transforming the parameters on the TO scale to the FROM scale.
		##   When {symmetric}=TRUE, Q is minimized, but when
		##   {symmetric}=FALSE, only Q1 will be minimized. 
		
		##   Identify optional arguments that might be passed
		##   to the various functions used to compute response probabilities
		dots <- list(...)
		if (length(dots$D)) D <- dots$D else D <- 1
		if (length(dots$D.drm)) D.drm <- dots$D.drm else D.drm <- D
		if (length(dots$D.gpcm)) D.gpcm <- dots$D.gpcm else D.gpcm <- D
		if (length(dots$D.grm)) D.grm <- dots$D.grm else D.grm <- D
		if (length(dots$incorrect)) incorrect <- dots$incorrect else incorrect <- FALSE
		if (length(dots$catprob)) catprob <- dots$catprob else catprob <- TRUE
		
		##   The theta values that will be used to compute response probabilities for Q1
		theta.t <- as.matrix(weights.t[[1]])
		
		##   The theta values that will be used to compute response probabilities for Q1
		theta.f <- as.matrix(weights.f[[1]])
		
		##   This object is used to store the TO scale parameters
		##   after they have been transformed to the FROM scale
		to.f <- to
		
		##   This object is used to store the FROM scale parameters
		##   after they have been transformed to the TO scale
		from.t <- from
		
		pm <- as.poly.mod(length(to@cat),to@model,to@items)
		
		if (dimensions==1) {
			##   Assign the vector of startvals to meaningful objects
			
			##   This parameterization is for the Rasch model
			##   where only the difficulties are transformed
			if (length(startvals)==1) {
				alpha <- 1
				beta <- startvals[1]
				
			##   This parameterization is for all other models
			} else {
				alpha <- startvals[1]
				beta <- startvals[2]
			}
			
			##   Identify the nrm and mcm items
			mcnr <- NULL
			for (i in 1:length(to@model)) {
				if (to@model[i]=="mcm"|to@model[i]=="nrm") mcnr <- c(mcnr, to@items[[i]])
			}
			
			##   Transform the FROM scale parameters onto the TO scale
			##   using the linking constants alpha and beta
			from.t@a <- as.matrix(from@a/alpha)
			from.t@b <- as.matrix(alpha * from@b + beta)
			from.t@b[mcnr,] <- from@b[mcnr,]-(beta/alpha)*from@a[mcnr,]
			
			##   Transform the FROM scale parameters onto the TO scale
			##   using the linking constants alpha and beta
			if (symmetric==TRUE) { 
				to.f@a <- as.matrix(alpha*to@a)
				to.f@b <- as.matrix((to@b-beta)/alpha)
				to.f@b[mcnr,] <- to@b[mcnr,]+(beta*to@a[mcnr,])
				
				##   Compute response probabilities using the TO scale parameters that 
				##   have been transformed onto the FROM scale
				prob.to.f <- .Mixed(to.f, theta.f, D.drm, D.gpcm, D.grm, incorrect, catprob)$p
				
				##   Compute response probabilities using the untransformed FROM scale parameters
				prob.from <- .Mixed(from, theta.f, D.drm, D.gpcm, D.grm, incorrect, catprob)$p
			}
		} else {
		
			if (length(startvals)==dimensions) {
				##   This parameterization is for the multidimensional Rasch case
				A <- diag(rep(1,dimensions))
				m <- startvals
				
			} else {
				##   This parameterization is used for all other models
				##   although it differs depending on the specified dilation
				if (dilation=="ODL") {
					A <- matrix(startvals[1:(dimensions^2)],dimensions,dimensions)
					m <- startvals[(dimensions^2+1):length(startvals)]
					
				##   For the LL and MIN approaches, the starting values need
				##   to be rotated using the orthogonal rotation matrix T
				##   prior to transforming any of the parameters. This makes it
				##   so that parameters only capture changes in variability
				} else if (dilation=="LL") {
					A <- diag(rep(startvals[1],dimensions))%*%T
					m <- startvals[2:length(startvals)]
					
				} else if (dilation=="MIN") {
					A <- diag(startvals[1:dimensions])%*%T
					m <- startvals[(dimensions+1):length(startvals)]
				} 
				
			}
			
			##   Transform the FROM scale parameters onto the TO scale
			##   using the linking constants A and m
			from.t@a <- from@a%*%ginv(A)
			from.t@b <- from@b-matrix(from.t@a %*% m, nrow(from@b),ncol(from@b))
			
			##   Transform the TO scale parameters onto the FROM scale
			##   using the linking constants A and m
			if (symmetric==TRUE) {
				to.f@a <- to@a%*%A
				to.f@b <- to@b+matrix(to@ a%*% m, nrow(to@b),ncol(to@b))
				
				##   Compute response probabilities using the TO scale parameters that 
				##   have been transformed onto the FROM scale
				prob.to.f <- .Mixed(to.f, theta.f, D.drm, D.gpcm, D.grm, incorrect, catprob)$p
				
				##   Compute response probabilities using the untransformed FROM scale parameters
				prob.from <- .Mixed(from, theta.f, D.drm, D.gpcm, D.grm, incorrect, catprob)$p
			}
		}
		
		##   Compute response probabilities using the untransformed TO scale parameters
		tmp.to <- .Mixed(to, theta.t, D.drm, D.gpcm, D.grm, incorrect, catprob)
		
		##   Use the information in the {irt.prob} object tmp.to to 
		##   determine the number of columns in the matrix of probabilities
		##   associated with each item (i.e., the object p.cat)
		prob.to <- tmp.to$p
		p.cat <- tmp.to$p.cat
		
		##   Compute response probabilities using the FROM scale parameters that 
		##   have been transformed onto the TO scale
		prob.from.t <- .Mixed(from.t, theta.t, D.drm, D.gpcm, D.grm, incorrect, catprob)$p
		
		##   Create a vector identifying the response model associated
		##   with each column of the matrix of probabilities
		p.mod <- rep(tmp.to$p.mod,p.cat)
		
		##   Initialize the vector for the scoring function
		##   that will be used for the Stocking-Lord method
		scr <- NULL
		
		##   Create a set of default values for the scoring function
		##   that range from 1 to Kj for K columns of probabilities for item j
		for (h in 1:length(p.cat)) {
			scr <- c(scr,seq(1,p.cat[h]))
		}
		
		##   Create the equivalent of a scoring function for the Haebara method
		##   This will be modified to exclude nrm and mcm items (if applicable)
		scr.hb <- rep(1,length(scr))
		
		##   Use this if one of the default values (1 or 2) is used for the {score} argument
		if (length(score)==1) {
		
			##   Use this if the lowest category should have a scoring weight of zero
			if (score==2) {
				scr <- scr-1
				
				##   If the argument {incorrect} equals FALSE, either explicitly or
				##   by not including the argument, there will only be one column of 
				##   probabilities for each dichotomous item. This formulation of the
				##   scoring function will set the weight for these columns equal to 
				##   zero. As such, set the weights for these items equal to one
				cat <- rep(p.cat,p.cat)
				scr[cat==1] <- 1
			}
			
			
		##   Use this if the researcher supplies a vector of score weights for
		##   all of the columns in p
		} else {
			if (length(score)==length(scr)) {
				scr <- score
			} else {
				warning("The length of {score} does not match the number of response probabilities. Score was set to 1")
			}
		}
		
		##   If nrm and mcm items are supposed to be excluded from the
		##   estimation of linking constants, set the probabilities for all responses
		##   for all of these items on both tests to zero
		if (mn.exclude %in% c(2,3)) {
			scr[p.mod=="nrm"|p.mod=="mcm"] <- 0
			scr.hb[p.mod=="nrm"|p.mod=="mcm"] <- 0
		}
		
		##   Compute sum of squares difference for the 
		##   Haebara and Stocking-Lord methods
		W1 <- as.vector(weights.t[[2]])
		W2 <- as.vector(weights.f[[2]])
		if (transform=="HB") {
			L1 <- ncol(prob.to)*sum(W1)
			Q1 <- W1*(prob.to-prob.from.t)^2
			Q1 <- sum(Q1 %*% scr.hb)
			if (symmetric==TRUE) {
				L2 <- ncol(prob.to)*sum(W2)
				Q2 <- W2*(prob.from-prob.to.f)^2
				Q2 <- sum(Q2 %*% scr.hb)
				SS <- (Q1/L1)+(Q2/L2)
			} else {
				SS <- Q1/L1
			}
		} else if (transform=="SL") {
			L1 <- sum(W1)
			TCC.to <- prob.to %*% scr
			TCC.from <- prob.from.t %*% scr
			Q1 <-sum(W1*(TCC.to - TCC.from)^2)
			if (symmetric==TRUE) {
				L2 <- sum(W2)
				TCC.to <- prob.to.f %*% scr
				TCC.from <- prob.from %*% scr
				Q2 <-sum(W2*(TCC.from- TCC.to)^2)
				SS <- (Q1/L1)+(Q2/L2)
			} else {
				SS <- Q1/L1
			}
		} 
		return(SS)
	}
	
	
	
	
	##   This function is used to compute either the  orthogonal or 
	##   oblique rotation matrix for the multidimensional approaches
	.Rotate <- function(to, from, orthogonal) {
		if (orthogonal==TRUE) {
			##   Orthogonal procrustes rotation
			M <- svd(t(to)%*%from)
			T <- M$v%*%t(M$u)
			
		} else {
			##   Oblique procrustes rotation
			T <- ginv(t(from) %*% from) %*% t(from) %*% to
		}
		
		TF <- T
		mn <- sum((to-from%*%ginv(T))^2)
		for (i in 1:3) {
			##   Due to indeterminacies in the sign of eigenvectors (or the orientation), 
			##   adjust the rotation matrix to get the matrix that truly minimizes the 
			##   difference between the two sets of slopes
			if (i==1) T1 <- t(T)
			if (i==2) T1 <- T*-1
			if (i==3) T1 <- t(T*-1)
			tmp <- sum((to-from%*%ginv(T1))^2)
			if (mn>tmp) {
				mn <- tmp
				TF <- T1
			}
		}
		return(T)
	}
	
	
	
	
	##   This is a stripped down version of the {mixed} function
	.Mixed <- function(x, theta, D.drm, D.gpcm, D.grm, incorrect, catprob) {
		
		mod <- x@model
		p <- NULL
		p.cat <- NULL
		p.mod <- NULL
		for (i in 1:length(mod)) {
			if (mod[i]=="drm") tmp <- .Drm(x, theta, D.drm, incorrect)
			if (mod[i]=="gpcm") tmp <- .Gpcm(x, theta, D.gpcm)
			if (mod[i]=="grm") tmp <- .Grm(x, theta, catprob, D.grm)
			if (mod[i]=="mcm") tmp <- .Mcm(x, theta)
			if (mod[i]=="nrm") tmp <- .Nrm(x, theta)
			p <- cbind(p, tmp$p)
			p.cat <- c(p.cat,tmp$cat)
			p.mod <- c(p.mod,tmp$mod)
		}
		return(list(p=p,p.cat=p.cat,p.mod=p.mod))
	}
	
	
	
	##   This is a stripped down version of the {drm} function
	.Drm <- function(x, theta, D.drm, incorrect) {
		
		dimensions <- x@dimensions
		items <- x@items$drm
		n <- length(items)
		a <- as.matrix(x@a[items,1:dimensions])
		if (n==1) a <- t(a)
		b <- x@b[items,1]
		c <- x@c[items,1]
		
		p <- NULL
		for (i in 1:length(b)) {
			if (dimensions==1) {
				cp <- c[i]+(1-c[i])/(1+exp(-D.drm*a[i]*(theta-b[i])))
			} else {
				cp <- c[i]+(1-c[i])/(1+exp(-(theta %*% a[i,]+b[i])))
				
			}
			if (incorrect==TRUE) p <- cbind(p,(1-cp),cp) else p <- cbind(p,cp)
		}
		if (incorrect==TRUE) cat <- rep(2,n) else cat <- rep(1,n)
		mod <- rep("drm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	
	
	##   This is a stripped down version of the {gpcm} function
	.Gpcm <- function(x, theta, D.gpcm) {
		
		dimensions <- x@dimensions
		items <- x@items$gpcm
		n <- length(items)
		a <- as.matrix(x@a[items,1:dimensions])
		b <- as.matrix(x@b[items,])
		if (length(items)==1) {
			a <- t(a)
			b <- t(b)
		}
		cat <- x@cat[items]
		
		p <- NULL 
		for (i in 1:nrow(b)) {
			dif <- 0 
			den <- NULL 
			
			for (k in 0:(cat[i]-1)) {
				if (k>=1) dif <- dif+b[i,k]
				if (dimensions==1) {
					d <- exp(D.gpcm*a[i]*(k*theta-dif))
				} else {
					d <- exp(k*(theta %*% a[i,])+dif)
				}
				den <- cbind(den, d)
			}
			den <- apply(den,1,sum)
			
			dif <- 0
			for (k in 0:(cat[i]-1)) {
				if (k>=1) dif <- dif+b[i,k]
				if (dimensions==1) {
					cp <- (exp(D.gpcm*a[i]*(k*theta-dif)))/den
				} else {
					cp <- exp(k*(theta %*% a[i,])+dif)/den
				}
				p <- cbind(p,cp)
			}
		}
		mod <- rep("gpcm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	
	
	##   This is a stripped down version of the {grm} function
	.Grm <- function(x, theta, catprob, D.grm) {
	
		dimensions <- x@dimensions
		items <- x@items$grm
		n <- length(items)
		a <- as.matrix(x@a[items,1:dimensions])
		b <- as.matrix(x@b[items,])
		if (length(items)==1) {
			a <- t(a)
			b <- t(b)
		}
		cat <- x@cat[items]
		p <- NULL 
		
		##   Compute category probabilities
		if (catprob==TRUE) { 
			for (i in 1:n) {
				ct <- cat[i]-1
				if (dimensions==1) {
					cp <- 1-1/(1+exp(-D.grm*a[i]*(theta-b[i,1])))
				} else {
					cp <- 1-1/(1+exp(-(theta %*% a[i,]+b[i,1])))
				}
				p <- cbind(p, cp)
				
				for (k in 1:ct) {
					if (dimensions==1) {
						if (k<ct) {
							cp <- (1/(1+exp(-D.grm*a[i]*(theta-b[i,k]))))-(1/(1+exp(-D.grm*a[i]*(theta-b[i,k+1]))))
						} else if (k==ct) {
							cp <- 1/(1+exp(-D.grm*a[i]*(theta-b[i,k])))
						}
					} else {
						if (k<ct) {
							cp <- (1/(1+exp(-(theta %*% a[i,]+b[i,k]))))-(1/(1+exp(-(theta %*% a[i,]+b[i,k+1]))))
						} else if (k==ct) {
							cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
						}
					}
					p <- cbind(p, cp)
				}
			}
			
		# Compute cumulative probabilities
		} else if (catprob==FALSE) { 
			for (i in 1:n) {
				for (k in 1:(cat[i]-1)) {
					if (dimensions==1) {
						cp <- 1/(1+exp(-D.grm*a[i]*(theta-b[i,k])))
					} else {
						cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
					}
					p <- cbind(p, cp)
				}
			}
		}
		if (catprob==FALSE) cat <- cat-1
		mod <- rep("grm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	
	
	##   This is a stripped down version of the {mcm} function
	.Mcm <- function(x, theta) {
	
		dimensions <- x@dimensions
		items <- x@items$mcm
		n <- length(items)
		a <- as.matrix(x@a[items,]) 
		b <- as.matrix(x@b[items,])
		c <- as.matrix(x@c[items,])
		if (length(items)==1) {
			a <- t(a)
			b <- t(b)
			c <- t(c)
		}
		cat <- x@cat[items]
		
		p <- NULL 
		for (i in 1:n) {
			den <- NULL
			a1 <- a[i,][!is.na(a[i,])]
			b1 <- b[i,][!is.na(b[i,])]
			c1 <- c[i,][!is.na(c[i,])]
			for (k in 1:cat[i]) {
				tmp <- (k-1)*dimensions
				tmp1 <- tmp+dimensions
				d <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])
				den <- cbind(den, d)
			}
			den <- apply(den,1,sum)
			for (k in 1:cat[i]) {
				tmp <- (k-1)*dimensions
				tmp1 <- tmp+dimensions
				if (k==1) {
					cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k]))/den
				} else {
					cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k])+c1[k-1]*(exp((theta %*% a1[1:dimensions])+b1[1])))/den
				}
				p <- cbind(p,cp)
			}
		}
		mod <- rep("mcm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	
	
	##   This is a stripped down version of the {nrm} function
	.Nrm <- function(x, theta) {
	
		dimensions <- x@dimensions
		items <- x@items$nrm
		n <- length(items)
		a <- as.matrix(x@a[items,]) 
		b <- as.matrix(x@b[items,])
		if (length(items)==1) {
			a <- t(a)
			b <- t(b)
		}
		cat <- x@cat[items]
		
		p <- NULL 
		for (i in 1:n) {
			den <- NULL
			a1 <- a[i,][!is.na(a[i,])]
			b1 <- b[i,][!is.na(b[i,])]
			for (k in 1:cat[i]) {
				tmp <- (k-1)*dimensions
				tmp1 <- tmp+dimensions
				d <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])
				den <- cbind(den, d)
			}
			den <- apply(den,1,sum)
			for (k in 1:cat[i]) {
				tmp <- (k-1)*dimensions
				tmp1 <- tmp+dimensions
				cp <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])/den
				p <- cbind(p,cp)
			}
		}
		mod <- rep("nrm",n)
		return(list(p=p,cat=cat,mod=mod))
	}
	
	
	
	##   The bult-in SD function R uses n-1 in the denominator
	##   This function uses n in the denominator
	.sd <- function(x) {
		z <- x[!is.na(x)]
		out <- sqrt(sum((z-mean(z))^2)/length(z))
		return(out)
	}
	
	
	
	##   This function computes descriptive statistics for the common item parameters
	.Descriptives <- function(a1, a2, b1, b2, c1, c2, pm, cat, mn.exclude, dimensions) {
	
	if (dimensions==1) {
		##   Transform the category parameters for mcm and nrm items in group 1
		if (ncol(a1)!=ncol(b1)) {
			b1r <- as.matrix(-b1/matrix(a1,nrow(a1),ncol(b1)))
		} else {
			b1r <- as.matrix(-b1/a1)
		}
		
		##   Transform the category parameters for mcm and nrm items in group 2
		if (ncol(a2)!=ncol(b2)) {
			b2r <- as.matrix(-b2/matrix(a2,nrow(a2),ncol(b2)))
		} else {
			b2r <- as.matrix(-b2/a2)
		}
		
		##   In the above transformations it is possible for some slope parameters
		##   to equal zero. As such the transformation will return values equal to
		##   infinity or negative infinity. These values should be set to NA.
		b1r[b1r== Inf] <- NA
		b2r[b2r== Inf] <- NA
		b1r[b1r== -Inf] <- NA
		b2r[b2r== -Inf] <- NA
		
		##   Initialize an object to store the descriptives
		descrip <- NULL
		
		##   Initialize an object to identify MCM or NRM items
		mcnr <- NULL
		
		##   Vector of response models for the common items
		pm.mod <- pm@model
		
		##   List of items associated with each model in pm.mod
		pm.it <- pm@items
		
		##   Extract the item parameters for each model in each group
		for (j in 1:length(pm.mod)) {
			a1k <- a1[pm.it[[j]],]
			a2k <- a2[pm.it[[j]],]
			b1k <- b1[pm.it[[j]],]
			b2k <- b2[pm.it[[j]],]
			c1k <- c1[pm.it[[j]],]
			c2k <- c2[pm.it[[j]],]
			b1kr <- b1r[pm.it[[j]],]
			b2kr <- b2r[pm.it[[j]],]
			
			##   Identify the number of each type of parameter, and compute the 
			##   means and SDs of each parameter type for both groups
			a <- c(length(a1k[!is.na(a1k)]),mean(a1k,na.rm=TRUE),mean(a2k,na.rm=TRUE),.sd(a1k),.sd(a2k))
			b <- c(length(b1k[!is.na(b1k)]),mean(b1k,na.rm=TRUE),mean(b2k,na.rm=TRUE),.sd(b1k),.sd(b2k))
			c <- c(length(c1k[!is.na(c1k)]),mean(c1k,na.rm=TRUE),mean(c2k,na.rm=TRUE),.sd(c1k),.sd(c2k))
			bd <- c(length(b1kr[!is.na(b1kr)]),mean(b1kr,na.rm=TRUE),mean(b2kr,na.rm=TRUE),.sd(b1kr),.sd(b2kr))
			
			##   Compile these descriptives for each model
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
		
		##   Check to see if all of the common items are modeled using the 
		##   NRM or MCM. If so, do not exclude any items from the computation
		##   of the descriptives (which will ultimately be used with the moment methods)
		if (sum((pm.mod=="nrm"|pm.mod=="mcm")+0)==length(pm.mod)) mn.exclude <- 0
		
		
		if (mn.exclude %in% c(0,2)) {
			b1[mcnr,] <- b1r[mcnr,] 
			b2[mcnr,] <- b2r[mcnr,] 
		} else {
			##   Remove NRM and/or MCM item parameters
			if (!is.null(mcnr)) {
				a1 <- a1[-mcnr,]
				a2 <- a2[-mcnr,]
				b1 <- b1[-mcnr,]
				b2 <- b2[-mcnr,]
			}
		}
		
		##   Identify the number of each type of parameter, and compute the 
		##   means and SDs of each parameter type for both groups
		a <- c(length(a1[!is.na(a1)]),mean(a1,na.rm=TRUE),mean(a2,na.rm=TRUE),.sd(a1),.sd(a2))
		b <- c(length(b1[!is.na(b1)]),mean(b1,na.rm=TRUE),mean(b2,na.rm=TRUE),.sd(b1),.sd(b2))
		
		##   Compile these descriptives for all "included" common items
		desall <- round(cbind(a,b),6)
		colnames(desall) <- c("a","b")
		rownames(desall) <- c("N Pars:","Mean: To","Mean: From","SD: To","SD: From")
		descrip[[length(descrip)+1]] <- desall
		names(descrip) <- c(pm.mod,"all")
		
	} else if (dimensions>1) {
		##   Initialize an object to store the descriptives
		descrip <- NULL
		
		##   Initialize an object to identify MCM or NRM items
		mcnr <- NULL
		
		##   Vector of response models for the common items
		pm.mod <- pm@model
		
		##   List of items associated with each model in pm.mod
		pm.it <- pm@items
		
		a1.all <- a2.all <- vector("list",dimensions)
		
		##   Extract the item parameters for each model in each group
		for (j in 1:length(pm.mod)) {
			a1k <- a1[pm.it[[j]],]
			a2k <- a2[pm.it[[j]],]
			b1k <- b1[pm.it[[j]],]
			b2k <- b2[pm.it[[j]],]
			c1k <- c1[pm.it[[j]],]
			c2k <- c2[pm.it[[j]],]
			catk <- cat[pm.it[[j]]]
			
			a <- NULL
			if (pm.mod[j]=="nrm"|pm.mod[j]=="mcm") {
				mcatk <- max(catk)
				for (k in 1:dimensions) {
					tmp1 <- a1k[,(((k-1)*mcatk)+1):(k*mcatk)]
					tmp2 <- a2k[,(((k-1)*mcatk)+1):(k*mcatk)]
					a <- cbind(a,c(length(tmp1[!is.na(tmp1)]),mean(tmp1,na.rm=TRUE),mean(tmp2,na.rm=TRUE),.sd(tmp1),.sd(tmp2)))
					a1.all[[k]] <- c(a1.all[[k]],as.vector(tmp1))
					a2.all[[k]] <- c(a2.all[[k]],as.vector(tmp2))
				}
			} else {
				for (k in 1:dimensions) {
					tmp1 <- a1k[,k]
					tmp2 <- a2k[,k]
					a <- cbind(a,c(length(tmp1[!is.na(tmp1)]),mean(tmp1,na.rm=TRUE),mean(tmp2,na.rm=TRUE),.sd(tmp1),.sd(tmp2)))
					a1.all[[k]] <- c(a1.all[[k]],as.vector(tmp1))
					a2.all[[k]] <- c(a2.all[[k]],as.vector(tmp2))
				}
			}
			
			##   Identify the number of each type of parameter, and compute the 
			##   means and SDs of each parameter type for both groups
			b <- c(length(b1k[!is.na(b1k)]),mean(b1k,na.rm=TRUE),mean(b2k,na.rm=TRUE),.sd(b1k),.sd(b2k))
			c <- c(length(c1k[!is.na(c1k)]),mean(c1k,na.rm=TRUE),mean(c2k,na.rm=TRUE),.sd(c1k),.sd(c2k))
			
			##   Compute MDISC and MDIF for all common items for both groups
			mdisc1 <- sqrt(apply(a1k^2,1,sum,na.rm=TRUE))
			mdisc2 <- sqrt(apply(a2k^2,1,sum,na.rm=TRUE))
			mdif1 <- -b1k/mdisc1
			mdif2 <- -b2k/mdisc2
			
			##   Identify the number MDISC and MDIF parameters, and compute the 
			##   means and SDs of each for both groups
			mdc <- c(length(mdisc1[!is.na(mdisc1)]),mean(mdisc1,na.rm=TRUE),mean(mdisc2,na.rm=TRUE),.sd(mdisc1),.sd(mdisc2))
			mdf <- c(length(mdif1[!is.na(mdif1)]),mean(mdif1,na.rm=TRUE),mean(mdif2,na.rm=TRUE),.sd(mdif1),.sd(mdif2))
			
			##   Compile these descriptives for each model
			if (pm.mod[j]=="drm" ) {
				des <- data.frame(cbind(a,b,c,mdc,mdf))
				colnames(des) <- c(paste("a",1:dimensions,sep=""),"d","c","MDISC","MDIF")
			} else if (pm.mod[j]=="mcm") {
				des <- data.frame(cbind(a,b,c,mdc,mdf))
				colnames(des) <- c(paste("a",1:dimensions,sep=""),"d","c","MDISC","MDIF")
			} else {
				des <- data.frame(cbind(a,b,mdc,mdf))
				colnames(des) <- c(paste("a",1:dimensions,sep=""),"d","MDISC","MDIF")
			}
			
			if (dim.flag==TRUE) des <- des[,-2]  ##  WHAT IS THIS?
			colnames(des)[1] <- "a"
			rownames(des) <- c("N Pars:","Mean: To","Mean: From","SD: To","SD: From")
			descrip[[j]] <- round(des,6)
		}
		
		a <- NULL
		for (k in 1:dimensions) {
			tmp1 <- a1.all[[k]]
			tmp2 <- a2.all[[k]]
			a <- cbind(a,c(length(tmp1[!is.na(tmp1)]),mean(tmp1,na.rm=TRUE),mean(tmp2,na.rm=TRUE),.sd(tmp1),.sd(tmp2)))
		}
		b <- c(length(b1[!is.na(b1)]),mean(b1,na.rm=TRUE),mean(b2,na.rm=TRUE),.sd(b1),.sd(b2))
		mdisc1 <- sqrt(apply(a1^2,1,sum,na.rm=TRUE))
		mdisc2 <- sqrt(apply(a2^2,1,sum,na.rm=TRUE))
		mdif1 <- -b1/mdisc1
		mdif2 <- -b2/mdisc2
		mdc <- c(length(mdisc1[!is.na(mdisc1)]),mean(mdisc1,na.rm=TRUE),mean(mdisc2,na.rm=TRUE),.sd(mdisc1),.sd(mdisc2))
		mdf <- c(length(mdif1[!is.na(mdif1)]),mean(mdif1,na.rm=TRUE),mean(mdif2,na.rm=TRUE),.sd(mdif1),.sd(mdif2))
		
		desall <- round(cbind(a,b,mdc,mdf),6)
		colnames(desall) <- c(paste("a",1:dimensions,sep=""),"d","MDISC","MDIF")
		rownames(desall) <- c("N Pars:","Mean: To","Mean: From","SD: To","SD: From")
		descrip[[length(descrip)+1]] <- desall
		names(descrip) <- c(pm.mod,"all")
		}
	
	return(descrip)
	}
	
	
	
	
	### START PLINK ###
	ng <- x@groups
	md <- max(x@dimensions)
	 if (missing(grp.names)) grp.names <- names(x@pars)
	
	# Create or recode the values in dim.order
	# dim.order.RM is the originally specified dim.order
	# The columns correspond to the total number of dimensions across groups
	# The rows correspond to the groups
	if (md>1) {
		if (is.null(dim.order)) {
			dim.order.RM <- NULL
			dim.order <- matrix(NA,ng,md)
			for (i in 1:ng) {
				# Set the default ordering to be the same as the ordering of supplied slope parameters
				dim.order[i,1:x@dimensions[i]] <- 1:x@dimensions[i]
			}
		} else {
			# If ones are used as placeholders for common dimensions
			if (sum(dim.order, na.rm=TRUE)==sum(x@dimensions)) {
				dim.order.RM <- dim.order
				for (i in 1:ng) {
					dim.order[i,!is.na(dim.order[i,])] <- 1:x@dimensions[i]
				}
			} else {
				dim.order.RM <- dim.order
			}
		}
	}
	
	# Extract common items
	com <- vector("list",ng-1)
	catci <- vector("list",ng-1)
	pm <- vector("list",ng-1)
	for (i in 1:(ng-1)) {
		com[[i]] <- vector("list",2)
		if (ng==2) {
			it.com <- x@common
			com[[i]][[1]] <- x@pars[[i]][it.com[,1],] # Common item parameterss for the lower group
			com[[i]][[2]] <- x@pars[[i+1]][it.com[,2],] # Common item parameters for the higher group
			catci[[i]] <- x@cat[[i]][it.com[,1]] # Category vector for the common items
		} else if (ng>2) {
			it.com <- x@common[[i]]
			com[[i]][[1]] <- x@pars[[i]][it.com[,1],] # Common item parameters for the lower group
			com[[i]][[2]] <- x@pars[[i+1]][it.com[,2],] # Common item parameters for the higher group
			catci[[i]] <- x@cat[[i]][it.com[,1]] # Category vector for the common items
		}
		if (i<base.grp) com[[i]] <- com[[i]][c(2,1)] # The first element is the "To" set and the second is the "From" set
		mod <-  x@poly.mod[[i]]@model
		items <- x@poly.mod[[i]]@items
		#Extract the common item numbers for the "To" group to facilitate the creation of a poly.mod object for the common items
		if (ng==2) it.com <- x@common[,1] else it.com <- x@common[[i]][,1] 
		poly <- mod1 <- NULL
	
		# Identify the common items associated with each model
		# Not all models need to have common items
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
	
	# Extract the common item parameters 
	link.out <- vector("list",ng-1)
	for (i in 1:(ng-1)) {
		tmp1 <- sep.pars(com[[i]][[1]],catci[[i]],pm[[i]],x@dimensions[i],x@location[i]) # "To" parameters
		tmp2 <- sep.pars(com[[i]][[2]],catci[[i]],pm[[i]],x@dimensions[i+1],x@location[i+1]) # "From" parameters
		tmp.md <- max(c(x@dimensions[i],x@dimensions[i+1]))
		
		# Identify common dimensions
		dim.flag <- FALSE
		if (tmp.md>1) {
			do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
			do2 <- dim.order[i+1,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
			tmp1@a <- as.matrix(as.matrix(tmp1@a)[,do1])
			tmp2@a <- as.matrix(as.matrix(tmp2@a)[,do2])
			tmp1@dimensions <- tmp2@dimensions <- dimensions <- length(do1)
			for (j in i:(i+1)) {
				# If the parameters for a given group are unidimensional
				if (x@dimensions[j]==1) {
					for (k in 1:length(pm[[i]]@model)) {
						if (pm[[i]]@model[k]=="drm"|pm[[i]]@model[k]=="grm"|pm[[i]]@model[k]=="gpcm") {
							tmp1@b[pm[[i]]@items[[k]],] <- as.matrix(tmp1@b[pm[[i]]@items[[k]],]*-as.vector(tmp1@a)[pm[[i]]@items[[k]]])
							tmp2@b[pm[[i]]@items[[k]],] <- as.matrix(tmp2@b[pm[[i]]@items[[k]],]*-as.vector(tmp2@a)[pm[[i]]@items[[k]]])
						}
					}
				}
			}
							
			# If there is only one common dimension
			if (dimensions==1) {
				tmp1@a <- cbind(tmp1@a,0)
				tmp2@a <- cbind(tmp2@a,0)
				tmp1@dimensions <- tmp2@dimensions <- dimensions <- 2
				dim.flag <- TRUE
			}
		} else {
			dimensions <- 1
		}
		
		a1 <- tmp1@a
		a2 <- tmp2@a
		b1 <- tmp1@b
		b2 <- tmp2@b
		c1 <- tmp1@c
		c2 <- tmp2@c
		
		if (missing(weights.t)) wgt.t <- as.weight(dimensions=dimensions) else wgt.t <- weights.t 
		if (missing(weights.f)) wgt.f <- wgt.t else wgt.f <- weights.f 
		
		#if (is.list(weights[[1]])) wgt <- weights[[i]] else 
		
		if (missing(method)) {
			if (dimensions==1) method <- c("MM","MS","HB","SL") else method <- c("RM","HB","SL") 
		} else {
			method <- toupper(method)
			if (sum(method%in% c("MM","MS","HB","SL","RM") )==0) {
				warning("No appropriate method was selected. All methods will be used")
				if (dimensions==1) method <- c("MM","MS","HB","SL") else method <- c("RM","HB","SL") 
			}
		}
		
		descrip <- .Descriptives(a1,a2,b1,b2,c1,c2,pm[[i]],catci[[i]],mn.exclude,dimensions)
		
		chk.sl <- TRUE
		if (sum((mod=="nrm"|mod=="mcm")+0)==length(mod)) {
			chk.sl <- FALSE
			meth <- NULL
			if ("MS" %in% method) meth <- "Mean/Sigma"
			if ("MM" %in% method) {
				if (length(meth)) meth <- "Mean/Sigma and Mean/Mean" else meth <- "Mean/Mean"
			}
			if (length(meth)) cat(paste("All items were used to compute the",meth,"linking constants\n"))
		}
		
		# Perform the Calibration
		rasch.flag <- FALSE
		constants <- list(NULL)
		it <- con <- obj <- NULL
		
		if (dimensions==1) {
			# Compute the mean/mean and mean/sigma constants
			mm <- ms <- NULL
			A1 <- descrip$all[3,1]/descrip$all[2,1]
			A2 <- descrip$all[4,2]/descrip$all[5,2]
			#A1 <- mean(a2,na.rm=TRUE)/mean(a1,na.rm=TRUE)
			#A2 <- .sd(b1)/.sd(b2)
			if (rasch.flag==TRUE) A2 <- a1[1]
			B1 <- descrip$all[2,2]-A1*descrip$all[3,2]
			B2 <- descrip$all[2,2]-A2*descrip$all[3,2]
			#B1 <- mean(b1,na.rm=TRUE)-A1*mean(b2,na.rm=TRUE)
			#B2 <- mean(b1,na.rm=TRUE)-A2*mean(b2,na.rm=TRUE)
			mm <- round(c(A1,B1),6)
			ms <- round(c(A2,B2),6)
			names(mm) <- names(ms) <- c("A","B")
			if ("MM" %in% method) constants$MM <- mm
			if ("MS" %in% method) constants$MS <- ms
			dilation <- "N/A"
			T <- NA
			
			# Determine starting values for the characteristic curve methods
			if (missing(startvals)) {
				startvals <- ms
			} else {
				if (is.character(startvals)) {
					if (toupper(startvals)=="MM") startvals <- mm
					if (toupper(startvals)=="MS") startvals <- ms
				}
			}
			tmp <- c(a1,a2)
			if (length(tmp[tmp==1])==length(tmp)) {
				startvals <- startvals[2]
				rasch.flag <- TRUE
			}
			lower <- startvals-0.25
			upper <- startvals+0.25
		} else {
			# Direct method
			tmp <- as.vector(a1)
			if (length(tmp[tmp==1])==length(tmp)) {
				A <- diag(rep(1,dimensions))
				rasch.flag <- TRUE
			} else {
				# This actually estimates the inverse of A
				A <- .Rotate(a1,a2,FALSE)
			}
			
			tmp.b1 <- as.vector(b1)[!is.na(as.vector(b1))]
			tmp.b2 <- as.vector(b2)[!is.na(as.vector(b2))]
			tmp.a2 <- NULL
			for (j in 1:dimensions) {
				tmp.a2a <- (!is.na(b2))*a2[,j]
				tmp.a2a[is.na(b2)] <- NA
				tmp.a2a <- as.vector(tmp.a2a)[!is.na(as.vector(tmp.a2a))]
				tmp.a2 <- cbind(tmp.a2,tmp.a2a)
			}
			y <- tmp.b2-tmp.b1 
			X <- tmp.a2%*%A
			m <- ginv(t(X)%*%X)%*%t(X)%*%y
			m <- as.vector(m)
			A <- ginv(A)
			rownames(A) <- rep("",dimensions)
			colnames(A) <- c(rep("",dimensions-1),"A")
			names(m) <- paste("m",1:dimensions,sep="")
			if ("RM" %in% method) {
				if (dim.flag==FALSE) {
					constants$RM <- list(A=round(A,6),m=round(m,6))
				} else {
					constants$RM <- c(round(A[1,1],6),round(m[1],6))
					names(constants$RM) <- c("A","m")
				}
			}
			
			# Compute the rotation matrix for certain characteristic curve methods
			if (rasch.flag==TRUE) T <- diag(rep(1,dimensions)) else T <- .Rotate(a1,a2,TRUE)
			
			# Determine starting values for the characteristic curve methods
			if (rasch.flag==TRUE) {
				if (missing(startvals)) startvals <- m else startvals <- startvals[((length(startvals)-dimensions)+1):length(startvals)]
			} else {
				if (missing(startvals)) {
					if (dilation=="ODL") {
						startvals <- c(as.vector(A),m) 
						lower <- startvals-0.25
						upper <- startvals+0.25
					} else if (dilation=="LL") {
						startvals <- c(1,m)
						lower <- c(-2,m-0.25)
						upper <- c(2,m+0.25)
					} else if (dilation=="MIN") {
						startvals <- c(rep(1,dimensions),m)
						lower <- c(rep(-2,dimensions),m-0.25)
						upper <- c(rep(2,dimensions),m+0.25)
					} 
				} else {
					if (dilation!="ODL") {
						if ((length(startvals)-dimensions)!=dimensions^2) stop(paste("The number of elements must equal",dimensions^2))
					}
					lower <- startvals-0.5
					upper <- startvals+0.5
				}
			}
			
		}
		
		# Characteristic curve methods
		
		# Haebara Method
		if ("HB" %in% method) {
			hb <- nlminb(startvals, .CC, to=tmp1, from=tmp2, dimensions=dimensions, weights.t=wgt.t, weights.f=wgt.f, score=score, transform="HB", symmetric=symmetric, mn.exclude=mn.exclude, dilation=dilation, T=T, lower=lower, upper=upper, ...)
			if (dimensions==1) { 
				if (rasch.flag==TRUE) hb$par <- c(1,hb$par)
				names(hb$par) <- c("A","B")
				constants$HB <- round(hb$par,6)
			} else {
				if (dilation=="ODL") {
					if (rasch.flag==TRUE) A <- diag(rep(1,dimensions)) else A <- matrix(hb$par[1:(dimensions^2)],dimensions,dimensions)
					m <- hb$par[-c(1:(dimensions^2))]
				} else if (dilation=="LL") {
					if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(rep(hb$par[1],dimensions))
					m <- hb$par[-1]
					A <- T%*%K
				} else if (dilation=="MIN") {
					if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(hb$par[1:dimensions])
					m <- hb$par[-c(1:dimensions)]
					A <- T%*%K
				} 
				rownames(A) <-rep("",dimensions)
				colnames(A) <- c(rep("",dimensions-1),"A")
				names(m) <- paste("m",1:dimensions,sep="")
				if (dilation=="ODL") {
					constants$HB <- list(A=round(A,6), m=round(m,6))
				} else {
					rownames(T) <- rownames(K) <- rep("",dimensions)
					colnames(T) <- c(rep("",dimensions-1),"T")
					colnames(K) <- c(rep("",dimensions-1),"K")
					constants$HB <- list(T=round(T,6), K=round(K,6), A=round(A,6), m=round(m,6))
				}
			}
			if (dim.flag==TRUE) {
				constants$HB <- c(constants$HB$A[1,1],constants$HB$m[1])
				names(constants$HB) <- c("A","m")
			}
			it <- c(it, HB=hb$iterations)
			con <- c(con, HB=hb$message)
			obj <- c(obj, HB=hb$objective)
		}
		
		
		# Stocking-Lord Method
		if ("SL" %in% method) {
			if (chk.sl==TRUE) {
				sl <- nlminb(startvals, .CC, to=tmp1, from=tmp2, dimensions=dimensions, weights.t=wgt.t, weights.f=wgt.f, score=score, transform="SL", symmetric=symmetric, mn.exclude=mn.exclude, dilation=dilation, T=T, lower=lower, upper=upper, ...)
				if (dimensions==1) { 
					if (rasch.flag==TRUE) sl$par <- c(1,sl$par)
					names(sl$par) <- c("A","B")
					constants$SL <- round(sl$par,6)
				} else {
					if (dilation=="ODL") {
						if (rasch.flag==TRUE) A <- diag(rep(1,dimensions)) else A <- matrix(sl$par[1:(dimensions^2)],dimensions,dimensions)
						m <- sl$par[-c(1:(dimensions^2))]
					} else if (dilation=="LL") {
						if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(rep(sl$par[1],dimensions))
						m <- sl$par[-1]
						A <- T%*%K
					} else if (dilation=="MIN") {
						if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(sl$par[1:dimensions])
						m <- sl$par[-c(1:dimensions)]
						A <- T%*%K
					} 
					rownames(A) <-rep("",dimensions)
					colnames(A) <- c(rep("",dimensions-1),"A")
					names(m) <- paste("m",1:dimensions,sep="")
					if (dilation=="ODL") {
						constants$SL <- list(A=round(A,6), m=round(m,6))
					} else {
						rownames(T) <- rownames(K) <- rep("",dimensions)
						colnames(T) <- c(rep("",dimensions-1),"T")
						colnames(K) <- c(rep("",dimensions-1),"K")
						constants$SL <- list(T=round(T,6), K=round(K,6), A=round(A,6), m=round(m,6))
					}
				}
				if (dim.flag==TRUE) {
					constants$SL <- c(constants$SL$A[1,1],constants$SL$m[1])
					names(constants$SL) <- c("A","m")
				}
				it <- c(it, SL=sl$iterations)
				con <- c(con, SL=sl$message)
				obj <- c(obj, SL=sl$objective)
			}
		}
		
		constants[[1]] <- NULL
		if (is.null(it)) it <- 0
		if (is.null(con)) con <- "N/A"
		if (is.null(obj)) obj <- 0
		
		if (i==base.grp) {
			nms <- paste(grp.names[i+1],"/",grp.names[i],"*",sep="")
		} else if (i < base.grp) {
			if ((i+1)==base.grp) {
				nms <- paste(grp.names[i],"/",grp.names[i+1],"*",sep="")
			} else {
				nms <- paste(grp.names[i],"/",grp.names[i+1],sep="")
			}
		} else if (i > base.grp) {
			nms <- paste(grp.names[i+1],"/",grp.names[i],sep="")
		}
		
		link.out[[i]] <- new("link", constants=constants, descriptives=descrip, iterations=it, objective=obj, convergence=con, base.grp=base.grp, grp.names=nms, include.mcm.nrm=mn.exclude,n=tmp1@n, mod.lab=tmp1@mod.lab, dilation=dilation)
		
	}
	
	if (!missing(rescale)) {
		if ((toupper(rescale)%in%method)==FALSE) {
			if (dimensions==1) {
				if ("SL" %in% method) rescale <- "SL" else rescale <- method[length(method)]
			} else {
				if ("RM" %in% method) rescale <- "RM" else rescale <- method[length(method)]
			}
			warning(paste("No linking constants were computed for the rescale method you selected. The parameters will be rescaled using the ",rescale," method",sep=""))
		}
		tmp.con <- vector("list",ng)
		j <- 1
		for (i in 1:ng) {
			if (i==base.grp) {
				if (x@dimensions[i]==1) tmp.con[[i]] <- c(1,0) else tmp.con[[i]] <- list(A=diag(rep(1,x@dimensions[i])),m=rep(0,x@dimensions[i]))
			} else {
				tmp.con[[i]] <- eval(parse(text=paste("link.out[[",j,"]]@constants$",rescale,sep="")))
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
			tmp <- sep.pars(x@pars[[i]],x@cat[[i]],x@poly.mod[[i]],x@dimensions[i],x@location[i])
			if (!missing(ability)) tmpa <- ability[[i]] 
			dimensions <- x@dimensions[i]
			tmp1 <- tmp
			if (i==base.grp) {
				tmp.j <- i
			} else if (i < base.grp) {
				tmp.j <- i:(base.grp-1)
			} else if (i > base.grp) {
				tmp.j <- i:(base.grp+1)
			}
			# Perform j iterations to place the parameters for the given group on the base group scale
			for (j in tmp.j) {
				if (dimensions==1) {
					tmp1@a <- as.matrix(tmp1@a/tmp.con[[j]][1])
					tmp1@b <-  as.matrix(tmp.con[[j]][1]*tmp1@b + tmp.con[[j]][2])
					tmp1@b[mcnr,] <- tmp@b[mcnr,]-(tmp.con[[j]][2]/tmp.con[[j]][1])*tmp@a[mcnr,]
					tmp@b[mcnr,] <- tmp1@b[mcnr,]
					if (!missing(ability)) tmpa <- tmp.con[[j]][1]*tmpa + tmp.con[[j]][2]
				} else {
					if (i==base.grp) {
						do1 <- dim.order[j,!is.na(dim.order[j,])]
					} else if (i < base.grp) {
						do1 <- dim.order[j,!is.na(dim.order[j,]) &!is.na(dim.order[j+1,])]
					} else if (i > base.grp) {
						do1 <- dim.order[j,!is.na(dim.order[j,]) &!is.na(dim.order[j-1,])]
					}
					
					if (length(do1)==1) {
						tmp1@a[,do1] <- tmp1@a[,do1]*ginv(tmp.con[[j]][1])
						tmp1@b <- tmp1@b-matrix(tmp1@a[,do1]*tmp.con[[j]][2],nrow(tmp1@b),ncol(tmp1@b))
						if (!missing(ability)) tmpa <- tmp.con[[j]][1]*tmpa + tmp.con[[j]][2]
					} else {
						tmp1@a[,do1] <- tmp1@a[,do1]%*%ginv(tmp.con[[j]]$A)
						tmp1@b <- tmp1@b-matrix(tmp1@a[,do1]%*%tmp.con[[j]]$m,nrow(tmp1@b),ncol(tmp1@b))
						if (!missing(ability)) tmpa <- t(tmp.con[[j]]$A%*%t(tmpa[,do1]) + tmp.con[[j]]$m)
					}
				}
			}
			out.pars[[i]] <- tmp1
			names(out.pars)[[i]] <- grp.names[i]
			if (!missing(ability)) {
				out.ability[[i]] <- tmpa
				names(out.ability)[[i]] <- grp.names[i]
			}
		}
		# Substitute the non-scaled parameters for the common items into the set of rescaled item parameters
		if (rescale.com==FALSE) {
			com <- x@common
			if (is.matrix(com)) com <- list(com)
			for (i in (base.grp-1):1) {
				if (base.grp==1) break
				b.min <- min(c(ncol(out.pars[[i]]@b),ncol(out.pars[[i+1]]@b)))
				c.min <- min(c(ncol(out.pars[[i]]@c),ncol(out.pars[[i+1]]@c)))
				
				out.pars[[i]]@b[com[[i]][,1],1:b.min] <- out.pars[[i+1]]@b[com[[i]][,2],1:b.min]
				out.pars[[i]]@c[com[[i]][,1],1:c.min] <- out.pars[[i+1]]@c[com[[i]][,2],1:c.min]
				if (dimensions==1) {
					out.pars[[i]]@a[com[[i]][,1],] <- out.pars[[i+1]]@a[com[[i]][,2],]
				} else {
					do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					do2 <- dim.order[i+1,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					out.pars[[i]]@a[com[[i]][,1],do1] <- out.pars[[i+1]]@a[com[[i]][,2],do2]
					
				}
			}
			
			for (i in base.grp:(ng-1)) {
				if (base.grp==ng) break
				b.min <- min(c(ncol(out.pars[[i]]@b),ncol(out.pars[[i+1]]@b)))
				c.min <- min(c(ncol(out.pars[[i]]@c),ncol(out.pars[[i+1]]@c)))
				
				out.pars[[i+1]]@b[com[[i]][,2],1:b.min] <- out.pars[[i]]@b[com[[i]][,1],1:b.min]
				out.pars[[i+1]]@c[com[[i]][,2],1:c.min] <- out.pars[[i]]@c[com[[i]][,1],1:c.min]
				if (dimensions==1) {
					out.pars[[i+1]]@a[com[[i]][,2],] <- out.pars[[i]]@a[com[[i]][,1],]
				} else {
					do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					do2 <- dim.order[i+1,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					out.pars[[i+1]]@a[com[[i]][,2],do2] <- out.pars[[i]]@a[com[[i]][,1],do1]
				}
			}
			
		}
		
	} else {
		if (!missing(ability)) {
			if (dimensions==1) {
				if ("SL" %in% method) rsc <- "SL" else rsc<- method[length(method)]
			} else {
				if ("RM" %in% method) rsc <- "RM" else rsc<- method[length(method)]
			}
			cat(paste("No rescale method was identified. Ability parameters will be rescaled using the ",rsc," method.\n",sep=""))
			tmp.con <- vector("list",ng)
			j <- 1
			for (i in 1:ng) {
				if (i==base.grp) {
					if (dimensions==1) tmp.con[[i]] <- c(1,0) else tmp.con[[i]] <- list(A=diag(rep(1,dimensions)),m=rep(0,dimensions))
				} else {
					tmp.con[[i]] <- eval(parse(text=paste("link.out[[",j,"]]@constants$",rsc,sep="")))
					j <- j+1
				}
			}
			
			out.ability <- vector("list",ng)
			for (i in 1:ng) {
				tmpa <- ability[[i]]
				if (i==base.grp) {
					tmp.j <- i
				} else if (i < base.grp) {
					tmp.j <- i:(base.grp-1)
				} else if (i > base.grp) {
					tmp.j <- i:(base.grp+1)
				}
				# Perform j iterations to place the ability estimates for the given group on the base group scale
				for (j in tmp.j) {
					
					if (dimensions==1) {
						tmpa <- tmp.con[[j]][1]*tmpa + tmp.con[[j]][2]
					} else {
						if (i==base.grp) {
							do1 <- dim.order[j,!is.na(dim.order[j,])]
						} else if (i < base.grp) {
							do1 <- dim.order[j,!is.na(dim.order[j,]) &!is.na(dim.order[j+1,])]
						} else if (i > base.grp) {
							do1 <- dim.order[j,!is.na(dim.order[j,]) &!is.na(dim.order[j-1,])]
						}
						
						if (length(do1)==1) {
							tmpa <- tmp.con[[j]][1]*tmpa + tmp.con[[j]][2]
						} else {
							tmpa <- t(tmp.con[[j]]$A%*%t(tmpa[,do1]) + tmp.con[[j]]$m)
						}
					}
				}
				out.ability[[i]] <- tmpa
				names(out.ability)[[i]] <- grp.names[i]
			}
		}
	}
	
	if (ng==2) {
		if (!missing(ability)) {
			if (!missing(rescale)) {
				return(list(link=link.out[[1]],pars=combine.pars(out.pars,x@common,grp.names),ability=out.ability))
			} else {
				return(list(link=link.out[[1]],ability=out.ability))
			}
		} else {
			if (!missing(rescale)) {
				return(list(link=link.out[[1]],pars=combine.pars(out.pars,x@common,grp.names)))
			} else {
				return(link.out[[1]])
			}
		}
	} else {
		if (!missing(ability)) {
			if (!missing(rescale)) {
				return(list(link=link.out,pars=combine.pars(out.pars,x@common,grp.names),ability=out.ability))
			} else {
				return(list(link=link.out,ability=out.ability))
			}
		} else {
			if (!missing(rescale)) {
				return(list(link=link.out,pars=combine.pars(out.pars,x@common,grp.names)))
			} else {
				return(link.out)
			}
		}
	}
})