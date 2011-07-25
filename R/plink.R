##   This function estimates the linking constants between two or
##   more groups of tests and can rescale item and/or ability parameters

setGeneric("plink", function(x, common, rescale, ability, method, weights.t, weights.f, startvals, exclude, score=1, base.grp=1, symmetric=FALSE, rescale.com=TRUE, grp.names=NULL, dilation="LL",  dim.order=NULL, ...) standardGeneric("plink"))



setMethod("plink", signature(x="list", common="list"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, exclude, score, base.grp, symmetric, rescale.com,  grp.names, dilation, dim.order, ...) {

	x <- combine.pars(x, common, ...)
	callGeneric()
	
})



setMethod("plink", signature(x="list", common="matrix"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, exclude, score, base.grp, symmetric, rescale.com, grp.names, dilation, dim.order, ...) {

	x <- combine.pars(x, common, ...)
	callGeneric()
	
})



setMethod("plink", signature(x="list", common="data.frame"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, exclude, score, base.grp, symmetric, rescale.com, grp.names, dilation, dim.order, ...) {
	
	x <- combine.pars(x, common, ...)
	callGeneric()
	
})



setMethod("plink", signature(x="irt.pars", common="ANY"), function(x, common, rescale, ability, method, weights.t, weights.f, startvals, exclude, score, base.grp, symmetric, rescale.com, grp.names, dilation, dim.order, ...) {


	##################################################################
	##      Function that will be minimized for the characteristic curve methods
	##                       (both unidimensional and multidimensional)
	##################################################################
	
	.CC <- function(sv, to, from, dimensions, weights.t, weights.f, sc, transform, symmetric, dilation, T, ...) {
		##   In general, the criterion that will minimized is Q = Q1+Q2
		##   where Q1 corresponds to the differences in probabilities after
		##   transforming the parameters on the FROM scale to the TO scale
		##   and Q2 corresponds to the differences in probabilities after
		##   transforming the parameters on the TO scale to the FROM scale.
		##   When {symmetric}=TRUE, Q is minimized, but when
		##   {symmetric}=FALSE, only Q1 will be minimized. 
		
		##  Identify optional arguments that might be passed
		##  to the various functions used to compute response probabilities
		dots <- list(...)
		if (length(dots$D)) D <- dots$D else D <- 1
		if (length(dots$D.drm)) D.drm <- dots$D.drm else D.drm <- D
		if (length(dots$D.gpcm)) D.gpcm <- dots$D.gpcm else D.gpcm <- D
		if (length(dots$D.grm)) D.grm <- dots$D.grm else D.grm <- D
		if (length(dots$catprob)) catprob <- dots$catprob else catprob <- TRUE
		
		if (transform=="SL") incorrect <- TRUE else incorrect <- FALSE
		
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
			##   Assign the vector of startvals (sv) to meaningful objects
			
			##   This parameterization is for the Rasch model
			##   where only the difficulties are transformed
			if (length(sv)==1) {
				alpha <- 1
				beta <- sv[1]
				
			##   This parameterization is for all other models
			} else {
				alpha <- sv[1]
				beta <- sv[2]
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
		
			if (length(sv)==dimensions) {
				##   This parameterization is for the multidimensional Rasch case
				A <- diag(rep(1,dimensions))
				m <- sv
				
			} else {
				##   This parameterization is used for all other models
				##   although it differs depending on the specified dilation
				if (dilation=="ODL") {
					A <- matrix(sv[1:(dimensions^2)],dimensions,dimensions)
					m <- sv[(dimensions^2+1):length(sv)]
					
				##   For the LL and MIN approaches, the starting values need
				##   to be rotated using the orthogonal rotation matrix T
				##   prior to transforming any of the parameters. This makes it
				##   so that parameters only capture changes in variability
				} else if (dilation=="LL") {
					A <- diag(rep(sv[1],dimensions))%*%T
					m <- sv[2:length(sv)]
					
				} else if (dilation=="MIN") {
					A <- diag(sv[1:dimensions])%*%T
					m <- sv[(dimensions+1):length(sv)]
				} 
				
			}
			
			##   Transform the FROM scale parameters onto the TO scale
			##   using the linking constants A and m
			from.t@a <- from@a%*%ginv(A)
			if (dilation=="ODL") {
				from.t@b <- from@b-matrix(from.t@a %*% m, nrow(from@b),ncol(from@b))
			} else {
				from.t@b <- from@b-matrix(from@a %*% ginv(T) %*% m, nrow(from@b),ncol(from@b))
			}
			
			##   Transform the TO scale parameters onto the FROM scale
			##   using the linking constants A and m
			if (symmetric==TRUE) {
				to.f@a <- to@a %*% A
				if (dilation=="ODL") {
					to.f@b <- to@b+matrix(to@a %*% ginv(A) %*% m, nrow(to@b),ncol(to@b))
				} else {
					to.f@b <- to@b+matrix(to@a %*% ginv(T) %*% m, nrow(to@b),ncol(to@b))
				}
				
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
		
			##   Adjust {scr} for the MCM so that the 'do not know' category has a weight of zero
			if (tmp.to$p.mod[h]=="mcm") {
				scr <- c(scr,seq(0,p.cat[h]-1))
			} else {
				scr <- c(scr,seq(1,p.cat[h]))
			}
		}
		
		##   Use this if one of the default values (1 or 2) is used for the {score} argument
		if (length(sc)==1) {
		
			##   Use this if the lowest category should have a scoring weight of zero
			if (sc==1) {
				scr <- scr-1
				
				##   If the argument {incorrect} equals FALSE, either explicitly or
				##   by not including the argument, there will only be one column of 
				##   probabilities for each dichotomous item. This formulation of the
				##   scoring function will set the weight for these columns equal to 
				##   zero. As such, set the weights for these items equal to one
				cat <- rep(p.cat,p.cat)
				scr[cat==1] <- 1
				
				##   Make sure MCM 'do not know' categories still have a weight of zero
				scr[scr<0] <- 0
			}
			
			
			
		##   Use this if the researcher supplies a vector of score weights for
		##   all of the columns in p
		} else {
			if (length(sc)==length(scr)) {
				scr <- sc
			} else {
				warning("The length of {score} does not match the number of response probabilities. Score was set to 1")
			}
		}
		
		##   Compute sum of squares difference for the 
		##   Haebara and Stocking-Lord methods
		W1 <- as.vector(weights.t[[2]])
		W2 <- as.vector(weights.f[[2]])
		if (transform=="HB") {
			L1 <- ncol(prob.to)*sum(W1)
			Q1 <- W1*(prob.to-prob.from.t)^2
			Q1 <- sum(apply(Q1,1,sum))
			if (symmetric==TRUE) {
				L2 <- ncol(prob.to)*sum(W2)
				Q2 <- W2*(prob.from-prob.to.f)^2
				Q2 <- sum(apply(Q2,1,sum))
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
			
			##   Check to see if X'X is positive definite
			tmp <- eigen(t(from) %*% from)$values
			if (min(tmp)<0) warning("The rotation matrix is not positive definite. Be hesitant in trusting these results")
		}
		
		TF <- T
		mn <- sum((to-from%*%ginv(T))^2)
		##   Due to indeterminacies in the sign of eigenvectors, 
		##   adjust the rotation matrix to get the matrix that truly minimizes the 
		##   difference between the two sets of slopes
		tmp <- sum((to-from%*%ginv(T*-1))^2)
		if (mn>tmp) {
			mn <- tmp
			T <- T*-1
		}
		
		return(T)
	}
	
	
	###########################################################
	##             Functions used to compute response probabilities
	###########################################################
	
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
	
	
	########################################################
	##             Function used to computes descriptive statistics 
	##                        for the common item parameters
	########################################################
	
	.Descriptives <- function(a1, a2, b1, b2, c1, c2, pm, cat, dimensions) {
	
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
				colnames(des) <- c("a","b","(-b/a)","c")
				mcnr <- c(mcnr, pm.it[[j]])
			} else if (pm.mod[j]=="nrm") {
				des <- data.frame(cbind(a,b,bd))
				colnames(des) <- c("a","b","(-b/a)")
				mcnr <- c(mcnr, pm.it[[j]])
			} else {
				des <- data.frame(cbind(a,b))
				colnames(des) <- c("a","b")
			}
			rownames(des) <- c("N Pars:","Mean: To","Mean: From","SD: To","SD: From")
			descrip[[j]] <- round(des,6)
		}
		
		##   Identify the number of each type of parameter, and compute the 
		##   means and SDs of each parameter type for both groups
		a <- c(length(a1[!is.na(a1)]),mean(a1,na.rm=TRUE),mean(a2,na.rm=TRUE),.sd(a1),.sd(a2))
		b <- c(length(b1[!is.na(b1)]),mean(b1,na.rm=TRUE),mean(b2,na.rm=TRUE),.sd(b1),.sd(b2))
		if (length(c1[!is.na(c1)])) {
			c <- c(length(c1[!is.na(c1)]),mean(c1,na.rm=TRUE),mean(c2,na.rm=TRUE),.sd(c1),.sd(c2))
		}
		
		##   Compile these descriptives for all "included" common items
		if (c[1]>0) {
			desall <- round(cbind(a,b,c),6)
			colnames(desall) <- c("a","b","c")
		} else {
			desall <- round(cbind(a,b),6)
			colnames(desall) <- c("a","b")
		}
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
			if (length(pm.it[[j]])==1) {
				a1k <- t(a1k)
				a2k <- t(a2k)
			}
			
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
			
			##   If there is only one common dimension, eliminate the dummy slope from the descriptives
			if (dim.flag==TRUE) des <- des[,-2]  
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
		if (length(c1[!is.na(c1)])) {
			c <- c(length(c1[!is.na(c1)]),mean(c1,na.rm=TRUE),mean(c2,na.rm=TRUE),.sd(c1),.sd(c2))
		}
		mdisc1 <- sqrt(apply(a1^2,1,sum,na.rm=TRUE))
		mdisc2 <- sqrt(apply(a2^2,1,sum,na.rm=TRUE))
		mdif1 <- -b1/mdisc1
		mdif2 <- -b2/mdisc2
		mdc <- c(length(mdisc1[!is.na(mdisc1)]),mean(mdisc1,na.rm=TRUE),mean(mdisc2,na.rm=TRUE),.sd(mdisc1),.sd(mdisc2))
		mdf <- c(length(mdif1[!is.na(mdif1)]),mean(mdif1,na.rm=TRUE),mean(mdif2,na.rm=TRUE),.sd(mdif1),.sd(mdif2))
		
		if (c[1]>0) {
			desall <- round(cbind(a,b,c,mdc,mdf),6)
			colnames(desall) <- c(paste("a",1:dimensions,sep=""),"d","c","MDISC","MDIF")
		} else {
			desall <- round(cbind(a,b,mdc,mdf),6)
			colnames(desall) <- c(paste("a",1:dimensions,sep=""),"d","MDISC","MDIF")
		}
		rownames(desall) <- c("N Pars:","Mean: To","Mean: From","SD: To","SD: From")
		descrip[[length(descrip)+1]] <- desall
		names(descrip) <- c(pm.mod,"all")
		}
	
	return(descrip)
	}
	
	
	
	############################################
	##                        START PLINK
	############################################
	
	##   Number of groups
	ng <- x@groups
	
	##   Make {score} a list with length equal to the number of groups minus one
	if(!is.list(score)) {
		tmp <- vector("list",ng-1)
		for (i in 1:(ng-1)) {
			tmp[[i]] <- score
		}
		score <- tmp
	} else {
		if (length(score)!=(ng-1)) {
			tmp <- vector("list",ng-1)
			for (i in 1:(ng-1)) {
				tmp[[i]] <- 1
			}
			score <- tmp
			
			warning("The number of list elements in {score} does not correspond to the number of groups minus one. Score was set to 1 for all cases")
		}
	}
	
	##   Make {startvals} a list with length equal to the number of groups minus one
	if (missing(startvals)) {
		##   Make a list with NULL values
		##   Default starting values will be used
		startvals <- vector("list",ng-1)
	} else {
		if(!is.list(startvals)) {
			tmp <- vector("list",ng-1)
			for (i in 1:(ng-1)) {
				tmp[[i]] <- startvals
			}
			startvals <- tmp
		} else {
			if (length(startvals)!=(ng-1)) {
				startvals <- vector("list",ng-1)
				warning("The number of list elements in {startvals} does not correspond to the number of groups minus one. The default values were used.")
			}
		}
	}
	
	##   Maximum number of dimensions across groups
	md <- max(x@dimensions)
	
	##   Create group names (if necessary)
	 if (missing(grp.names)) grp.names <- names(x@pars)
	
	##   Create or recode the values in {dim.order}. The object
	##   dim.order.RM is the originally specified dim.order
	
	##   The columns correspond to the total number of dimensions across groups
	##   The rows correspond to the groups
	
	##   If there is more than one dimension for any group
	if (md>1) {
		if (is.null(dim.order)) {
			dim.order.RM <- NULL
			
			##   Initialize an object to identify the common dimensions between groups
			dim.order <- matrix(NA,ng,md)
			
			##   Loop through all of the groups
			for (i in 1:ng) {
				##   Set the default ordering of factors to be the same as 
				##   the ordering of supplied slope parameters. If there are
				##   different numbers of dimensions in different groups
				##   the right-most unmatched columns will be treated as the 
				##   unique factors for the given group
				dim.order[i,1:x@dimensions[i]] <- 1:x@dimensions[i]
			}
			
		##   If the {dim.order} was specified
		} else {
			##   If ones are used as placeholders for common dimensions
			##   create a set of incremental values 
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
	
	##   Check to see if {exclude} is formatted properly (if applicable)
	if (!missing(exclude)) {
		if (is.list(exclude)) {
			if (length(exclude)!=ng) stop("The {exclude} argument must be a list with length equal to the number of groups")
		}
	}
	
	
	##   Initialize an object to store the common item
	##   parameters for each pair of adjacent groups
	com <- vector("list",ng-1)
	
	##   Initialize an object to store the number of response categories 
	##   for the common items for each pair of adjacent groups
	catci <- vector("list",ng-1)
	
	##   Initialize an object to store the poly.mod objects for
	##   the common items for each pair of adjacent groups
	pm <- vector("list",ng-1)
	
	##   Loop through all of the groups (minus 1)
	for (i in 1:(ng-1)) {
	
		##   If there are only two groups total
		if (ng==2) {
			##   Get the common item matrix
			it.com <- x@common
			
		##   If there are more than two groups
		} else if (ng>2) {
			##   Get the common item matrix for the given pair of groups
			it.com <- x@common[[i]]
		}
		
		##   Modify {it.com} to account for any excluded items from {exclude}
		if (!missing(exclude)) {
			
			##   Exclude all items associated with the specified models
			if (is.character(exclude)) {
				for (j in exclude) {
					items <- eval(parse(text=paste("x@poly.mod[[i]]@items$",j,sep="")))
					it.com <- it.com[(it.com[,1] %in% items)==FALSE,]
				}
			} else {
				##   Items associated with the given models should be excluded
				##   for the given pair of tests based on {exclude} for the lower group
				if (is.character(exclude[[i]])) {
					tmp <- suppressWarnings(as.numeric(exclude[[i]]))
					items.char <- exclude[[i]][is.na(tmp)]
					for (j in items.char) {
						items <- eval(parse(text=paste("x@poly.mod[[i]]@items$",j,sep="")))
						it.com <- it.com[(it.com[,1] %in% items)==FALSE,]
					}
					
					##   Check to see if there are additional items that should be excluded
					items.num <- tmp[!is.na(tmp)]
					it.com <- it.com[(it.com[,1] %in% items.num)==FALSE,]
				} else {
					it.com <- it.com[(it.com[,1] %in% exclude[[i]])==FALSE,]
				}
				
				##   Items associated with the given models should be excluded
				##   for the given pair of tests based on {exclude} for the lower group
				if (is.character(exclude[[i+1]])) {
					tmp <- suppressWarnings(as.numeric(exclude[[i+1]]))
					items.char <- exclude[[i+1]][is.na(tmp)]
					for (j in items.char) {
						items <- eval(parse(text=paste("x@poly.mod[[i+1]]@items$",j,sep="")))
						it.com <- it.com[(it.com[,2] %in% items)==FALSE,]
					}
					
					##   Check to see if there are additional items that should be excluded
					items.num <- tmp[!is.na(tmp)]
					it.com <- it.com[(it.com[,2] %in% items.num)==FALSE,]
				} else {
					it.com <- it.com[(it.com[,2] %in% exclude[[i+1]])==FALSE,]
				}
			}
			
			if (ng==2) {
				x@common <- it.com
			} else {
				x@common[[i]] <- it.com
			}
		}
		
		
		##   Initialize a list to store the common item parameters
		##   for each group of the adjacent paired groups
		com[[i]] <- vector("list",2)
		
		##   Extract the common item parameters for the lower group
		com[[i]][[1]] <- x@pars[[i]][it.com[,1],]
		if (is.vector(com[[i]][[1]])) com[[i]][[1]] <- t(com[[i]][[1]][!is.na(com[[i]][[1]])]) 
		
		##   Extract the common item parameters for the higher group
		com[[i]][[2]] <- x@pars[[i+1]][it.com[,2],]
		if (is.vector(com[[i]][[2]])) com[[i]][[2]] <- t(com[[i]][[2]][!is.na(com[[i]][[2]])]) 
		
		##   Extract a vector of the number of response categories 
		##   for the common items for the given pair of groups
		catci[[i]] <- x@cat[[i]][it.com[,1]]
		
		
		##   The ordering of parameters in com[[i]] should be such that
		##   the first element is the "To" set and the second is the "From" set
		##   In cases where the lower group comes before the base.group
		##   the ordering of the list elements must be switched. For example
		##   if there are two groups, G1 and G2 and G2 is the base group
		##   we want to estimate constants to put the G1 parameters on the 
		##   G2 scale. When we initially populate the object {com}, the first
		##   element includes the parameters for G1 and the second for G2.
		##   because G2 is the "To" scale, these elements need to be 
		##   re-ordered so that the first list element includes the parameters
		##   for G2 and the second includes the parameters for G1
		
		##   Perform the switch (if necessary)
		if (i<base.grp) com[[i]] <- com[[i]][c(2,1)] 
		
		##   Identify all of the item response models used 
		##   (for both common and unique items)
		mod <-  x@poly.mod[[i]]@model
		
		##   Identify the common and unique items
		##   associated with each item response model
		items <- x@poly.mod[[i]]@items
		
		##   Extract the common item numbers for the "To" group to 
		##   facilitate the creation of a poly.mod object for the common items
		if (ng==2) {
			it.com <- x@common[,1] 
		} else {
			it.com <- x@common[[i]][,1] 
		}
		
		##   Initialize an object to store the item response
		##   models used for the common items
		mod1 <- NULL
		
		##   Initialize an object to store the item numbers associated with
		##   the item response models used for the common items
		poly <- vector("list",5)
	
		##   Identify the common items associated with each model
		##   Not all models need to have common items
		
		##   Initialize an object to increment with the 
		##   included item response models
		step <- 1
		for (k in 1:length(mod)) {
			##  Determine the number of common items
			##   associated with the given model
			tmp <- seq(1,length(it.com))[it.com %in% items[[k]]]
			
			##   If no common items are associated with this model
			if (length(tmp)==0) {
				next 
				
			##   If there are common items are associated with this model
			}else {
				##   Add the given model
				mod1 <- c(mod1,mod[k])
				
				##   Add the items associated with the given model
				poly[[k]] <- tmp
				#poly[[step]] <- tmp
				step <- step+1
			}
		}
		
		poly <- poly[unlist(lapply(poly,length))>0]
		# Create poly.mod object for the common items
		pm[[i]] <- as.poly.mod(length(it.com),mod1,poly) }
		
		
	##   Initialize an object to store the linking constants
	##   and common item descriptives for each pair
	##   of adjacent groups
	link.out <- vector("list",ng-1)
	
	##   Separate the common item parameters to get the specific parameters
	##   Loop through all of the groups
	for (i in 1:(ng-1)) {
	
		##   The "To" parameters
		tmp1 <- sep.pars(com[[i]][[1]],catci[[i]],pm[[i]],x@dimensions[i],x@location[i])
		
		##   The "From" parameters
		tmp2 <- sep.pars(com[[i]][[2]],catci[[i]],pm[[i]],x@dimensions[i+1],x@location[i+1])
		
		##   The maximum number of dimensions across both groups
		##   (the dimensions can differ from group to group)
		tmp.md <- max(c(x@dimensions[i],x@dimensions[i+1]))
		
		##############################################
		##                Identify common dimensions
		##############################################
		
		##   Create a flag to identify whether (in the multidimensional
		##   case) there is only one common dimension
		dim.flag <- FALSE
		
		##   If there are any groups with more than one dimension
		if (tmp.md>1) {
			
			##   Identify the common dimensions in the lower group
			do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
			
			##   Identify the common dimensions in the higher group
			do2 <- dim.order[i+1,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
			
			##   Extract the slope parameters for the common dimensions
			tmp1@a <- as.matrix(as.matrix(tmp1@a)[,do1])
			tmp2@a <- as.matrix(as.matrix(tmp2@a)[,do2])
			
			##   Set the number of dimensions for both groups equal
			##   to the number of common dimensions
			tmp1@dimensions <- tmp2@dimensions <- dimensions <- length(do1)
			
			##   Check to see if the parameters for a given pair of groups
			##   are unidimensional (i.e., parameterized using a slope/difficulty
			##   rather than a slope/intercept as commonly used in the MD case)
			for (j in i:(i+1)) {
				##   If the parameters for a given group use a slope/difficulty parameterization,
				##   reformat them to use a slope/intercept parameterization
				if (x@dimensions[j]==1) {
					for (k in 1:length(pm[[i]]@model)) {
						if (pm[[i]]@model[k]=="drm"|pm[[i]]@model[k]=="grm"|pm[[i]]@model[k]=="gpcm") {
							tmp1@b[pm[[i]]@items[[k]],] <- as.matrix(tmp1@b[pm[[i]]@items[[k]],]*-as.vector(tmp1@a)[pm[[i]]@items[[k]]])
							tmp2@b[pm[[i]]@items[[k]],] <- as.matrix(tmp2@b[pm[[i]]@items[[k]],]*-as.vector(tmp2@a)[pm[[i]]@items[[k]]])
						}
					}
				}
			}
							
			##   If there is only one common dimension add a vector of 
			##   zeros to the matrix of slope parameters so that the MD
			##   linking methods can be used (this should not affect the estimation)
			if (dimensions==1) {
				tmp1@a <- cbind(tmp1@a,0)
				tmp2@a <- cbind(tmp2@a,0)
				tmp1@dimensions <- tmp2@dimensions <- dimensions <- 2
				dim.flag <- TRUE
			}
		} else {
			dimensions <- 1
		}
		
		
		##   Extract the common item parameters for each group
		a1 <- tmp1@a
		a2 <- tmp2@a
		b1 <- tmp1@b
		b2 <- tmp2@b
		c1 <- tmp1@c
		c2 <- tmp2@c
		
		##   Identify the weights to be used for the given group for the TO scale
		if (missing(weights.t)) {
			if (dimensions>2) {
				s <- matrix(.6,dimensions,dimensions)
				diag(s) <- rep(1,dimensions)
				th <- mvrnorm(1000,rep(0,dimensions),s)
				colnames(th) <- paste("theta",1:dimensions,sep="")
				wt <- 1
				for (j in 1:dimensions) {
					wt <- wt*dnorm(th[,j])
				}
				wgt.t <- list(points=th,weights=wt)
			} else {
				wgt.t <- as.weight(dimensions=dimensions)
			}
		} else {
			if (is.list(weights.t[[1]])) wgt.t <- weights.t[[i]] else wgt.t <- weights.t
		}
		
		##   Identify the weights to be used for the given group for the FROM scale
		if (missing(weights.f)) {
			wgt.f <- wgt.t 
		} else {
			if (is.list(weights.f[[1]])) wgt.f <- weights.f[[i]] else wgt.f <- weights.f
		}
		
		##   Use all methods (if none are specified)
		if (missing(method)) {
			if (dimensions==1) method <- c("MM","MS","HB","SL") else method <- "LS"
		} else {
			method <- toupper(method)
			if (sum(method%in% c("MM","MS","HB","SL","LS") )==0) {
				warning("No appropriate method was selected. All default methods will be used")
				if (dimensions==1) method <- c("MM","MS","HB","SL") else method <- "LS"
			}
		}
		
		##   Compute the descriptive statistics for reporting
		descrip <- .Descriptives(a1,a2,b1,b2,c1,c2,pm[[i]],catci[[i]],dimensions)
		
		##   Check to see if the {rescale} method is included in the {method} argument
		##   If not, add it to {method}
		if (!missing(rescale)) {
			if ((toupper(rescale)%in%method)==FALSE) {
				method <- c(method,toupper(rescale))
			}
		}
		
		
		########################################
		##            Perform the Calibration
		########################################
		
		##   Initialize an object to store the linking constants for each method
		constants <- list(NULL)
		
		##   Initialize objects to store the iteration, convergence, and criterion
		##   value information for the characteristic curve methods
		it <- con <- obj <- NULL
		
		##   Unidimensional case
		if (dimensions==1) {
		
			##   Initialize an object to check whether the Rasch model
			##   or PCM (with slopes equal to 1) is being used
			tmp <- c(a1,a2)
			if (length(tmp[tmp==1])==length(tmp)) {
				rasch.flag <- TRUE
			} else {
				rasch.flag <- FALSE
			}
			
			##   Compute linking constants using the moment methods
			mm <- ms <- NULL
			
			##   Mean/Mean - A
			A1 <- descrip$all[3,1]/descrip$all[2,1]
			
			##   Mean/Sigma  - A
			A2 <- descrip$all[4,2]/descrip$all[5,2]
			
			##   If the Rasch model or PCM is used, set the A2 constant equal to A1
			if (rasch.flag==TRUE) A2 <- A1
			
			##   Mean/Mean - B
			B1 <- descrip$all[2,2]-A1*descrip$all[3,2]
			
			##   Mean/Sigma - B
			B2 <- descrip$all[2,2]-A2*descrip$all[3,2]
			
			##   Format the constants for the moment methods
			##   to be included in the output
			mm <- round(c(A1,B1),6)
			ms <- round(c(A2,B2),6)
			names(mm) <- names(ms) <- c("A","B")
			
			##   Add these constants to the output
			if ("MM" %in% method) constants$MM <- mm
			if ("MS" %in% method) constants$MS <- ms
			
			##   These values are only needed for the multidimensional case, but
			##   they are required for the criterion function for the characteristic curve methods
			dilation <- "N/A"
			T <- NA
			
			##   Determine starting values for the characteristic curve methods
			if (is.null(startvals[[i]])) {
				##   Use the Mean/Sigma values
				sv <- mm
			} else {
				if (is.character(startvals[[i]])) {
					if (toupper(startvals[[i]])=="MM") sv <- mm
					if (toupper(startvals[[i]])=="MS") sv <- ms
				} else {
					sv <- startvals[[i]]
				}
			}
			
			if (rasch.flag==TRUE) {
				##   Only include starting values for the B constant
				sv <- sv[2]
			}
			
		##   Multidimensional case
		} else {
		
			##   Initialize an object to check whether the multidimensional 
			##   Rasch model or MPCM (with slopes equal to 1) is being used
			tmp <- as.vector(a1)
			if (length(tmp[tmp==1])==length(tmp)) {
				rasch.flag <- TRUE
			} else {
				rasch.flag <- FALSE
			}
			
			##   Reformat b1 and b2 as vectors.  This is necessary
			##   when there are polytomous items
			tmp.b1 <- as.vector(b1)[!is.na(as.vector(b1))]
			tmp.b2 <- as.vector(b2)[!is.na(as.vector(b2))]
			
			##   When transforming b1 and b2 into vectors
			##   there will be more rows than in the original matrices
			##   if there are polytomous items.  As such, the matrix
			##   of slope parameters needs to be adjusted so that
			##   appropriate slopes are matched up with each of 
			##   the reformatted b parameters.  For the LS methods
			##   it is only necessary to re-specify the matrix of 
			##   slopes for the "From" group
			tmp.a2 <- NULL
			for (j in 1:dimensions) {
				tmp.a2a <- (!is.na(b2))*a2[,j]
				tmp.a2a[is.na(b2)] <- NA
				tmp.a2a <- as.vector(tmp.a2a)[!is.na(as.vector(tmp.a2a))]
				tmp.a2 <- cbind(tmp.a2,tmp.a2a)
			}
			
			##   Compute linking constants using the least squares (LS) method (oblique or orthogonal)
			##   Create the rotation matrix
			tmp <- as.vector(a1)
			if (rasch.flag==TRUE) {
				A <- diag(rep(1,dimensions))
			} else {
				##   Oblique procrustes method
				if (dilation=="ODL") {
					# This actually estimates the inverse of A
					A <- .Rotate(a1,a2,FALSE)
					
					##   Compute the translation vector 
					y <- tmp.b2-tmp.b1 
					X <- tmp.a2%*%A
					m <- ginv(t(X)%*%X)%*%t(X)%*%y
					
					m <- as.vector(m)
					A <- ginv(A)
					rownames(A) <- rep("",dimensions)
					colnames(A) <- c(rep("",dimensions-1),"A")
					names(m) <- paste("m",1:dimensions,sep="")
					
					if ("LS" %in% method) {
						##   When there are two or more common dimensions between the groups
						if (dim.flag==FALSE) {
							constants$LS <- list(A=round(A,6),m=round(m,6))
							
						##   When there is only one common dimension between the groups
						} else {
							constants$LS <- c(round(A[1,1],6),round(m[1],6))
							names(constants$LS) <- c("A","m")
						}
					}
					
				##   Orthogonal procrustes method
				} else {
				
					a1c <- a1-matrix(apply(a1,2,mean,na.rm=TRUE),nrow(a1),ncol(a1),byrow=TRUE)
					a2c <- a2-matrix(apply(a2,2,mean,na.rm=TRUE),nrow(a2),ncol(a2),byrow=TRUE)
					T <- .Rotate(a1,a2,TRUE)
					
					##   Estimate the translation vector 
					y <- tmp.b2-tmp.b1 
					X <- tmp.a2%*%T
					m <- ginv(t(X)%*%X)%*%t(X)%*%y
					
					##   Estimate the dilation parameter(s)
					if (dilation=="LL") {
						k <- sum(diag(t(T)%*%t(a2c)%*%a1c))/sum(diag(t(a2c)%*%a2c))
						K <- diag(rep(k,dimensions))
					} else if (dilation=="MIN") {
						K <- diag(diag(t(a1c)%*%a2c%*%T))%*%solve(diag(diag(t(T)%*%t(a2c)%*%a2c%*%T)))
					}
					
					K <- ginv(K)
					m <- as.vector(m)
					rownames(T) <- rownames(K) <- rep("",dimensions)
					colnames(T) <- c(rep("",dimensions-1),"T")
					colnames(K) <- c(rep("",dimensions-1),"K")
					names(m) <- paste("m",1:dimensions,sep="")
					
					if ("LS" %in% method) {
						##   When there are two or more common dimensions between the groups
						if (dim.flag==FALSE) {
							constants$LS <- list(T=round(T,6),K=round(K,6),m=round(m,6))
							
						##   When there is only one common dimension between the groups
						} else {
							constants$LS <- c(round(T[1,1]*K[1,1],6),round(m[1],6))
							names(constants$LS) <- c("A","m")
						}
					}
				}
			}
			
			
			if ("HB" %in% method | "SL" %in% method) {
				##   Determine starting values for the characteristic curve methods
				if (rasch.flag==TRUE) {
					if (is.null(startvals[[i]])) {
						##   Use the values from the LS method
						sv <- m 
					} else {
						sv <- startvals[[i]][((length(startvals[[i]])-dimensions)+1):length(startvals[[i]])]
					}
				} else {
					if (is.null(startvals[[i]])) {
					
						##   Oshima, Davey, & Lee method
						if (dilation=="ODL") {
							##   Use the values from the LS method
							sv <- c(as.vector(A),m) 
							
						##   Li & Lissitz method
						} else if (dilation=="LL") {
							##   Use the values from the LS method
							sv <- c(K[1,1],m)
							
						##   Min method
						} else if (dilation=="MIN") {
							##   Use the values from the LS method
							sv <- c(diag(K),m)
						
						} 
					} else {
					
						sv <- startvals[[i]]
						
						##   Check to see if the correct number of starting values have been specified
						if (dilation=="ODL") {
							if ((length(sv))!=dimensions+dimensions^2) {
								sv <- c(as.vector(A),m) 
								warning(paste("The number of elements in {startvals} for groups",i,"and",i+1,"does not equal",dimensions+dimensions^2,". The default values were used>"))
							}
						} else if (dilation=="LL") {
							if ((length(sv))!=dimensions+1) {
								sv <- c(K[1,1],m)
								warning(paste("The number of elements in {startvals} for groups",i,"and",i+1,"does not equal",dimensions+1,". The default values were used>"))
							}
						} else if (dilation=="MIN") {
							if ((length(sv))!=dimensions*2) {
								sv <- c(diag(K),m)
								warning(paste("The number of elements in {startvals} for groups",i,"and",i+1,"does not equal",dimensions*2,". The default values were used>"))
							}
						}
					}
				}
			}
		}
		
		
		##   Characteristic curve methods
		
		##   Haebara Method
		if ("HB" %in% method) {
		
			##   Estimate the linking constants
			hb <- nlminb(sv, .CC, to=tmp1, from=tmp2, dimensions=dimensions, weights.t=wgt.t, weights.f=wgt.f, sc=score[[i]], transform="HB", symmetric=symmetric, dilation=dilation, T=T, ...)
			
			##   Reformat the information from the estimation (as necessary)
			##   to prepare it for being output
			if (dimensions==1) { 
				if (rasch.flag==TRUE) hb$par <- c(1,hb$par)
				names(hb$par) <- c("A","B")
				constants$HB <- round(hb$par,6)
			} else {
				if (dilation=="ODL") {
					##   Reformat the estimated values for the rotation matrix as a matrix
					if (rasch.flag==TRUE) A <- diag(rep(1,dimensions)) else A <- matrix(hb$par[1:(dimensions^2)],dimensions,dimensions)
					
					##   Extract the constants for the translation vector
					m <- hb$par[-c(1:(dimensions^2))]
					
				} else if (dilation=="LL") {
					##   Reformat the single dilation parameter as a diagonal matrix
					if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(rep(hb$par[1],dimensions))
					
					##   Extract the constants for the translation vector
					m <- hb$par[-1]
					
				} else if (dilation=="MIN") {
					##   Reformat the dilation parameters as a diagonal matrix
					if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(hb$par[1:dimensions])
					
					##   Extract the constants for the translation vector
					m <- hb$par[-c(1:dimensions)]
				} 
				
				names(m) <- paste("m",1:dimensions,sep="")
				if (dilation=="ODL") {
					rownames(A) <-rep("",dimensions)
					colnames(A) <- c(rep("",dimensions-1),"A")
					constants$HB <- list(A=round(A,6), m=round(m,6))
				} else {
					rownames(T) <- rownames(K) <- rep("",dimensions)
					colnames(T) <- c(rep("",dimensions-1),"T")
					colnames(K) <- c(rep("",dimensions-1),"K")
					constants$HB <- list(T=round(T,6), K=round(K,6), m=round(m,6))
				}
			}
			
			##   When there is only one common dimension between the groups
			if (dim.flag==TRUE) {
				constants$HB <- c(constants$HB$A[1,1],constants$HB$m[1])
				names(constants$HB) <- c("A","m")
			}
			
			it <- c(it, HB=hb$iterations)
			con <- c(con, HB=hb$message)
			obj <- c(obj, HB=hb$objective)
		}
		
		
		##   Stocking-Lord Method
		if ("SL" %in% method) {
			
			##   Estimate the linking constants
			sl <- nlminb(sv, .CC, to=tmp1, from=tmp2, dimensions=dimensions, weights.t=wgt.t, weights.f=wgt.f, sc=score[[i]], transform="SL", symmetric=symmetric, dilation=dilation, T=T, ,...)
			
			##   Reformat the information from the estimation (as necessary)
			##   to prepare it for being output
			if (dimensions==1) { 
				if (rasch.flag==TRUE) sl$par <- c(1,sl$par)
				names(sl$par) <- c("A","B")
				constants$SL <- round(sl$par,6)
			} else {
				if (dilation=="ODL") {
					##   Reformat the estimated values for the rotation matrix as a matrix
					if (rasch.flag==TRUE) A <- diag(rep(1,dimensions)) else A <- matrix(sl$par[1:(dimensions^2)],dimensions,dimensions)
					
					##   Extract the constants for the translation vector
					m <- sl$par[-c(1:(dimensions^2))]
					
				} else if (dilation=="LL") {
					##   Reformat the single dilation parameter as a diagonal matrix
					if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(rep(sl$par[1],dimensions))
					
					##   Extract the constants for the translation vector
					m <- sl$par[-1]
					
				} else if (dilation=="MIN") {
					##   Reformat the dilation parameters as a diagonal matrix
					if (rasch.flag==TRUE) K <- diag(rep(1,dimensions)) else K <- diag(sl$par[1:dimensions])
					
					##   Extract the constants for the translation vector
					m <- sl$par[-c(1:dimensions)]
				} 
				
				names(m) <- paste("m",1:dimensions,sep="")
				if (dilation=="ODL") {
					rownames(A) <-rep("",dimensions)
					colnames(A) <- c(rep("",dimensions-1),"A")
					constants$SL <- list(A=round(A,6), m=round(m,6))
				} else {
					rownames(T) <- rownames(K) <- rep("",dimensions)
					colnames(T) <- c(rep("",dimensions-1),"T")
					colnames(K) <- c(rep("",dimensions-1),"K")
					constants$SL <- list(T=round(T,6), K=round(K,6), m=round(m,6))
				}
			}
			
			##   When there is only one common dimension between the groups
			if (dim.flag==TRUE) {
				constants$SL <- c(constants$SL$A[1,1],constants$SL$m[1])
				names(constants$SL) <- c("A","m")
			}
			
			it <- c(it, SL=sl$iterations)
			con <- c(con, SL=sl$message)
			obj <- c(obj, SL=sl$objective)
		}
		
		constants[[1]] <- NULL
		if (is.null(it)) it <- 0
		if (is.null(con)) con <- "N/A"
		if (is.null(obj)) obj <- 0
		
		##  Reset the formal arguments for nlminb
		formals(nlminb)$control <- list()
		
		##   Create labels identifying the pairs of groups for the linking constants
		##   Include an asterisk identifying the base group
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
		
		##   Create the {link} object containing the linking constants
		##   and common item descriptive statistics for the pair of adjacent groups
		link.out[[i]] <- new("link", constants=constants, descriptives=descrip, iterations=it, objective=obj, convergence=con, base.grp=base.grp, grp.names=nms, n=tmp1@n, mod.lab=tmp1@mod.lab, dilation=dilation)
		
	}
	
	########################################
	##            Rescale the Parameters
	########################################
	
	if (!missing(rescale)) {
		##   Change the method used to rescale the item parameters (if necessary)
		if ((toupper(rescale)%in%method)==FALSE) {
			if (dimensions==1) {
				if ("SL" %in% method) rescale <- "SL" else rescale <- method[length(method)]
			} else {
				if ("LS" %in% method) rescale <- "LS" else rescale <- method[length(method)]
			}
			warning(paste("No linking constants were computed for the rescale method you selected. The parameters will be rescaled using the ",rescale," method",sep=""))
		}
		
		##   Initialize a list to store the specific linking constants
		##   that will be used to rescale all of the item parameters
		tmp.con <- vector("list",ng)
		
		##   Initialize an object to increment the list element for tmp.con
		j <- 1
		
		##  Loop through all of the groups to compile the linking constants
		for (i in 1:ng) {
			##   For the base group, set the linking constants to 1 and 0 for A and B 
			##   respectively in the unidimensional case, or to an identity matrix and 
			##   a vector of zeros for the rotation matrix and translation vector in the 
			##   multidimensional case
			if (i==base.grp) {
				if (x@dimensions[i]==1) {
					tmp.con[[i]] <- c(1,0) 
				} else {
					if (dilation=="ODL") {
						tmp.con[[i]] <- list(A=diag(rep(1,x@dimensions[i])),m=rep(0,x@dimensions[i]))
					} else {
						tmp.con[[i]] <- list(T=diag(rep(1,x@dimensions[i])), K=diag(rep(1,x@dimensions[i])),m=rep(0,x@dimensions[i]))
					}
				}
			
			##   Extract the specific linking constants for the given pair of tests
			} else {
				tmp.con[[i]] <- eval(parse(text=paste("link.out[[",j,"]]@constants$",rescale,sep="")))
				j <- j+1
			}
		}
		
		##   Initialize objects to store the rescaled item and ability parameters
		out.pars <- vector("list",ng)
		out.ability <- vector("list",ng)
		
		
		##   Loop through all of the groups and rescale the item parameters
		##   and ability parameters (if applicable)
		for (i in 1:ng) {
		
			##   poly.mod object for unique and common items for a given group
			pm <- x@poly.mod[[i]]@model
			
			##   Identify NRM and MCM items
			mcnr <- NULL
			for (k in 1:length(pm)) {
				if (pm[k]=="mcm"|pm[k]=="nrm") mcnr <- c(mcnr, x@poly.mod[[i]]@items[[k]])
			}
			
			##   Separate out the item parameters for the given group
			tmp <- sep.pars(x@pars[[i]],x@cat[[i]],x@poly.mod[[i]],x@dimensions[i],x@location[i])
			
			##   Extract ability estimates for the given group (if applicable)
			if (!missing(ability)) tmpa <- ability[[i]] 
			
			##   Number of dimensions for the given group
			dimensions <- x@dimensions[i]
			
			##   Initialize a set of indexing values to loop over to
			##   place a given set of parameters on the base scale
			tmp1 <- tmp
			if (i==base.grp) {
				tmp.j <- i
			} else if (i < base.grp) {
				tmp.j <- i:(base.grp-1)
			} else if (i > base.grp) {
				tmp.j <- i:(base.grp+1)
			}
			
			##   Perform j iterations to place the parameters for the given group on the base group scale
			for (j in tmp.j) {
				if (dimensions==1) {
					tmp1@a <- as.matrix(tmp1@a/tmp.con[[j]][1])
					tmp1@b <-  as.matrix(tmp.con[[j]][1]*tmp1@b + tmp.con[[j]][2])
					tmp1@b[mcnr,] <- tmp@b[mcnr,]-(tmp.con[[j]][2]/tmp.con[[j]][1])*tmp@a[mcnr,]
					tmp@b[mcnr,] <- tmp1@b[mcnr,]
					if (!missing(ability)) tmpa <- tmp.con[[j]][1]*tmpa + tmp.con[[j]][2]
				} else {
					##   Identify the dimensions for the given group that should be transformed
					if (i==base.grp) {
						do1 <- dim.order[i,!is.na(dim.order[i,])]
						do4 <- 1:length(do1)
					} else if (i < base.grp) {
						do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[j+1,])]
						do2 <- dim.order[j+1,!is.na(dim.order[i,]) &!is.na(dim.order[j+1,])]
						do3 <- dim.order[j+1,!is.na(dim.order[j,]) & !is.na(dim.order[j+1,])]
						do4 <- c(1:length(do3))[(do3%in%do2)==TRUE]
					} else if (i > base.grp) {
						do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[j-1,])]
						do2 <- dim.order[j-1,!is.na(dim.order[i,]) &!is.na(dim.order[j-1,])]
						do3 <- dim.order[j-1,!is.na(dim.order[j,]) & !is.na(dim.order[j-1,])]
						do4 <- c(1:length(do3))[(do3%in%do2)==TRUE]
					}
					
					if (length(do1)==0) {
						next
					} else if (length(do1)==1) {
						if (is.vector(tmp.con[[j]])) {
							tmp1@a[,do1] <- tmp1@a[,do1]*ginv(tmp.con[[j]][1])
							tmp1@b <- tmp1@b-matrix(tmp1@a[,do1]*tmp.con[[j]][2],nrow(tmp1@b),ncol(tmp1@b))
							if (!missing(ability)) tmpa[,do1] <- tmp.con[[j]][1]*tmpa[,do1] + tmp.con[[j]][2]
						} else {
							if (dilation=="ODL") {
								tmp1@a[,do1] <- tmp1@a[,do1]%*%ginv(tmp.con[[j]]$A[do4,do4])
								tmp1@b <- tmp1@b-matrix(tmp1@a[,do1]%*%tmp.con[[j]]$m[do4],nrow(tmp1@b),ncol(tmp1@b))
								if (!missing(ability)) tmpa[,do1] <- t(tmp.con[[j]]$A[do4,do4]%*%t(tmpa[,do1]) + tmp.con[[j]]$m[do4])
							} else {
								tmp1a <- tmp1@a
								tmp1a[,do1] <- tmp1a[,do1]%*%ginv(tmp.con[[j]]$T[do4,do4])
								tmp1@a[,do1] <- tmp1@a[,do1]%*%ginv(tmp.con[[j]]$T[do4,do4])%*%ginv(tmp.con[[j]]$K[do4,do4])
								tmp1@b <- tmp1@b-matrix(tmp1a[,do1]%*%tmp.con[[j]]$m[do4],nrow(tmp1@b),ncol(tmp1@b))
								if (!missing(ability)) tmpa[,do1] <- t(tmp.con[[j]]$T[do4,do4]%*%tmp.con[[j]]$K[do4,do4]%*%t(tmpa[,do1]) + as.vector(tmp.con[[j]]$K[do4,do4]%*%tmp.con[[j]]$m[do4]))
							}
						}
					} else {
						if (dilation=="ODL") {
							tmp1@a[,do1] <- tmp1@a[,do1]%*%ginv(tmp.con[[j]]$A[do4,do4])
							tmp1@b <- tmp1@b-matrix(tmp1@a[,do1]%*%tmp.con[[j]]$m[do4],nrow(tmp1@b),ncol(tmp1@b))
							if (!missing(ability)) tmpa[,do1] <- t(tmp.con[[j]]$A[do4,do4]%*%t(tmpa[,do1]) + tmp.con[[j]]$m[do4])
						} else {
							tmp1a <- tmp1@a
							tmp1a[,do1] <- tmp1a[,do1]%*%ginv(tmp.con[[j]]$T[do4,do4])
							tmp1@a[,do1] <- tmp1@a[,do1]%*%ginv(tmp.con[[j]]$T[do4,do4])%*%ginv(tmp.con[[j]]$K[do4,do4])
							tmp1@b <- tmp1@b-matrix(tmp1a[,do1]%*%tmp.con[[j]]$m[do4],nrow(tmp1@b),ncol(tmp1@b))
							if (!missing(ability)) tmpa[,do1] <- t(tmp.con[[j]]$T[do4,do4]%*%tmp.con[[j]]$K[do4,do4]%*%t(tmpa[,do1]) + as.vector(tmp.con[[j]]$K[do4,do4]%*%tmp.con[[j]]$m[do4]))
						}
					}
				}
			}
			
			##   Compile the rescaled item parameters
			out.pars[[i]] <- tmp1
			names(out.pars)[[i]] <- grp.names[i]
			
			##   Compile the rescaled ability estimates
			if (!missing(ability)) {
				out.ability[[i]] <- tmpa
				names(out.ability)[[i]] <- grp.names[i]
			}
		}
		
		##   Substitute the non-scaled parameters for the common items into the set of rescaled item parameters
		if (rescale.com==FALSE) {
		
			##   Extract the common item matrix/matrices
			com <- x@common
			if (is.matrix(com)) com <- list(com)
			
			##   Loop through all of the groups lower than the base group
			for (i in (base.grp-1):1) {
				if (base.grp==1) break
				
				##   There may be different numbers of columns for the matrices
				##   of b and c parameters for the two groups. Figure out the minimum
				b.min <- min(c(ncol(out.pars[[i]]@b),ncol(out.pars[[i+1]]@b)))
				c.min <- min(c(ncol(out.pars[[i]]@c),ncol(out.pars[[i+1]]@c)))
				
				##   Replace the b and c parameters for the group further from the base 
				##   group with the parameters from the group closer to the base group
				##   (i.e., the "To" scale parameters in each of the adjacent pairs used
				##   when estimating the linking constants)
				out.pars[[i]]@b[com[[i]][,1],1:b.min] <- out.pars[[i+1]]@b[com[[i]][,2],1:b.min]
				out.pars[[i]]@c[com[[i]][,1],1:c.min] <- out.pars[[i+1]]@c[com[[i]][,2],1:c.min]
				
				##   Replace the slope parameters for the group further from the base 
				##   group with the parameters from the group closer to the base group
				if (dimensions==1) {
					out.pars[[i]]@a[com[[i]][,1],] <- out.pars[[i+1]]@a[com[[i]][,2],]
				} else {
					##   Make the replacement only for the common dimensions
					do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					do2 <- dim.order[i+1,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					out.pars[[i]]@a[com[[i]][,1],do1] <- out.pars[[i+1]]@a[com[[i]][,2],do2]
					
				}
			}
			
			##   Loop through all of the groups higher than the base group
			for (i in base.grp:(ng-1)) {
				if (base.grp==ng) break
				
				##   There may be different numbers of columns for the matrices
				##   of b and c parameters for the two groups. Figure out the minimum
				b.min <- min(c(ncol(out.pars[[i]]@b),ncol(out.pars[[i+1]]@b)))
				c.min <- min(c(ncol(out.pars[[i]]@c),ncol(out.pars[[i+1]]@c)))
				
				##   Replace the b and c parameters for the group further from the base 
				##   group with the parameters from the group closer to the base group
				##   (i.e., the "To" scale parameters in each of the adjacent pairs used
				##   when estimating the linking constants)
				out.pars[[i+1]]@b[com[[i]][,2],1:b.min] <- out.pars[[i]]@b[com[[i]][,1],1:b.min]
				out.pars[[i+1]]@c[com[[i]][,2],1:c.min] <- out.pars[[i]]@c[com[[i]][,1],1:c.min]
				
				##   Replace the slope parameters for the group further from the base 
				##   group with the parameters from the group closer to the base group
				if (dimensions==1) {
					out.pars[[i+1]]@a[com[[i]][,2],] <- out.pars[[i]]@a[com[[i]][,1],]
				} else {
					##   Make the replacement only for the common dimensions
					do1 <- dim.order[i,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					do2 <- dim.order[i+1,!is.na(dim.order[i,]) &!is.na(dim.order[i+1,])]
					out.pars[[i+1]]@a[com[[i]][,2],do2] <- out.pars[[i]]@a[com[[i]][,1],do1]
				}
			}
			
		}
		
	} else {
		if (!missing(ability)) {
			##   Change the method used to rescale the item parameters (if necessary)
			if (dimensions==1) {
				if ("SL" %in% method) rsc <- "SL" else rsc<- method[length(method)]
			} else {
				if ("LS" %in% method) rsc <- "LS" else rsc<- method[length(method)]
			}
			cat(paste("No rescale method was identified. Ability parameters will be rescaled using the ",rsc," method.\n",sep=""))
			
			
			##   Initialize a list to store the specific linking constants
			##   that will be used to rescale all of the ability estimates
			tmp.con <- vector("list",ng)
			
			##   Initialize an object to increment the list element for tmp.con
			j <- 1
			
			##  Loop through all of the groups to compile the linking constants
			for (i in 1:ng) {
				##   For the base group, set the linking constants to 1 and 0 for A and B 
				##   respectively in the unidimensional case, or to an identity matrix and 
				##   a vector of zeros for the rotation matrix and translation vector in the 
				##   multidimensional case
				if (i==base.grp) {
					if (dimensions==1) {
						tmp.con[[i]] <- c(1,0) 
					} else {
						if (dilation=="ODL") {
							tmp.con[[i]] <- list(A=diag(rep(1,dimensions)),m=rep(0,dimensions))
						} else {
							tmp.con[[i]] <- list(T=diag(rep(1,dimensions)),K=diag(rep(1,dimensions)),m=rep(0,dimensions))
						}
					}
				} else {
					##   Extract the specific linking constants for the given pair of tests
					tmp.con[[i]] <- eval(parse(text=paste("link.out[[",j,"]]@constants$",rsc,sep="")))
					j <- j+1
				}
			}
			
			##   Initialize an object to store the rescaled ability estimates
			out.ability <- vector("list",ng)
			
			##   Loop through all of the groups and rescale the ability estimates
			for (i in 1:ng) {
				
				##   Number of dimensions for the given group
				dimensions <- x@dimensions[i]
				
				##   Initialize a set of indexing values to loop over to
				##   place a given set of parameters on the base scale
				tmpa <- ability[[i]]
				if (i==base.grp) {
					tmp.j <- i
				} else if (i < base.grp) {
					tmp.j <- i:(base.grp-1)
				} else if (i > base.grp) {
					tmp.j <- i:(base.grp+1)
				}
				
				##   Perform j iterations to place the ability estimates for the given group on the base group scale
				for (j in tmp.j) {
					
					if (dimensions==1) {
						tmpa <- tmp.con[[j]][1]*tmpa + tmp.con[[j]][2]
					} else {
						##   Identify the dimensions for the given group that should be transformed
						if (i==base.grp) {
							do1 <- dim.order[j,!is.na(dim.order[j,])]
						} else if (i < base.grp) {
							do1 <- dim.order[j,!is.na(dim.order[j,]) &!is.na(dim.order[j+1,])]
						} else if (i > base.grp) {
							do1 <- dim.order[j,!is.na(dim.order[j,]) &!is.na(dim.order[j-1,])]
						}
						
						if (length(do1)==1) {
							tmpa[,do1] <- tmp.con[[j]][1]*tmpa[,do1] + tmp.con[[j]][2]
						} else {
							if (dilation=="ODL") {
								tmpa[,do1] <- t(tmp.con[[j]]$A%*%t(tmpa[,do1]) + tmp.con[[j]]$m)
							} else {
								tmpa[,do1] <- t(tmp.con[[j]]$T%*%tmp.con[[j]]$K%*%t(tmpa[,do1]) + as.vector(tmp.con[[j]]$K%*%tmp.con[[j]]$m))
							}
						}
					}
				}
				
				##   Compile the rescaled ability estimates
				out.ability[[i]] <- tmpa
				names(out.ability)[[i]] <- grp.names[i]
			}
		}
	}
	
	##   Include a location parameter if the original set of item parameters had one
	if (!missing(rescale)) {
		for (i in 1:length(out.pars)) {
			out.pars[[i]] <- sep.pars(as.irt.pars(out.pars[[i]]),loc.out=x@location[i])
		}
	}
	
	
	##   Combine the {link} object and rescaled item parameters/ability estimates (as necessary) to be output
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
