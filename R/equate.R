##   This function conducts IRT true-score and observed-score equating

setGeneric("equate", function(x, method="TSE", true.scores, ts.low=TRUE, base.grp=1, score=1, startval, weights1, weights2, syn.weights, ...) standardGeneric("equate"))



setMethod("equate", signature(x="list"), function(x, method, true.scores, ts.low, base.grp, score, startval, weights1, weights2, syn.weights, ...) {

	##   Extract the rescaled item parameters from the object output by {plink}
	if (length(x$pars)) {
		x<- x$pars
		callGeneric()
	} else {
		stop("There were no parameters in {x}, re-run plink and specify and argument for {rescale} then try again.")
	}
	
})



setMethod("equate", signature(x="irt.pars"), function(x, method, true.scores, ts.low, base.grp, score, startval, weights1, weights2, syn.weights, ...) {


	##      Function that will be minimized for the True Score Equating
	.TSE <- function(startval, truescore, pars, sc, ...) {
		
		##   Compute response probabilities for the base group first
		prob <- .Mixed(pars, startval, ...)
		
		scr <- .scr(prob, sc)
		scr1 <- .scr(prob, sc=1)
		
		##   Compute the TCC for the given startval
		tcc <- prob$p %*% scr
		
		##   Compute the derivative of the TCC for the given startval
		tcc.der <- prob$p1 %*% scr1
		
		new.ts <- startval-((truescore-tcc)/-tcc.der)
		
		return(new.ts)
	}

	
	##   This is a modified version of the {mixed} function
	.Mixed <- function(x, theta, ...) {
		
		##   Identify optional arguments that might be passed
		##   to the various functions used to compute response probabilities
		dots <- list(...)
		if (length(dots$D)) D <- dots$D else D <- 1
		if (length(dots$D.drm)) D.drm <- dots$D.drm else D.drm <- D
		if (length(dots$D.gpcm)) D.gpcm <- dots$D.gpcm else D.gpcm <- D
		if (length(dots$D.grm)) D.grm <- dots$D.grm else D.grm <- D
		if (length(dots$catprob)) catprob <- dots$catprob else catprob <- TRUE
		
		mod <- x@model
		p <- p1 <- NULL
		p.cat <- NULL
		p.mod <- NULL
		for (i in 1:length(mod)) {
			if (mod[i]=="drm") tmp <- .Drm(x, theta, D.drm, incorrect=TRUE)
			if (mod[i]=="gpcm") tmp <- .Gpcm(x, theta, D.gpcm)
			if (mod[i]=="grm") tmp <- .Grm(x, theta, catprob, D.grm)
			if (mod[i]=="nrm") tmp <- .Nrm(x, theta)
			if (mod[i]=="mcm") tmp <- .Mcm(x, theta)
			p <- cbind(p, tmp$p)
			p1 <- cbind(p1, tmp$p1)
			p.cat <- c(p.cat,tmp$cat)
			p.mod <- c(p.mod,tmp$mod)
		}
		return(list(p=p,p1=p1,p.cat=p.cat,p.mod=p.mod))
	}
	
	
	
	##   This is a modified version of the {drm} function
	.Drm <- function(x, theta, D.drm, incorrect) {
		
		items <- x@items$drm
		n <- length(items)
		a <- x@a[items,1]
		b <- x@b[items,1]
		c <- x@c[items,1]
		
		p <- p1 <- NULL
		for (i in 1:length(b)) {
			##   Compute the response probabilities
			cp <- c[i]+(1-c[i])/(1+exp(-D.drm*a[i]*(theta-b[i])))
			
			##   Compute the derivative of the response probabilities
			cp1 <- (D.drm*a[i]*(1-c[i])*exp(D.drm*a[i]*(theta-b[i])))/(1+exp(D.drm*a[i]*(theta-b[i])))^2
			
			if (incorrect==TRUE) {
				p <- cbind(p,(1-cp),cp) 
				p1 <- cbind(p1,(1-cp1),cp1) 
			} else {
				p <- cbind(p,cp)
				p1 <- cbind(p1,cp1)
			}
		}
		
		if (incorrect==TRUE) cat <- rep(2,n) else cat <- rep(1,n)
		mod <- rep("drm",n)
		return(list(p=p,p1=p1,cat=cat,mod=mod))
	}
	
	
	
	##   This is a modified version of the {gpcm} function
	.Gpcm <- function(x, theta, D.gpcm) {
		
		items <- x@items$gpcm
		n <- length(items)
		a <- x@a[items,1]
		b <- as.matrix(x@b[items,])
		if (length(items)==1) b <- t(b)
		cat <- x@cat[items]
		
		p <- p1 <- NULL 
		for (i in 1:nrow(b)) {
			dif <- 0 
			den <- den.der <- NULL 
			
			for (k in 0:(cat[i]-1)) {
				if (k>=1) dif <- dif+b[i,k]
				d <- exp(D.gpcm*a[i]*(k*theta-dif))
				d1 <- D.gpcm*a[i]*k*d
				
				den <- cbind(den, d)
				den.der <- cbind(den.der,d1)
			}
			den <- apply(den,1,sum)
			den.der <- apply(den.der,1,sum)
			
			dif <- 0
			for (k in 0:(cat[i]-1)) {
				if (k>=1) dif <- dif+b[i,k]
				
				##   Compute the category probabilities
				cp <- (exp(D.gpcm*a[i]*(k*theta-dif)))/den
				
				##   Compute the derivative of the category probabilities
				cp1 <- (den*D.gpcm*a[i]*k*(exp(D.gpcm*a[i]*(k*theta-dif)))-(exp(D.gpcm*a[i]*(k*theta-dif)))*den.der)/den^2
				
				p <- cbind(p,cp)
				p1 <- cbind(p1,cp1)
			}
		}
		mod <- rep("gpcm",n)
		return(list(p=p,p1=p1,cat=cat,mod=mod))
	}
	
	
	
	##   This is a modified version of the {grm} function
	.Grm <- function(x, theta, catprob, D.grm) {
		
		items <- x@items$grm
		n <- length(items)
		a <- x@a[items,1]
		b <- as.matrix(x@b[items,])
		if (length(items)==1) b <- t(b)
		cat <- x@cat[items]
		
		p <- p1 <- NULL 
		
		##   Compute category probabilities
		if (catprob==TRUE) { 
			for (i in 1:n) {
				ct <- cat[i]-1
				
				##   Compute probability for the lowest category
				cp <- 1-1/(1+exp(-D.grm*a[i]*(theta-b[i,1])))
				
				##   Compute the derivative of the probability for the lowest category
				cp1 <- -(D.grm*a[i]*exp(D.grm*a[i]*(theta-b[i,1])))/(1+exp(D.grm*a[i]*(theta-b[i,1])))^2
				
				p <- cbind(p, cp)
				p1 <- cbind(p1,cp1)
				
				for (k in 1:ct) {
					if (k<ct) {
						
						##   Compute the probability for the given category
						cp <- (1/(1+exp(-D.grm*a[i]*(theta-b[i,k]))))-(1/(1+exp(-D.grm*a[i]*(theta-b[i,k+1]))))
						
						##   Compute the derivative of the probability for the given category
						cpa <- (D.grm*a[i]*exp(D.grm*a[i]*(theta-b[i,k])))/(1+exp(D.grm*a[i]*(theta-b[i,k])))^2
						cpb <- (D.grm*a[i]*exp(D.grm*a[i]*(theta-b[i,k+1])))/(1+exp(D.grm*a[i]*(theta-b[i,k+1])))^2
						
						cp1 <- cpa-cpb
						
						
					} else if (k==ct) {
					
						##   Compute the probability for the highest category
						cp <- 1/(1+exp(-D.grm*a[i]*(theta-b[i,k])))
						
						##   Compute the derivative of the probability for the highest category
						cp1 <- (D.grm*a[i]*exp(D.grm*a[i]*(theta-b[i,k])))/(1+exp(D.grm*a[i]*(theta-b[i,k])))^2
						
					}
					
					p <- cbind(p, cp)
					p1 <- cbind(p1,cp1)
				}
			}
			
		# Compute cumulative probabilities
		} else if (catprob==FALSE) { 
			for (i in 1:n) {
				for (k in 1:(cat[i]-1)) {
					
					##   Compute the probability
					cp <- 1/(1+exp(-D.grm*a[i]*(theta-b[i,k])))
					
					##   Compute the derivative of the probability
					cp1 <- (D.grm*a[i]*exp(D.grm*a[i]*(theta-b[i,k])))/(1+exp(D.grm*a[i]*(theta-b[i,k])))^2
					
					p <- cbind(p, cp)
					p1 <- cbind(p1,cp1)
				}
			}
		}
		if (catprob==FALSE) cat <- cat-1
		mod <- rep("grm",n)
		return(list(p=p,p1=p1,cat=cat,mod=mod))
	}
	
	
	##   This is a modified version of the {nrm} function
	.Nrm <- function(x, theta) {
	
		items <- x@items$nrm
		n <- length(items)
		a <- as.matrix(x@a[items,]) 
		b <- as.matrix(x@b[items,])
		if (length(items)==1) {
			a <- t(a)
			b <- t(b)
		}
		cat <- x@cat[items]
		
		p <- p1 <- NULL 
		
		for (i in 1:n) {
			den <- den.der <- NULL
			a1 <- a[i,][!is.na(a[i,])]
			b1 <- b[i,][!is.na(b[i,])]
			for (k in 1:cat[i]) {
				d <- exp((theta*a1[k])+b1[k])
				d1 <- a1[k]*d
				
				den <- cbind(den, d)
				den.der <- cbind(den.der, d1)
			}
			
			den <- apply(den,1,sum)
			den.der <- apply(den.der,1,sum)
			
			for (k in 1:cat[i]) {
				##   Compute category probabilities
				cp <- exp((theta*a1[k])+b1[k])/den
				
				##   Compute derivative of category probabilities
				cp1 <- (den*a1[k]*exp((theta*a1[k])+b1[k])-exp((theta*a1[k])+b1[k])*den.der)/den^2
				
				p <- cbind(p,cp)
				p1 <- cbind(p1,cp1)
				
			}
		}
		mod <- rep("nrm",n)
		return(list(p=p,p1=p1,cat=cat,mod=mod))
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
		
		p <- p1 <- NULL 
		for (i in 1:n) {
			den <- den.der <- NULL
			a1 <- a[i,][!is.na(a[i,])]
			b1 <- b[i,][!is.na(b[i,])]
			c1 <- c[i,][!is.na(c[i,])]
			for (k in 1:cat[i]) {
				d <- exp((theta*a1[k])+b1[k])
				d1 <- a1[k]*d
				
				den <- cbind(den, d)
				den.der <- cbind(den.der, d1)
				
			}
			den <- apply(den,1,sum)
			den.der <- apply(den.der,1,sum)
			
			for (k in 2:cat[i]) {
			
				cp <- (exp((theta*a1[k])+b1[k])+c1[k-1]*(exp((theta*a1[1])+b1[1])))/den
				
				##   Compute derivative of category probabilities
				cp1 <- (den*a1[k]*exp((theta*a1[k])+b1[k])+c1[k-1]*a1[1]*(exp((theta*a1[1])+b1[1]))-exp((theta*a1[k])+b1[k])+c1[k-1]*(exp((theta*a1[1])+b1[1]))*den.der)/den^2
				
				p <- cbind(p,cp)
				p1 <- cbind(p1,cp1)
			}
		}
		mod <- rep("mcm",n)
		return(list(p=p,p1=p1,cat=cat-1,mod=mod))
	}
	
	
	##   This function is used to find the theta value that
	##   produces a TCC that is very close to the specified true score
	new.startval <- function(pars, ts, sc, ...) {
	
		##   Generate an initial set of theta values
		th <- seq(-6,6,.05)
		
		repeat {
		
			##   Create the TCCs
			prob <- .Mixed(pars, th,...)
			scr <- .scr(prob, sc)
			prob <- prob$p %*% scr
			
			##   Set the new start value to the theta that corresponds
			##   to the TCC value closest to the true score
			startval <-th[abs(ts-prob)==min(abs(ts-prob))]
			
			##   If the current range of theta values does not produce
			##   a TCC close enough to the true score, update the
			##   range of theta values
			if (min(abs(ts-prob))>1) {
			
				##   If the absolute value of theta is 20 logits, exit the function
				if (abs(startval)>=20) {
					startval <- -999
					break
				}
				
				##   Update the range of theta values
				if (startval<0) {
					th <- seq(startval-10,startval,0.05)
				} else {
					th <- seq(startval,startval+10,0.05)
				}
			} else {
				break
			}
		}
		
		return(startval)
	}
	
	
	##   This function is used to create the score weights for creating the TCCs
	.scr <- function(prob, sc) {
	
	##   prob - output from running .Mixed
	##   sc - {score} argument in the equate function
	
		##   Use the information in the {irt.prob} object tmp.to to 
		##   determine the number of columns in the matrix of probabilities
		##   associated with each item (i.e., the object p.cat)
		p.cat <- prob$p.cat
		
		##   Create a vector identifying the response model associated
		##   with each column of the matrix of probabilities
		p.mod <- rep(prob$p.mod,p.cat)
		
		##   Initialize the vector for the scoring function
		##   that will be used for the Stocking-Lord method
		scr <- NULL
		
		##   Create a set of default values for the scoring function
		##   that range from 1 to Kj for K columns of probabilities for item j
		for (h in 1:length(p.cat)) {
			scr <- c(scr,seq(1,p.cat[h]))
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
			} else if (sc>2) {
				warning("The value specified for {score} is invalid. Score was set to 1")
				scr <- scr-1
				cat <- rep(p.cat,p.cat)
				scr[cat==1] <- 1
			}
			
		##   Use this if the researcher supplies a vector of score 
		##   weights for all of the columns in p
		} else {
			if (length(sc)==length(scr)) {
				scr <- sc
			} else {
				warning("The length of {score} does not match the number of response probabilities. Score was set to 1")
			}
		}
		
		return(scr)
	}
	
	
	##########   START EQUATE  ##########
	
	##   Number of groups
	ng <- x@groups
	
	##   Make {score} a list with length equal to the number of groups
	if(!is.list(score)) {
		tmp <- vector("list",ng)
		for (i in 1:ng) {
			tmp[[i]] <- score
		}
		score <- tmp
	}
	
	##   TRUE SCORE EQUATING
	if ("TSE" %in% method) {
		
		##   Separate the item parameters for each of the groups
		pars <- sep.pars(x)
		
		if (ng==1) pars <- list(pars)
		
		##   Set a starting value for the optimization (if necessary)
		if (missing(startval)) startval <- -4
		
		##   Initialize object to hold all of the equated true scores
		tse.out <- NULL
		
		##   Generate a set of observed scores that fall within 
		##   the range of true scores
		ts.low.flag <- FALSE
		if (missing(true.scores)) {
			tmp.min <- sum(pars[[base.grp]]@c, na.rm=TRUE)
			tmp.max <- sum(pars[[base.grp]]@cat-1)
			true.scores <- ceiling(tmp.min):(tmp.max)
			ts.low.flag <- TRUE
		}
		
		##   Initialize an object to indentify the lowest true score
		##   value with a corresponding true score for another test
		start <- 0
		
		##   Initialize an object to indicate that the remaining set of true
		##   scores in {ts}have no corresponding theta values (in a reasonable range)
		na.flag <- FALSE
		
		##   Loop through all true scores in the base group
		for (ts in true.scores) {
			
			##   Initialize an object that will be updated if the convergence
			##   goes awry and an attempt is made to get a value closer
			##   to the final theta value (this type of attempt will only
			##   be made twice)
			flag <- 0
			
			if (na.flag==TRUE) {
				##   All of the remaining true scores in {ts}will
				##   have theta values equal to NA
				tse.out <- rbind(tse.out, c(NA,ts))
				next
			}
			
			##   Initialize the number of iterations
			iter <- 1
			
			##   Use the Newton-Raphson method to find the
			repeat{
				theta <- .TSE(startval, ts, pars[[base.grp]], score[[base.grp]], ...)
				
				##   Check to see if the convergence has gone awry
				if (theta==Inf|theta==-Inf|is.nan(theta)) {
				
					##   Attempt to get a starting value closer to the final theta value
					startval <- new.startval(pars[[base.grp]], ts, score[[base.grp]], ...)
					flag <- flag+1
					
					##   When no theta value can be found for the specified true score
					if (startval==-999|flag==2) {
						##   Set the associated theta value to missing
						theta <- NA
						
						##   Check to see if some other theta value has already been found
						##   If so, make it so all remaining theta values will be missing
						##   (this assumes that the remaining true scores are higher than
						##   the current true score)
						if (start>0) na.flag <- TRUE
						
						##   Reset the starting value
						startval <- -6
						
						##   Move to the next true score
						break
					}
					
				} else {
					##   Check to see if the new estimate is sufficiently
					##   close to the old estimate of theta
					if (abs(startval-theta)<1e-10) {
					
						##   Update the start object
						if (start==0) start <- ts-1
						
						##   Move to the next true score
						break 
						
					##   If the convergence criterion is not met
					} else {
						##   Update the iteration number and the 
						##   current value of theta
						iter <- iter+1
						startval <- theta
					}
					
					##   Check to see if the maximum number of iterations has been reached
					if (iter==100) {
						cat(paste("Maximum iterations reached for true score:",ts,"\n"))
						
						##   Move to the next true score
						break
					}
				}
			}
			
			##   Compile the matrix of equated true scores
			tse.out <- rbind(tse.out, c(theta,ts))
		}
		
		if (ng>1) {
			##   Compute equated true scores for all other groups
			for (i in c(1:ng)[-base.grp]) {
			
				##   Compute response probabilities for 
				prob <- .Mixed(pars[[i]], tse.out[,1], ...)
				
				scr <- .scr(prob, score[[i]])
				
				##   Compute the TCC for the given startval
				tcc <- prob$p %*% scr
				
				tse.out <- cbind(tse.out, round(tcc,6))
			}
			colnames(tse.out) <- c("theta",names(x@pars)[base.grp],names(x@pars)[-base.grp])
		} else {
			colnames(tse.out) <- c("theta","x")
		}
		
		##   Use Kolen's (1981) approach to interpolate values
		##   for equated true scores in the range of observed scores
		##   (from one to the value below the lowest estimated true score)
		if (ng>1) {
		
			if (ts.low==TRUE & ts.low.flag==TRUE) {
				
				##   Compute the sum of the guessing parameters for the base group
				tmp <- sum(pars[[base.grp]]@c, na.rm=TRUE)
				
				##   Initialize an object to store the sum of the guessing
				##   parameters for each group
				sum.c <- NULL
				
				##   Compute the sum of the guessing parameters for the given group
				for (i in c(1:ng)[-base.grp]) {
					sum.c <- c(sum.c, sum(pars[[i]]@c, na.rm=TRUE)/tmp)
				}
				sum.c[is.nan(sum.c)|sum.c==Inf|sum.c==-Inf] <- 1
				
				##   Initialize a vector of observed scores below the true score range
				tse.low <- rep(0,ng)
				
				for (i in 1:start) {
					tse.low <- rbind(tse.low, c(i,sum.c*i))
				}
				tse.low <- cbind(NA,tse.low)
				row.names(tse.low) <- NULL
				colnames(tse.low) <- colnames(tse.out)
				
				
				##   Check to see if there are any rows with NAs for theta after {start}
				tmp <- tse.out[tse.out[,2]>start,]
				if (nrow(tmp[is.na(tmp[,1]),]) >0){
					cat(paste("The maximum possible score is ",max(tse.out[,2])+1,", but theta equivalents cannot be computed for true scores greater than ",
					max(tse.out[!is.na(tse.out[,1]),2]),"\n",sep=""))
				}
				 
				tse.out <- tse.out[!is.na(tse.out[,1]),]
				
				##  Combine the lower estimated true scores with equated true scores
				tse.out <- rbind(tse.low, tse.out)
				tse.out <- tse.out[!is.nan(tse.out[,3]),]
				
				##  There is a possibility that duplicate values are created. Remove them
				tmp <- table(tse.out[,2])
				tmp <- as.numeric(names(tmp[tmp>1]))
				if (length(tmp)) {
					tmp <- c(1:nrow(tse.out))[tse.out[,2] %in% tmp & is.na(tse.out[,1])]
					tse.out <- tse.out[-tmp,]
				}
			}
		}
		
	} 
	
	##   OBSERVED SCORE EQUATING
	if ("OSE" %in% method) {
	
		##   Separate the item parameters for each of the groups
		pars <- sep.pars(x)
		
		##   Identify the group names
		if (ng==1) {
			pars <- list(pars)
			nms <- "x"
		} else {
			nms <- names(x@pars)
		}
		
		##   Initialize an object to store the compound binomial/multinomial
		##   distributions for each group
		dist <- dist1 <- dist2 <- dist1b <- dist2b <- vector("list",ng-1)
		
		##   Loop through all of the groups
		for (i in 1:(ng-1)) {
		
			if (i>=base.grp) grp <- i+1 else grp <- i
			
			
			##   Identify the weights to be used for population 1
			if (missing(weights1)) {
				wgt1 <- as.weight(normal.wt=TRUE)
			} else {
				if (is.list(weights1[[1]])) wgt1 <- weights1[[i]] else wgt1 <- weights1
			}
			
			##   Identify the weights to be used for population 2
			if (missing(weights2)) {
				wgt2 <- wgt1
			} else {
				if (is.list(weights2[[1]])) wgt2 <- weights2[[i]] else wgt2 <- weights2
			}
			
			##   Extract the theta values that will be used for each population
			##   when computing the probabilities
			theta1 <- wgt1[[1]]
			theta2 <- wgt2[[1]]
			
			##   Identify the synthetic weights to be used
			if (missing(syn.weights)) {
				syn <- c(.5,.5)
			} else {
				if (is.list(syn.weights)) syn <- syn.weights[[i]] else syn <- syn.weights
			}
			
			##   Compute response probabilities
			prob.b1 <- .Mixed(pars[[base.grp]], theta1, incorrect=TRUE, ...)
			prob.b2 <- .Mixed(pars[[base.grp]], theta2, incorrect=TRUE, ...)
			
			prob1 <- .Mixed(pars[[grp]], theta1, incorrect=TRUE, ...)
			prob2 <- .Mixed(pars[[grp]], theta2, incorrect=TRUE, ...)
			
			##   Extract the probabilities
			pb1 <- prob.b1$p
			pb2 <- prob.b2$p
			
			p1 <- prob1$p
			p2 <- prob2$p
			
			##   Identify the number of columns of probabilities associated with each item
			pb.cat <- prob.b1$p.cat
			p.cat <- prob1$p.cat
			
			##   Identify the number of theta values used
			##   to compute the probabilities
			n1 <- nrow(p1)
			n2 <- nrow(p2)
			
			##   Maximum possible score
			 cat.b <- pars[[base.grp]]@cat
			 cat <- pars[[grp]]@cat
			
			##   Correct the maximum possible score if there are MCM items
			if ("mcm" %in% pars[[base.grp]]@model) cat.b[pars[[base.grp]]@items$mcm] <- cat.b[pars[[base.grp]]@items$mcm]-1
			if ("mcm" %in% pars[[grp]]@model) cat.b[pars[[grp]]@items$mcm] <- cat.b[pars[[grp]]@items$mcm]-1
			
			
			##   Initialize lists to store the distributions for each
			##   observed score for each group
			dist1b[[i]] <- vector("list",2)
			dist2b[[i]] <- vector("list",2)
			dist1[[i]] <- vector("list",2)
			dist2[[i]] <- vector("list",2)
			
			##   Create the compound binomial/multinomial distributions for the base group
			for (j in 1:length(cat.b)) {
			
				if (j==1) {
					dist1b[[i]][[1]] <- pb1[,1:pb.cat[j]]
					dist2b[[i]][[1]] <- pb2[,1:pb.cat[j]]
					
					##  Minimum score for all further observed scores
					min.scr <- 2
					
				} else {
					
					cols <- (sum(pb.cat[1:(j-1)])+1):sum(pb.cat[1:j])
					
					##   Identify the maximum score for the given set of items
					max.scr <- max(cols)
					
					##   Initialize an object for the distribution for the given set of items
					dist1b[[i]][[2]] <- matrix(0,n1,max.scr-min.scr+1)
					dist2b[[i]][[2]] <- matrix(0,n2,max.scr-min.scr+1)
					
					##   Create the distributions
					for (k in 1:pb.cat[j]) {
						tmp <- dist1b[[i]][[1]]*pb1[,cols[k]]
						dist1b[[i]][[2]][,k:(ncol(tmp)+k-1)] <- dist1b[[i]][[2]][,k:(ncol(tmp)+k-1)]+tmp
						
						tmp <- dist2b[[i]][[1]]*pb2[,cols[k]]
						dist2b[[i]][[2]][,k:(ncol(tmp)+k-1)] <- dist2b[[i]][[2]][,k:(ncol(tmp)+k-1)]+tmp
						
					}
					
					##   Update the minimum score
					min.scr <- min.scr+1
					
					##   Replace the distribution in the first list
					##   element with the current distribution
					dist1b[[i]][[1]] <- dist1b[[i]][[2]]
					dist2b[[i]][[1]] <- dist2b[[i]][[2]]
				}
			}
			
			##   Create the compound binomial/multinomial distributions for the focal group
			for (j in 1:length(cat)) {
			
				if (j==1) {
					dist1[[i]][[1]] <- p1[,1:p.cat[j]]
					dist2[[i]][[1]] <- p2[,1:p.cat[j]]
					
					##  Minimum score for all further observed scores
					min.scr <- 2
					
				} else {
					
					cols <- (sum(p.cat[1:(j-1)])+1):sum(p.cat[1:j])
					
					##   Identify the maximum score for the given set of items
					max.scr <- max(cols)
					
					##   Initialize an object for the distribution for the given set of items
					dist1[[i]][[2]] <- matrix(0,n1,max.scr-min.scr+1)
					dist2[[i]][[2]] <- matrix(0,n2,max.scr-min.scr+1)
					
					##   Create the distributions
					for (k in 1:p.cat[j]) {
						tmp <- dist1[[i]][[1]]*p1[,cols[k]]
						dist1[[i]][[2]][,k:(ncol(tmp)+k-1)] <- dist1[[i]][[2]][,k:(ncol(tmp)+k-1)]+tmp
						
						tmp <- dist2[[i]][[1]]*p2[,cols[k]]
						dist2[[i]][[2]][,k:(ncol(tmp)+k-1)] <- dist2[[i]][[2]][,k:(ncol(tmp)+k-1)]+tmp
						
					}
					
					##   Update the minimum score
					min.scr <- min.scr+1
					
					##   Replace the distribution in the first list
					##   element with the current distribution
					dist1[[i]][[1]] <- dist1[[i]][[2]]
					dist2[[i]][[1]] <- dist2[[i]][[2]]
				}
			}
			
			##   Extract the last list element
			##   This is the final distribution
			dist1b[[i]] <- t(dist1b[[i]][[2]]) %*% wgt1[[2]]
			dist2b[[i]] <- t(dist2b[[i]][[2]]) %*% wgt2[[2]]
			dist1[[i]] <- t(dist1[[i]][[2]]) %*% wgt1[[2]]
			dist2[[i]] <- t(dist2[[i]][[2]]) %*% wgt2[[2]]
			
			##   Create the synthetic distributions
			s1 <- dist1b[[i]]*syn[1]+dist2b[[i]]*syn[2]
			s2 <- dist1[[i]]*syn[1]+dist2[[i]]*syn[2]
			
			##   Compile the distributions to be output
			dist[[i]] <- list(cbind(0:(nrow(dist1b[[i]])-1),dist1b[[i]],dist2b[[i]],s1),cbind(0:(nrow(dist1[[i]])-1),dist1[[i]],dist2[[i]],s2))
			names(dist[[i]]) <- c(nms[base.grp],nms[grp])
			colnames(dist[[i]][[1]]) <- colnames(dist[[i]][[2]]) <- c("score","pop1","pop2","syn")
			
			##   Compute cumulative proportions
			F1 <- cumsum(s1)
			F2 <- cumsum(s2)
			
			##   Compute the percentile ranks for the base group
			for (j in 1:length(F1)) {
				if (j==1) {
					pr <- 0
				} else {
					pr <- c(pr, F1[j-1]+((j-1)-((j-1)-.5))*(F1[j]-F1[j-1]))
				}
			}
			
			bd <- dist[[i]][[1]]
			
			##   Initialize an object to store the equated scores
			if (i==1) OSE <- bd[,1]
			
			##   Conduct the equipercentile equating
			tmp.ose <- NULL
			for (j in 1:length(F1)) {
					tmp <- 1:length(F2)
					tmp <- min(tmp[F2>pr[j]])
				if (tmp==1) {
					tmp.ose <- c(tmp.ose, (pr[j]/F2[tmp])+(tmp-1.5))
				} else {
					tmp.ose <- c(tmp.ose, ((pr[j]-F2[tmp-1])/(F2[tmp]-F2[tmp-1]))+(tmp-1.5))
				}
			}
			OSE <- cbind(OSE,tmp.ose)
			colnames(OSE)[1] <- nms[base.grp]
			colnames(OSE)[ncol(OSE)] <- nms[grp]
		}
		
		##   Adjust the equated zero observed score for all groups
		OSE[1,] <- 0
		
		if (length(dist)==1) dist <- dist[[1]]
		ose.out <- list(scores=OSE,dist=dist)
	}
	
	##   Return the equated true scores and/or observed scores
	if ("TSE" %in% method) {
		if ("OSE" %in% method) {
			return(list(tse=tse.out,ose=ose.out))
		} else {
			return(tse.out)
		}
	} else {
		return(ose.out)
	}
	
})