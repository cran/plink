##   This function computes response probabilities for items
##   modeled using the graded response model and the 
##   multidimensional graded response model

setGeneric("grm", function(x, cat, theta, dimensions=1, catprob=FALSE, D=1, location=FALSE, items, information=FALSE, angle, ...) standardGeneric("grm"))


setMethod("grm", signature(x="matrix", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, items, information, angle, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("grm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, items, information, angle, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



setMethod("grm", signature(x="list", cat="numeric"), function(x, cat, theta, dimensions, catprob, D, location, items, information, angle, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])),"grm")
	x <- sep.pars(x, cat, poly.mod, dimensions, location, loc.out=FALSE, ...)
	callGeneric()
	
})



##   For this method the objects, cat, dimensions, and location are contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("grm", signature(x="irt.pars", cat="ANY"), function(x, cat, theta, dimensions, catprob, D, location, items, information, angle, ...) {
	
	##   Loop through all groups. In this scenario, a list of {irt.prob} objects will be returned
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], x@dimensions[i], x@location[i], loc.out=FALSE, ...)
			out[[i]] <- grm(tmp, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, x@dimensions, x@location, loc.out=FALSE, ...)
		callGeneric()
	}
	
})



##   For this method the objects, cat, dimensions, and location are contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("grm", signature(x="sep.pars", cat="ANY"), function(x, cat, theta, dimensions, catprob, D, location, items, information, angle, ...) {

	##   The equation to compute probabilities is not (actually) parameterized using
	##   the location/threshold-deviation formulation. As such, in instances where a location
	##   parameter is included, it is necessary to reformat the parameters appropriately
	if (x@location==TRUE) {
		pars <- list(x@a,x@b,x@c)
		pm <- as.poly.mod(x@n[1], x@model, x@items)
		x <- sep.pars(pars, x@cat, pm, x@dimensions, location=TRUE, loc.out=FALSE)
	}
	
	##   Number of dimensions
	dimensions <- x@dimensions
	
	##   Identify the grm items
	if (missing(items)) items <- 1:x@n[1]
	tmp.items <- x@items$grm
	items <- tmp.items[tmp.items%in%items]
	
	##   Number of items
	n <- length(items)
	
	##   Extract the grm items and rescale the slope parameters (if necessary)
	a <- as.matrix(x@a[items,1:dimensions])
	b <- as.matrix(x@b[items,])
	
	##   If there is only a single item, the matrices specified above will have
	##   the wrong orientation. For example, the threshold parameters for this item
	##   will be in different rows of the matrix b instead of being in a matrix
	##   with a single row and multiple columns. As such, these matrices need
	##   to be transposed
	if (n==1) {
		a <- t(a)
		b <- t(b)
	}
	pars <- list(a=a, b=b, c=x@c[items,])
	
	cat <- x@cat[items]
	
	##   Check to see if the argument {D.grm} was passed via the function {mixed}
	dots <- list(...)
	if (length(dots$D.grm)) D <- dots$D.grm
	
	##   Generate theta values if {theta} is missing
	##   Different values should be generated depending on the number of dimensions
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
		##   If the user (purposefully or accidentally) specifies {theta} as a matrix
		##   or a list instead of a vector for the unidimensional case, turn all of the 
		##   values into a vector
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
		##   If, in the multidimensional case, only a vector of theta values is 
		##   supplied, treat this as a vector for each dimension then create all
		##   permutations of these values. If {theta} is formatted as a matrix
		##   or list from the outset, just find the permutations
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
	
	if (length(x@model[x@model!="grm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the grm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
	##   Initialize object to hold the response probabilities
	p <- p1 <- NULL
	
	##   Determine angles for computing information (in the multidimensional case)
	if (information==TRUE) {
		if (dimensions>1) {
			if (missing(angle)) {
				angle <- list()
				for (i in 1:(dimensions-1)) {
					angle[[i]] <- seq(0,90,10)
				}
				ang <- expand.grid(angle)
				angle <- as.matrix(cbind(ang[,1],90-ang[,1],ang[,-1]))
			} else {
				if (is.vector(angle)) {
					angle1 <- angle
					angle <- list()
					for (i in 1:(dimensions-1)) {
						angle[[i]] <- angle1
					}
					ang <- expand.grid(angle)
					angle <- as.matrix(cbind(ang[,1],90-ang[,1],ang[,-1]))
				} else if (is.matrix(angle)) {
					if (ncol(angle)!=dimensions) {
						warning("The number of columns in {angle} does not match the number of dimensions in {x}. Default angles were used.")
						angle <- list()
						for (i in 1:(dimensions-1)) {
							angle[[i]] <- seq(0,90,10)
						}
						ang <- expand.grid(angle)
						angle <- as.matrix(cbind(ang[,1],90-ang[,1],ang[,-1]))
					}
				}
			}
			dcos <- cos((pi*angle)/180)
		}
	}
	
	##   Compute category probabilities
	if (catprob==TRUE) { 
		for (i in 1:n) {
			ct <- cat[i]-1
			info <- rep(0,nrow(theta))
			
			##   Compute the probabilities (and derivative if applicable) for the lowest category
			if (dimensions==1) {
				cp <- 1-1/(1+exp(-D*a[i]*(theta-b[i,1])))
				
				if (information==TRUE) {
					cp1 <- -(D*a[i]*exp(D*a[i]*(theta-b[i,1])))/(1+exp(D*a[i]*(theta-b[i,1])))^2
					cp2 <- -(D*a[i]^2*exp(D*a[i]*(theta-b[i,1])))/(1+exp(D*a[i]*(theta-b[i,1])))^2
					tmp.info <- ((cp1^2)/cp)-cp2
					
				}
				
			} else {
				##   In the multidimensional case D is typically set equal to 1
				a[i,] <- a[i,]*D
				
				cp <- 1-1/(1+exp(-(theta %*% a[i,]+b[i,1])))
				
				if (information==TRUE) {
					cp1 <- (-exp(theta %*% a[i,]+b[i,1])/(1+exp(theta %*% a[i,]+b[i,1]))^2)%*%(a[i,]%*%t(dcos))
					cp2 <- (-exp(theta %*% a[i,]+b[i,1])/(1+exp(theta %*% a[i,]+b[i,1]))^2)%*%(a[i,]%*%t(dcos))
					tmp.info <- (as.matrix(cp1^2)/matrix(cp,length(cp),nrow(angle)))-as.matrix(cp2)
				}
			}
			
			p <- cbind(p, cp)
			colnames(p)[ncol(p)] <- paste("item_",items[i],".0",sep="")
			if (information==TRUE) info <- info+tmp.info
			
			
			for (k in 1:ct) {
				if (dimensions==1) {
					if (k<ct) {
						cp <- (1/(1+exp(-D*a[i]*(theta-b[i,k]))))-(1/(1+exp(-D*a[i]*(theta-b[i,k+1]))))
						
						if (information==TRUE) {
							cpa <- (D*a[i]*exp(D*a[i]*(theta-b[i,k])))/(1+exp(D*a[i]*(theta-b[i,k])))^2
							cpb <- (D*a[i]*exp(D*a[i]*(theta-b[i,k+1])))/(1+exp(D*a[i]*(theta-b[i,k+1])))^2
							cp1 <- cpa-cpb
							
							cpa1 <- (D*a[i]^2*exp(D*a[i]*(theta-b[i,k])))/(1+exp(D*a[i]*(theta-b[i,k])))^2
							cpb1 <- (D*a[i]^2*exp(D*a[i]*(theta-b[i,k+1])))/(1+exp(D*a[i]*(theta-b[i,k+1])))^2
							cp2 <- cpa1-cpb1
							tmp.info <- ((cp1^2)/cp)-cp2
						}
					} else if (k==ct) {
						##   Compute the probabilities for the highest category
						cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
						
						if (information==TRUE) {
							cp1 <- (D*a[i]*exp(D*a[i]*(theta-b[i,k])))/(1+exp(D*a[i]*(theta-b[i,k])))^2
							cp2 <- (D*a[i]^2*exp(D*a[i]*(theta-b[i,k])))/(1+exp(D*a[i]*(theta-b[i,k])))^2
							tmp.info <- ((cp1^2)/cp)-cp2
						}
					}
				} else {
					if (k<ct) {
						cp <- (1/(1+exp(-(theta %*% a[i,]+b[i,k]))))-(1/(1+exp(-(theta %*% a[i,]+b[i,k+1]))))
						
						if (information==TRUE) {
							cpa <- (exp(theta %*% a[i,]+b[i,k])/(1+exp(theta %*% a[i,]+b[i,k]))^2)%*%(a[i,]%*%t(dcos))
							cpb <- (exp(theta %*% a[i,]+b[i,k+1])/(1+exp(theta %*% a[i,]+b[i,k+1]))^2)%*%(a[i,]%*%t(dcos))
							cp1 <- cpa-cpb
							
							cpa1 <- (exp(theta %*% a[i,]+b[i,k])/(1+exp(theta %*% a[i,]+b[i,k]))^2)%*%(a[i,]^2%*%t(dcos))
							cpb1 <- (exp(theta %*% a[i,]+b[i,k+1])/(1+exp(theta %*% a[i,]+b[i,k+1]))^2)%*%(a[i,]^2%*%t(dcos))
							cp2 <- cpa1-cpb1
							tmp.info <- (as.matrix(cp1^2)/matrix(cp,length(cp),nrow(angle)))-as.matrix(cp2)
						}
						
					} else if (k==ct) {
						##   Compute the probabilities for the highest category
						cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
						
						if (information==TRUE) {
							cp1 <- (exp(theta %*% a[i,]+b[i,k])/(1+exp(theta %*% a[i,]+b[i,k]))^2)%*%(a[i,]%*%t(dcos))
							cp2 <- (exp(theta %*% a[i,]+b[i,k])/(1+exp(theta %*% a[i,]+b[i,k]))^2)%*%(a[i,]^2%*%t(dcos))
							tmp.info <- (as.matrix(cp1^2)/matrix(cp,length(cp),nrow(angle)))-as.matrix(cp2)
						}
					}
				}
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",items[i],".",k,sep="")
				if (information==TRUE) info <- info+tmp.info
			}
			if (information==TRUE) p1 <- cbind(p1, as.vector(info))
		}
		
	# Compute cumulative probabilities
	} else if (catprob==FALSE) { 
		for (i in 1:n) {
			info <- rep(0,nrow(theta))
			
			if (dimensions>1) a[i,] <- a[i,]*D
			
			for (k in 1:(cat[i]-1)) {
				if (dimensions==1) {
					cp <- 1/(1+exp(-D*a[i]*(theta-b[i,k])))
					if (information==TRUE) {
						cp1 <- (D*a[i]*exp(D*a[i]*(theta-b[i,k])))/(1+exp(D*a[i]*(theta-b[i,k])))^2
						tmp.info <- (cp1^2)/(cp*(1-cp))
					}
				} else {
					cp <- 1/(1+exp(-(theta %*% a[i,]+b[i,k])))
					if (information==TRUE) {
						cp1 <- (exp(theta %*% a[i,]+b[i,k])/(1+exp(theta %*% a[i,]+b[i,k]))^2)%*%(a[i,]^2%*%t(dcos))
						tmp.info <- as.matrix(cp1^2)/matrix(cp*(1-cp),length(cp),nrow(angle))
					}
					
				}
				p <- cbind(p, cp)
				colnames(p)[ncol(p)] <- paste("item_",items[i],".",k,sep="")
				if (information==TRUE) info <- info+tmp.info
			}
			if (information==TRUE) p1 <- cbind(p1, as.vector(info))
		}
	}
	
	p <- data.frame(cbind(theta,p))
	if (catprob==FALSE) cat <- cat-1
	
	##   Create and return the irt.prob object
	if (information==TRUE) {
		if (dimensions>1) {
			th <- NULL
			for (i in 1:nrow(angle)) {
				th <- rbind(th, cbind(theta,matrix(angle[i,],nrow(theta),dimensions,byrow=TRUE)))
			}
			colnames(th) <- c(paste("theta",1:dimensions,sep=""),paste("angle",1:dimensions,sep=""))
			p1 <- data.frame(cbind(th,p1))
			names(p1)[-c(1:(2*dimensions))] <- paste("item_",items,sep="")
		} else {
			p1 <- data.frame(cbind(theta,p1))
			names(p1) <- c("theta",paste("item_",items,sep=""))
		}
		p <- new("irt.prob", prob=p, info=p1, p.cat=cat, mod.lab=x@mod.lab[x@model=="grm"], dimensions=dimensions, D=c(D.grm=D), pars=pars, model="grm", items=list(grm=1:n))
	} else {
		p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="grm"], dimensions=dimensions, D=c(D.grm=D), pars=pars, model="grm", items=list(grm=1:n))
	}
	
	return(p)
	
})