##   This function computes response probabilities for items
##   modeled using the nominal response model and the
##   multidimensional nominal response model

setGeneric("nrm", function(x, cat, theta, dimensions=1, items, information=FALSE, angle, ...) standardGeneric("nrm"))



setMethod("nrm", signature(x="matrix", cat="numeric"), function(x, cat, theta, dimensions, items, information, angle, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"nrm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
	
})



setMethod("nrm", signature(x="data.frame", cat="numeric"), function(x, cat, theta, dimensions, items, information, angle, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(x),"nrm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
})


setMethod("nrm", signature(x="list", cat="numeric"), function(x, cat, theta, dimensions, items, information, angle, ...) {
	
	if(!hasArg(poly.mod)) poly.mod <- as.poly.mod(nrow(as.matrix(x[[1]])),"nrm")
	x <- sep.pars(x, cat, poly.mod, dimensions, ...)
	callGeneric()
	
})



##   For this method the objects, cat and dimensions are contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("nrm", signature(x="irt.pars", cat="ANY"), function(x, cat, theta, dimensions, items, information, angle, ...) {
	
	##   Loop through all groups. In this scenario, a list of {irt.prob} objects will be returned
	if (x@groups>1) {
		out <- vector("list", x@groups)
		for (i in 1:x@groups) {
			tmp <- sep.pars(x@pars[[i]], x@cat[[i]], x@poly.mod[[i]], dimensions=x@dimensions[i], ...)
			out[[i]] <- nrm(tmp, ...)
		}
		names(out) <- names(x@pars)
		return(out)
	} else {
		x <- sep.pars(x@pars, x@cat, x@poly.mod, dimensions=x@dimensions, ...)
		callGeneric()
	}
	
})



##   For this method the objects, cat and dimensionsare contained in {x} 
##   As such, these arguments are treated as missing in the signature
setMethod("nrm", signature(x="sep.pars", cat="ANY"), function(x, cat, theta, dimensions, items, information, angle, ...) {

	##   Identify the nrm items
	if (missing(items)) items <- 1:x@n[1]
	tmp.items <- x@items$nrm
	items <- tmp.items[tmp.items%in%items]
	
	##   Number of items
	n <- length(items)
	
	dimensions <- x@dimensions
	
	##   Extract the mcm items
	##   When creating the object with the item parameters
	##   they should be grouped first by dimensions then by categories. 
	##   For example, for an item with 2 dimensions and 4 categories
	##   we would have (a11,a12,a13,a14,a21,a22,a23,a24)
	##   This grouping applies for both the slope parameters (named a)
	##   and the category parameters (named b)
	a <- x@a[items,] 
	b <- x@b[items,]
	
	##   If there is only a single item, the matrices specified above will have
	##   the wrong orientation. For example, the slope parameters for this item
	##   will be in different rows of the matrix a instead of being in a matrix
	##   with a single row and multiple columns. As such, these matrices need
	##   to be transposed
	if (n==1) {
		a <- t(a)
		b <- t(b)
	}
	
	pars <- list(a=a, b=b, c=x@c[items,])
	cat <- x@cat[items]
	
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
	
	if (length(x@model[x@model!="nrm"])) warning("{x} contains mixed format items. Probabilities will only be computed for the nrm polytomous items.\nTo compute probabilities for mixed format items, use the function {mixed}.\n")
	
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
	
	
	for (i in 1:nrow(b)) {
		##   Object for the denominator in the final NRM equation
		den <- NULL
		
		info <- rep(0,nrow(theta))
		
		##   Because of how the parameters are organized in {x}
		##   there may be NAs in various rows. Remove these NAs
		##   before computing the response probabilities
		a1 <- a[i,][!is.na(a[i,])]
		b1 <- b[i,][!is.na(b[i,])]
		
		##   Compute the denominator
		for (k in 1:cat[i]) {
			k1 <- seq(k,cat[i]*dimensions,cat[i])
			d <- exp((theta %*% a1[k1])+b1[k])
			den <- cbind(den, d)
		}
		
		den <- apply(den,1,sum)
		
		##   Compute the response probabilities
		for (k in 1:cat[i]) {
			k1 <- seq(k,cat[i]*dimensions,cat[i])
			cp <- exp((theta %*% a1[k1])+b1[k])/den
			
			if (information==TRUE) {
				tmp.p1 <- tmp.p2 <- tmp.p3 <- rep(0,nrow(theta))
				for (v in 1:cat[i]) {
					if (dimensions==1) {
						tmp.p1 <- tmp.p1+exp((theta*a1[v])+b1[v])*a1[v]
						tmp.p2 <- tmp.p2+exp((theta*a1[v])+b1[v])*(a1[k]-a1[v])
						tmp.p3 <- tmp.p3+exp((theta*a1[v])+b1[v])*(a1[k]^2-a1[v]^2)
					} else {
						v1 <- seq(v,cat[i]*dimensions,cat[i])
						tmp.p1 <- tmp.p1+exp((theta%*%a1[v1])+b1[v])%*%(a1[v1]%*%t(dcos))
						tmp.p2 <- tmp.p2+exp((theta%*%a1[v1])+b1[v])%*%(a1[k1]%*%t(dcos)-a1[v1]%*%t(dcos))
						tmp.p3 <- tmp.p3+exp((theta%*%a1[v1])+b1[v])%*%((a1[k1]%*%t(dcos))^2-(a1[v1]%*%t(dcos))^2)
					}
				}
				
				cp1 <- (as.vector(exp((theta%*%a1[k1])+b1[k]))*tmp.p2)/den^2
				cp2 <- ((den^2)*as.vector(exp((theta%*%a1[k1])+b1[k]))*tmp.p3-as.vector(exp((theta%*%a1[k1])+b1[k]))*tmp.p2*2*den*tmp.p1)/den^4
				tmp.info <- ((cp1^2)/as.vector(cp))-cp2
			}
			p <- cbind(p,cp)
			colnames(p)[ncol(p)] <- paste("item_",items[i],".",k,sep="")
			if (information==TRUE) info <- info+tmp.info
		}
		if (information==TRUE) p1 <- cbind(p1, as.vector(info))
	}
	p <- data.frame(cbind(theta,p))
	
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
		p <- new("irt.prob", prob=p, info=p1, p.cat=cat, mod.lab=x@mod.lab[x@model=="nrm"], dimensions=dimensions, D=c(D=1), pars=pars, model="nrm", items=list(nrm=1:n))
	} else {
		p <- new("irt.prob", prob=p, p.cat=cat, mod.lab=x@mod.lab[x@model=="nrm"], dimensions=dimensions, D=c(D=1), pars=pars, model="nrm", items=list(nrm=1:n))
	}
	return(p)
})
