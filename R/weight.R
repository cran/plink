as.weight <- function(n, theta, weight, quadrature=FALSE, normal.wt=FALSE, dimensions=1, ...) {
	if (dimensions<1) stop("You must specify at least one dimension")
	# Update/set additional arguments
	dots <- list(...)
	
	if (missing(theta)) {
		if (missing(weight)) {
			if (missing(n)) {
				if (dimensions==1) n <- 40 else n <- rep(7,dimensions)
			} else {
				if(length(n)!=dimensions) {
					if (length(n)==1) {
						n <- rep(n,dimensions)
					} else {
						stop("The number of elements in {n} must match the number of dimensions")
					}
				}
			}
			p <- w <- vector("list",length(n))
			if (quadrature==TRUE) {
				if (length(dots$dist)) dist <- dots$dist else dist <- "normal"
				if (length(dots$l)) l <- rep(dots$l,length.out=length(n)) else l <- rep(0,length(n))
				if (length(dots$u)) u <- rep(dots$u,length.out=length(n)) else u <- rep(1,length(n))
				
				if (length(dots$mu)) mu <- rep(dots$mu,length.out=length(n)) else mu <- rep(0,length(n))
				if (length(dots$sigma)) sigma <- rep(dots$sigma,length.out=length(n)) else sigma <- rep(1,length(n))
				
				if (length(dots$alpha)) alpha <- rep(dots$alpha,length.out=length(n)) else alpha <- rep(1,length(n))
				if (length(dots$beta)) beta <- rep(dots$beta,length.out=length(n)) else beta <- rep(1,length(n))
				
				for (i in 1:length(n)) {
					tmp <- gauss.quad.prob(n[i], dist, l[i], u[i], mu[i], sigma[i], alpha[i], beta[i])
					p[[i]] <- tmp$nodes
					w[[i]] <- tmp$weights
				}
			} else {
				if (length(dots$mean)) mean <- rep(dots$mean,length.out=length(n)) else mean <- rep(0,length(n))
				if (length(dots$sd)) sd <- rep(dots$sd,length.out=length(n)) else sd <- rep(1,length(n))
				for (i in 1:length(n)) {
					if (length(dots$mean) |length(dots$sd)) {
						p[[i]] <- seq(mean[i]-(3*sd[i]),mean[i]+(3*sd[i]),length.out=n[i])
					} else {
						p[[i]] <- seq(-4,4,length.out=n[i])
					}
					if (normal.wt==TRUE) w[[i]] <- dnorm(p[[i]], mean[i], sd[i]) else w[[i]] <- rep(1,n[i])
				}
			}
			
		} else {
			w <- weight
			if (is.list(weight)) {
				p <- vector("list",length(weight))
				for (i in 1:length(p)) {
					p[[i]] <- seq(-4,4,length.out=length(weight[[i]]))
				}
			} else {
				p <- seq(-4,4,length.out=length(weight))
			}
		}
	} else {
		if (missing(weight)) {
			p <- theta
			if (is.list(theta)) {
				nt <- length(theta)
				n <- unlist(lapply(theta,length))
				if (length(dots$mean)) mean <- rep(dots$mean,length.out=nt) else mean <- rep(0,nt)
				if (length(dots$sd)) sd <- rep(dots$sd,length.out=nt) else sd <- rep(1,nt)
				
				w <- vector("list", length(theta))
				for (i in 1:length(w)) {
					if (normal.wt==TRUE) w[[i]] <- dnorm(theta[[i]], mean[i],sd[i]) else w[[i]] <- rep(1,length(theta[[i]]))
				}
			} else {
				nt <- 1
				n <- length(theta)
				if (length(dots$mean)) mean <- dots$mean[1] else mean <- 0
				if (length(dots$sd)) sd <- dots$sd[1] else sd <- 1
				
				if (normal.wt==TRUE) w <- dnorm(theta,mean,sd) else w <- rep(1,length(theta))
			}
		} else {
			if (is.list(theta)) {
				for (i in 1:length(theta)) {
					if (nrow(as.matrix(theta[[i]]))!=nrow(as.matrix(weight[[i]])) ) {
						stop(paste("The theta values and weights in list element {",i,"} do not have the same length.",sep=""))
					}
				}
			} else {
				if (nrow(as.matrix(theta))!=nrow(as.matrix(weight)) ) {
					stop("theta and weight must have the same length")
				}
			}
			p <- theta
			w <- weight
		}
	} 
	
	if (is.list(p)) {
		if (length(p)==1) {
			p <- p[[1]]
			w <- w[[1]]
		} else {
			p <- expand.grid(p)
			names(p) <- paste("theta",1:ncol(p),sep="")
			tmp <- expand.grid(w)
			w <- tmp[,1]
			for (i in 2:length(n)) {
				w <- w*tmp[,i]
			}
		}
	} else {
		if (dimensions>1) {
			tmp.p <- tmp.w <- vector("list",dimensions)
			for (i in 1:dimensions) {
				tmp.p[[i]] <- p
				tmp.w[[i]] <- w
			}
			p <- expand.grid(tmp.p)
			names(p) <- paste("theta",1:ncol(p),sep="")
			tmp <- expand.grid(tmp.w)
			w <- tmp[,1]
			for (i in 2:dimensions) {
				w <- w*tmp[,i]
			}
		}
	}
	
	return(list(points=p,weights=w))
}