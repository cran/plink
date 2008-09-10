as.weight <- function(n=40, theta1=NULL, weight1=NULL, theta2=NULL, weight2=NULL, normal.wt=FALSE, ...) {
	if (is.null(theta1)) {
		if (is.null(weight1)) {
			dots <- list(...)
			if (length(dots$dist)) {
				tmp <- gauss.quad.prob(n, ...)
			} else {
				tmp <- gauss.quad.prob(n, dist="normal", ...)
			}
			p <- cbind(tmp$nodes,tmp$nodes)
			w <- cbind(tmp$weights,tmp$weights)
		} else {
			if (is.null(weight2)) {
				p <- seq(-4,4,length.out=length(weight1))
				p <- cbind(n,n)
				w <- cbind(weight1,weight1)
			} else {
				if (length(weight1)!=length(weight2)) {
					stop("weight1 and weight2 must have the same length")
				}
				n1 <- seq(-4,4,length.out=length(weight1))
				n2 <- seq(-4,4,length.out=length(weight2))
				p <- cbind(n1,n2)
				w <- cbind(weight1,weight2)
			}
		}
	} else {
		if (is.null(theta2)) {
			if (is.null(weight1)) {
				p <- cbind(theta1,theta1)
				if (normal.wt==TRUE) w <- dnorm(theta1) else w <- rep(1,length(theta1)/length(theta1))
				w <- cbind(w,w)
			} else {
				if (length(theta1)!=length(weight1)) {
					stop("theta1 and weight1 must have the same length")
				}
				p <- cbind(theta1,theta1)
				w <- cbind(weight1,weight1)
			}
		} else {
			if (length(theta1)!=length(theta2)) {
				stop("theta1 and theta2 must have the same length")
			}
			if (is.null(weight1)) {
				p <- cbind(theta1,theta2)
				if (normal.wt==TRUE) {
					w1 <- dnorm(theta1) 
					w2 <- dnorm(theta2)
				} else {
					w1 <- rep(1,length(theta1))
					w2 <- rep(1,length(theta2))
				}
				w <- cbind(w1,w2)
			} else {
				if (length(theta1)!=length(weight1)) {
					stop("theta1 and weight1 must have the same length")
				}
				if (is.null(weight2)) {
					p <- cbind(theta1,theta2)
					w1 <- weight1
					if (normal.wt==TRUE) w2 <- dnorm(theta2) else w2 <- rep(1,length(theta2))
					w <- cbind(w1,w2)
				} else {
					tmp <- (length(weight1) + length(weight2))/2
					if (length(theta1)!=tmp) {
						stop("theta1, weight1, and weight2 must have the same length")
					}
					p <- cbind(theta1,theta2)
					w <- cbind(weight1,weight2)
				}
			}
		}
	}
	return(list(points=p,weights=w))
}