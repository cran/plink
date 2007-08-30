as.weight <- function(theta1=NULL, weight1=NULL, theta2=NULL, weight2=NULL, normal.wt=FALSE) {
	if (is.null(theta1)) {
		if (is.null(weight1)) {
			n <- seq(-4,4,.05)
			if (normal.wt==TRUE) w <- dnorm(n) else w <- rep(1,length(n))
			n <- cbind(n,n)
			w <- cbind(w,w)
		} else {
			if (is.null(weight2)) {
				n <- seq(-4,4,length.out=length(weight1))
				n <- cbind(n,n)
				w <- cbind(weight1,weight1)
			} else {
				n1 <- seq(-4,4,length.out=length(weight1))
				n2 <- seq(-4,4,length.out=length(weight2))
				n <- cbind(n1,n2)
				w <- cbind(weight1,weight2)
			}
		}
	} else {
		if (is.null(theta2)) {
			if (is.null(weight1)) {
				n <- cbind(theta1,theta1)
				if (normal.wt==TRUE) w <- dnorm(theta1) else w <- rep(1,length(theta1))
				w <- cbind(w,w)
			} else {
				n <- cbind(theta1,theta1)
				w <- cbind(weight1,weight1)
			}
		} else {
			if (is.null(weight1)) {
				n <- cbind(theta1,theta2)
				if (normal.wt==TRUE) w1 <- dnorm(theta1) else w1 <- rep(1,length(theta1))
				if (normal.wt==TRUE) w2 <- dnorm(theta2) else w2 <- rep(1,length(theta2))
				w <- cbind(w1,w2)
			} else {
				if (is.null(weight2)) {
					n <- cbind(theta1,theta2)
					w1 <- weight1
					if (normal.wt==TRUE) w2 <- dnorm(theta2) else w2 <- rep(1,length(theta2))
					w <- cbind(w1,w2)
				} else {
					n <- cbind(theta1,theta2)
					w <- cbind(weight1,weight2)
				}
			}
		}
	}
	return(list(points=n,weights=w))
}
