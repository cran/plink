plot.irt.prob <- function(x, y, ..., type, separate, combine, items, item.names, item.nums, panels) {
	
	graphics.off()
	if (exists(".SavedPlots",where=1)) rm(".SavedPlots",pos=1)
	
	dimensions <- x@dimensions
	if (missing(type)) type <- "wireframe"
	if (missing(item.nums)) item.nums <- TRUE
	if(missing(items)) {
		items <- 1:length(x@p.cat)
	} else {
		tmp <- rep(1:length(x@p.cat),x@p.cat)
		tmp1 <- 1:length(tmp)
		tmp1 <- tmp1[tmp %in% items]
		x@prob <- x@prob[,c(1:dimensions,tmp1+dimensions)]
		x@p.cat <- x@p.cat[items]
	}
	
	if(!missing(separate)) {
		if (separate==TRUE) cat <- rep(1,ncol(x@prob)-x@dimensions) else cat <- x@p.cat
	} else {
		separate <- FALSE
	}
	
	if (!missing(combine)) {
		cat <- combine 
		if (sum(combine)<sum(x@p.cat)) {
			warning("{combine} did not identify all items. The specified subset will be plotted.")
			x@prob <- x@prob[,1:(sum(cat)+dimensions)]  
		} else if (sum(combine)>sum(x@p.cat)) {
			warning("{combine} identified too many items. The original item categories will be plotted.")
			cat <- x@p.cat
		} 
	} else {
		if(separate==FALSE) cat <- x@p.cat
	}
	
	theta <- as.matrix(x@prob[,1:dimensions])
	nt <- nrow(theta)
	ni <- length(cat)
	if (ncol(x@prob)>(dimensions+1)) {
		sx <- stack(x@prob[,-c(1:dimensions)]) 
	} else {
		sx <- data.frame(values=x@prob[,-c(1:dimensions)],ind=factor(rep(colnames(x@prob)[dimensions+1],nt)))
	}
	id <- NULL
	cid <- NULL
	
	for (i in 1:ni) {
		id <- c(id, rep(i,nt*cat[i]))
		for (j in 1:cat[i]) {
			cid <- c(cid, rep(j,nt))
		}
	}
	
	if (missing(item.names)) {
		if (!missing(combine)) {
		
		} else {
			if (separate==TRUE) {
				nms <- NULL
				for (i in 1:length(items)) {
					if (x@p.cat[i]==1) {
						nms <- c(nms,paste("Item",items[i]))
					} else {
						for (j in 1:x@p.cat[i]) {
							nms <- c(nms, paste("Item ",items[i],".",j,sep=""))
						}
					}
				}
				id <- factor(id,seq(1:ni),nms)
			} else {
				id <- factor(id,seq(1:ni),paste("Item",items))
			}
		}
	} else { 
		id <- factor(id,seq(1:ni),item.names)
	}
	if (sum(cat)>1) {
		tmp <- theta
		for (i in 2:sum(cat)) {
			tmp <- rbind(tmp,theta)
		}
		theta <- tmp
	}
	
	out <- cbind(theta,id,cid,sx)
	colnames(out)[1:dimensions] <- paste("theta",1:dimensions,sep="")
	if (dimensions>1) {
		tmp <- paste(names(out)[1:2],collapse="+")
		if (dimensions>2) {
			for (i in 1:(dimensions-2)) {
				out[,i+2] <- as.factor(out[,i+2])
				
			}
			tmp2 <- paste(c(names(out)[3:dimensions],"id"),collapse="+") 
		} else {
			tmp2 <- "id"
		}
		form<- as.formula(paste("values~",tmp,"|",tmp2,sep="")) 
	}
	
	# Custom strip
	if (dimensions>2) {
		vn <- c(paste("theta[",3:dimensions,"]",sep=""),"")
		strip.irt.prob <- function(which.given, which.panel, var.name, factor.levels, ...) {
			vn <- vn[which.given]
			fl <- factor.levels
			if (which.given<=dimensions-2) {
				expr <- paste(vn, "==", fl, collapse = ",")
				expr <- paste("expression(", expr, ")", sep = "")
				fl <- eval(parse(text = expr))
			}
			strip.default(which.given, which.panel, vn, fl, ...)
		}
	} else {
		strip.irt.prob <- TRUE
	}
	
	# Determine the number of panels to print
	if (dimensions<3) np <- ni else np <- ni*length(unique(out[,1]))^(dimensions-2)
	if (missing(panels)) {
		if (np<20) panels <- np else panels <- 20
	}
	
	cols <- c("lightpink1", "darkseagreen1", "burlywood1", "cadetblue3", "yellow", "darkorchid2", "coral", "seagreen4")
	if (dimensions<3) cols <- "lightblue" else cols <- c(cols[(dimensions-2):1],"lightblue")
	
	if (dimensions==1) {
		if (np>panels) {
			if (Sys.info()["sysname"] == "Windows") {
				cat("Use PgUp and PgDn to view different plot pages\n")
				windows(record=TRUE)
			}
			xyplot(values~theta1|id,out,type="l",
			as.table=TRUE,
			ylab="Probability",
			xlab=expression(paste(theta)),
			groups=cid,
			par.strip.text=list(cex=0.7),
			par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1) ) ,
			strip=strip.irt.prob,
			layout=c(0,panels),
			ylim=c(-0.05,1.05),...)
		} else {
			xyplot(values~theta1|id,out,type="l",
			as.table=TRUE,
			ylab="Probability",
			xlab=expression(paste(theta)),
			groups=cid,
			par.strip.text=list(cex=0.7),
			par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1) ) ,
			strip=strip.irt.prob,
			ylim=c(-0.05,1.05),...)
		}
	} else if (dimensions>1) {
		if (type=="vectorplot1"|type=="vectorplot2"|type=="vectorplot3") {
			if (dimensions==2) {
				if (!missing(items)) {
					a <- as.matrix(x@pars$a)[items,]
					b <- as.matrix(as.matrix(x@pars$b)[items,])
				} else {
					a <- as.matrix(x@pars$a)
					b <- as.matrix(x@pars$b)
				}
				num <- apply(b,1,sum,na.rm=T)
				mdisc <- sqrt(apply(a^2,1,sum,na.rm=T))
				den.con <- apply(!is.na(b),1,sum)
				dcos <- a/matrix(mdisc,length(mdisc),ncol(a))
				mdiff <- -apply(b,1,sum,na.rm=T)/(den.con*mdisc)
				
				# X and Y coordinates for each vector
				x1 <- dcos[,1]*mdiff
				y1 <- dcos[,2]*mdiff
				x2 <- mdisc*dcos[,1]+x1
				y2 <- mdisc*dcos[,2]+y1
				
				# X and Y coordinates for item numbers
				md <- mdiff+0.05
				md[mdiff<0] <- mdiff[mdiff<0]-0.05
				if (type!="vectorplot3")  md <- mdiff-0.05
				x3 <- dcos[,1]*md
				y3 <- dcos[,2]*md
				
				# Compute the reference composite
				ref <- abs(eigen(t(a)%*%a)$vectors[,1])
				
				par(mai=c(0.25,0.25,0.25,0.25))
				if (type=="vectorplot3") {
					xmin <- floor(min(x1))
					xmax <- ceiling(max(x1))
					ymin <- floor(min(y1))
					ymax <- ceiling(max(y1))
					if (xmin<0) {
						if (abs(min(x1)) < abs(xmin+.5)) xmin <- xmin+.5
					}
					if (xmax<0) {
						if (abs(max(x1)) < abs(xmax-.5)) xmax<- xmax-.5
					}
					if (ymin<0) {
						if (abs(min(y1)) < abs(ymin+.5)) ymin <- ymin+.5
					}
					if (ymax<0) {
						if (abs(max(y1)) < abs(ymax-.5)) ymax<- ymax-.5
					}
				} else {
					xmin <- floor(min(x1))
					xmax <- ceiling(max(x2))
					ymin <- floor(min(y1))
					ymax <- ceiling(max(y2))
					if (xmin<0) {
						if (abs(min(x1)) < abs(xmin+.5)) xmin <- xmin+.5
					}
					if (xmax<0) {
						if (abs(max(x2)) < abs(xmax-.5)) xmax<- xmax-.5
					}
					if (ymin<0) {
						if (abs(min(y1)) < abs(ymin+.5)) ymin <- ymin+.5
					}
					if (ymax<0) {
						if (abs(max(y2)) < abs(ymax-.5)) ymax<- ymax-.5
					}
				}
				xmm <- c(xmin,xmax)
				ymm <- c(ymin,ymax)
				
				dots <- list(...)
				if (length(dots$xlim)) xmm <- dots$xlim
				if (length(dots$ylim)) ymm <- dots$ylim
				
				plot(50,50, axes=FALSE,xlim=xmm,ylim=ymm,
				mar=rep(0.1,4),xlab="",ylab="")
				axis(1,pos=0,at=seq(xmm[1],xmm[2],0.5))
				axis(2,pos=0,at=seq(ymm[1],ymm[2],0.5))
				if (type=="vectorplot1") {
					arrows(x1,y1,x2,y2,length=.1)
				} else if (type=="vectorplot2") {
					arrows(x1,y1,x2,y2,length=.1)
					arrows(max(xmm[1],ymm[1])*ref[1],max(xmm[1],ymm[1])*ref[2],min(xmm[2],ymm[2])*ref[1],
						min(xmm[2],ymm[2])*ref[2],lwd=3,length=.2)
				} else if (type=="vectorplot3") {
					tmp <- rep(0,length(mdisc))
					arrows(tmp,tmp,x1,y1,length=.1)
				}
				if (item.nums==TRUE) text(x3,y3,1:length(mdisc), col=2)
				text(xmm[2],0.06, expression(paste(theta[1])))
				text(0.06,ymm[2], expression(paste(theta[2])))
			} else {
				stop("Vector plots are not supported for more than two dimensions")
			}
		} else {
			if (np>panels) {
				if (Sys.info()["sysname"] == "Windows") {
					cat("Use PgUp and PgDn to view different plot pages\n")
					windows(record=TRUE)
				}
				mlab <- "Probability"
				eval(parse(text=paste(type,"(form,out, as.table=TRUE, 
				zlab=list(label=mlab, rot=90),
				xlab=expression(paste(theta[1])),
				ylab=expression(paste(theta[2])),
				groups=cid,
				par.strip.text=list(cex=0.7),
				par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1) ) ,
				strip=strip.irt.prob,
				layout=c(0,panels),
				zlim=c(0,1),...)",sep="")))
			} else {
				mlab <- "Probability"
				eval(parse(text=paste(type,"(form,out,
				as.table=TRUE, 
				zlab=list(label=mlab, rot=90),
				xlab=expression(paste(theta[1])),
				ylab=expression(paste(theta[2])),
				groups=cid,
				par.strip.text=list(cex=0.7),
				par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1) ) ,
				strip=strip.irt.prob,
				layout=c(0,panels),
				zlim=c(0,1),...)",sep="")))
			}
		}
	}
}

plot.sep.pars <- function(x, y, ...) {
	x <- mixed(x, ...)
	plot(x, ...)
}

plot.irt.pars <- function(x, y, ...) {
	if (x@groups>1) {
		warning("There is more than one group in {x}. No plots were produced.")
	} else {
		dots <- list(...)
		if (length(dots$theta)) {
			theta <- dots$theta 
		} else {
			if (x@dimensions==1) {
				theta <- seq(-4,4,.05) 
			} else if (x@dimensions %in% 2:3) {
				theta <- seq(-4,4,.5)
			} else {
				theta <- -4:4
			}
		}
		if (length(dots$catprob)) catprob <- dots$catprob else catprob <- FALSE
		if (length(dots$incorrect)) incorrect <- dots$incorrect else incorrect <- FALSE
		if (length(dots$D)) D <- dots$D else D <- 1.7
		tmp <- mixed(x, theta=theta, catprob=catprob, location=x@location, incorrect=incorrect, D=D)
		plot(tmp, ...)
	}
}