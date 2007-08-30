summary.irt.pars <- function(object, descrip=FALSE, ...) {
	object <- sep.pars(object, ...)
	summary(object, descrip=descrip, ...)
}

summary.list <- function(object, descrip=FALSE, ...) {
	tmp <- object
	if (is.list(tmp[[1]])) {
		if (is.link(tmp[[1]][[1]])) {
			for (i in 1:length(tmp[[1]])) {
				obj <- tmp[[1]][[i]]
				if (is.null(obj@SL)) {
					tmpa <- rbind(obj@MM,obj@MS,obj@HB)
					rownames(tmpa) <- c("MM","MS","HB")
				} else {
					tmpa <- rbind(obj@MM,obj@MS,obj@HB,obj@SL)
					rownames(tmpa) <- c("MM","MS","HB","SL")
				}
				cat(names(tmp[[1]])[i],"\n")
				print(format(tmpa),quote=FALSE)
				cat("\n")
				if (descrip==TRUE) {
					print(obj@descriptives,quote=FALSE)
					cat("\n")
				}
			}
		} else {
			stop("The objects in the list are not of class {link}")
		}
	}  else {
		if (is.link(tmp[[1]]) & !is.link(tmp[[2]])) {
			obj <- tmp[[1]]
			if (is.null(obj@SL)) {
				tmpa <- rbind(obj@MM,obj@MS,obj@HB)
				rownames(tmpa) <- c("MM","MS","HB")
			} else {
				tmpa <- rbind(obj@MM,obj@MS,obj@HB,obj@SL)
				rownames(tmpa) <- c("MM","MS","HB","SL")
			}
			if (obj@base.grp==1) cat("Group1*/Group2\n") else cat("Group1/Group2*\n")
			print(format(tmpa),quote=FALSE)
			cat("\n")
			if (descrip==TRUE) {
				print(obj@descriptives,quote=FALSE)
				cat("\n")
			}
		} else if (is.link(tmp[[1]]) & is.link(tmp[[2]])) {
			for (i in 1:length(tmp)) {
				obj <- tmp[[i]]
				if (is.null(obj@SL)) {
					tmpa <- rbind(obj@MM,obj@MS,obj@HB)
					rownames(tmpa) <- c("MM","MS","HB")
				} else {
					tmpa <- rbind(obj@MM,obj@MS,obj@HB,obj@SL)
					rownames(tmpa) <- c("MM","MS","HB","SL")
				}
				cat(names(tmp)[i],"\n")
				print(format(tmpa),quote=FALSE)
				cat("\n")
				if (descrip==TRUE) {
					print(obj@descriptives,quote=FALSE)
					cat("\n")
				}
			}
		} else {
			for (i in 1:length(tmp)) {
				if (is.sep.pars(tmp[[i]])) {
					cat("Group: ",names(tmp[i]),"\n")
					object <- tmp[[i]]
				} else if (is.irt.pars(tmp[[i]])) {
					cat("Group: ",names(tmp[i]),"\n")
					object <- sep.pars(tmp[[i]])
				} else {
					stop("The objects in the list are not of class{sep.pars}, {irt.pars}, or {link}")
				}
				summary(object, descrip=descrip, ...)
				cat("\n")
			}
		}
	}
}

summary.sep.pars <- function(object, descrip=FALSE, ...) {
	.sd <- function(x) {
		z <- x[!is.na(x)]
		out <- sqrt(sum((z-mean(z,na.rm=T))^2)/length(z))
		return(out)
	}
	
	cat(paste("Number of Items:",object@n[1],"\n\n"))
	cat(paste("Number of Dichotomous Items:",object@n[2],"\n"))
	if (object@n[2]>=1) cat(paste("Dichotomous Model:",object@mod.lab[1],"\n\n")) else cat("\n")
	if (descrip==TRUE) {
		if (object@n[2]>1) {
			if (is.vector(object@a)) a <-object@a[object@items$drm] else a <- object@a[object@items$drm,1]
			if (is.vector(object@b)) b <-object@b[object@items$drm] else b <- object@b[object@items$drm,1]
			if (is.vector(object@c)) c <-object@c[object@items$drm] else c <- object@c[object@items$drm,1]
			a <- a[!is.na(a)]
			b <- b[!is.na(b)]
			c <- c[!is.na(c)]
			dt <- matrix(c(mean(a),mean(b),mean(c),.sd(a),.sd(b),.sd(c),min(a),min(b),min(c),max(a),max(b),max(c)),4,3,byrow=T)
			rownames(dt) <- c("Mean","SD","Min","Max")
			colnames(dt) <- c("a","b","c")
			dt <- round(dt,4)
			cat(paste(object@mod.lab[1],"Item Parameter Descriptives:\n"))
			print(format(dt),quote=F)
			cat("\n")
		} else if (object@n[2]==1) { 
			cat("There is only one item. No descriptives available.\n\n")
		}
	}
	
	cat(paste("Number of Polytomous Items:",object@n[3],"\n"))
	if (object@n[3]>=1) {
		if (object@model[1]=="drm") start <- 2 else start <- 1
		for (i in start:length(object@mod.lab)) {
			items <- object@items[[i]]
			if (length(object@model[object@model!="drm"])>1) {
				if (i==start) cat("\n")
				if (start==1) ml <- 3 else ml <- 2
				cat(paste(object@mod.lab[i],"\n"))
				cat(paste("Number of",object@mod.lab[i],"Items:",object@n[i+ml],"\n\n"))
				j <- i+ml
			} else {
				cat(paste("Polytomous Model:",object@mod.lab[2],"\n\n"))
				j <- 3
			}
			if (descrip==TRUE) {
				if (object@n[j]>1) {
					if (is.vector(object@a)) a <-object@a[items] else a <- object@a[items,]
					b <- object@b[items,]
					if (is.vector(object@c)) c <-object@c[items] else c <- object@c[items,]
					a <- a[!is.na(a)]
					#if (object@location==TRUE) b <- b[,1]
					b <- b[!is.na(b)]
					c <- c[!is.na(c)]
					if (object@mod.lab[i]=="mcm") {
						dt <- matrix(c(mean(a), mean(b), mean(c), .sd(a), .sd(b), .sd(c), min(a), min(b), min(c), max(a), max(b), max(c)),4,3,byrow=T)
						colnames(dt) <- c("a","c","d")
					} else {
						dt <- matrix(c(mean(a), mean(b), sd(a), sd(b), min(a), min(b), max(a), max(b)),4,2,byrow=T)
						if (object@mod.lab[i]=="nrm") colnames(dt) <- c("a","c") else colnames(dt) <- c("a","b")
					}
					rownames(dt) <- c("Mean","SD","Min","Max")
					dt <- round(dt,4)
					if (length(object@mod.lab[object@mod.lab!="drm"])>1) {
						cat(paste(object@mod.lab[i],"Item Parameter Descriptives:\n"))
					} else {
						cat("Polytomous Item Parameter Descriptives:\n")
					}
					print(format(dt),quote=F)
					cat("\n")
				} else if (object@n[j]==1) { 
					cat("There is only one item. No descriptives available.\n\n")
				}
			}
		}
	}
}

summary.link <- function(object, descrip=FALSE, ...) {
	if (is.null(object@SL)) {
		tmp <- rbind(object@MM,object@MS,object@HB)
		rownames(tmp) <- c("MM","MS","HB")
	} else {
		tmp <- rbind(object@MM,object@MS,object@HB,object@SL)
		rownames(tmp) <- c("MM","MS","HB","SL")
	}
	if (object@base.grp==1) cat("Group1*/Group2\n") else cat("Group1/Group2*\n")
	print(format(tmp),quote=FALSE)
	cat("\n")
	if (descrip==TRUE) print(object@descriptives,quote=FALSE)
}
