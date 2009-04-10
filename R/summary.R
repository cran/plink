summary.irt.pars <- function(object, ..., descrip=FALSE) {
	object <- sep.pars(object)
	summary(object, descrip=descrip)
}

summary.list <- function(object, ..., descrip=FALSE) {
	if (is.list(object[[1]])) {
		if (is.link(object[[1]][[1]])) {
			nms <- matrix(c("MM","MS","HB","SL","RM","Mean/Mean","Mean/Sigma","Haebara","Stocking-Lord","Reckase-Martineau"),5,2)
			for (i in 1:length(object[[1]])) {
				summary(object[[1]][[i]],descrip=descrip)
			}
		} else {
			stop("The objects in the list are not of class {link}")
		}
	}  else {
		if (is.link(object[[1]]) & !is.link(object[[2]])) {
			summary(object[[1]],descrip=descrip)
		} else if (is.link(object[[1]]) & is.link(object[[2]])) {
			for (i in 1:length(object)) {
				summary(object[[i]],descrip=descrip)
			}
		} else {
			for (i in 1:length(object)) {
				if (is.sep.pars(object[[i]])) {
					cat("--------",names(object[i]),"--------\n")
					obj <- object[[i]]
				} else if (is.irt.pars(object[[i]])) {
					cat("--------",names(object[i]),"--------\n")
					obj <- sep.pars(object[[i]])
				} else {
					stop("The objects in the list are not of class{sep.pars}, {irt.pars}, or {link}")
				}
				summary(obj, descrip=descrip)
				cat("\n")
			}
		}
	}
}

summary.sep.pars <- function(object, ..., descrip=FALSE){
	if (object@loc.out==TRUE) {
		pm <- as.poly.mod(object@n[1],object@model,object@items)
		object <- sep.pars(list(object@a,object@b,object@c),object@cat,pm,object@dimensions,TRUE,FALSE)
	}
	dimensions <- object@dimensions
	
	.sd <- function(x) {
		z <- x[!is.na(x)]
		out <- sqrt(sum((z-mean(z))^2)/length(z))
		return(out)
	}
	
	cat(paste("Total Number of Items:",object@n[1],"\n\n"))
	if (descrip==TRUE) {
		if (length(object@model)>1) {
			 a1 <- vector("list",dimensions)
			for (i in 1:length(object@mod.lab)) {
				items <- object@items[[i]]
				catk <- object@cat[items]
				
				if (object@model[i]=="nrm"|object@model[i]=="mcm") {
					mcatk <- max(catk)
					for (m in 1:dimensions) {
						a1[[m]] <- c(a1[[m]],as.vector(object@a[items,(((m-1)*mcatk)+1):(m*mcatk)]))
					}
				} else {
					for (m in 1:object@dimensions) {
						a1[[m]] <- c(a1[[m]],object@a[items,m])
					}
				}
			}
			a <- NULL
			for (i in 1:dimensions) {
				tmp <- a1[[i]]
				a <- cbind(a,c(length(tmp[!is.na(tmp)]),mean(tmp,na.rm=TRUE),.sd(tmp),min(tmp,na.rm=TRUE),max(tmp,na.rm=TRUE)))
			}
			b1 <- object@b
			c1 <- object@c
			b <- c(length(b1[!is.na (b1)]),mean(b1,na.rm=TRUE),.sd(b1),min(b1,na.rm=TRUE),max(b1,na.rm=TRUE))
			c <- c(length(c1[!is.na (c1)]),mean(c1,na.rm=TRUE),.sd(c1),min(c1,na.rm=TRUE),max(c1,na.rm=TRUE))
			mdisc <- sqrt(apply(as.matrix(object@a)^2,1,sum,na.rm=TRUE))
			mdif <- -b1/mdisc
			mdc <- c(length(mdisc[!is.na(mdisc)]),mean(mdisc,na.rm=TRUE),.sd(mdisc),min(mdisc,na.rm=TRUE),max(mdisc,na.rm=TRUE))
			mdf <- c(length(mdif[!is.na(mdif)]),mean(mdif,na.rm=TRUE),.sd(mdif),min(mdif,na.rm=TRUE),max(mdif,na.rm=TRUE))
			
			des <- cbind(a,b,c,mdc,mdf)
			rownames(des) <- c("N","Mean","SD","Min","Max")
			if (object@dimensions==1) {
				des <- des[,1:3]
				colnames(des) <- c("a","b","c")
			} else {
				colnames(des) <- c(paste("a",1:object@dimensions,sep=""),"d","c","MDISC","MDIF")
			}
			des <- round(des,4)
			cat("All Item Parameter Descriptives:\n")
			print(format(des,justify="right"),quote=FALSE)
			cat("\n")
		}
	}
	
	cat(paste("Number of Dichotomous Items:",object@n[2],"\n"))
	if (object@n[2]>=1) cat(paste("Dichotomous Model:",object@mod.lab[object@model=="drm"],"\n\n")) else cat("\n")
	if (descrip==TRUE) {
		if (object@n[2]>1) {
			a <- NULL
			for (j in 1:object@dimensions) {
				a1 <- object@a[object@items$drm,j]
				a <- cbind(a,c(length(a1[!is.na (a1)]), mean(a1,na.rm=TRUE),.sd(a1),min(a1,na.rm=TRUE),max(a1,na.rm=TRUE)))
			}
			b1 <- object@b[object@items$drm,]
			c1 <- object@c[object@items$drm,]
			b <- c(length(b1[!is.na (b1)]),mean(b1,na.rm=TRUE),.sd(b1),min(b1,na.rm=TRUE),max(b1,na.rm=TRUE))
			c <- c(length(c1[!is.na (c1)]),mean(c1,na.rm=TRUE),.sd(c1),min(c1,na.rm=TRUE),max(c1,na.rm=TRUE))
			mdisc <- sqrt(apply(as.matrix(object@a[object@items$drm,])^2,1,sum,na.rm=TRUE))
			mdif <- -b1/mdisc
			mdc <- c(length(mdisc[!is.na(mdisc)]),mean(mdisc,na.rm=TRUE),.sd(mdisc),min(mdisc,na.rm=TRUE),max(mdisc,na.rm=TRUE))
			mdf <- c(length(mdif[!is.na(mdif)]),mean(mdif,na.rm=TRUE),.sd(mdif),min(mdif,na.rm=TRUE),max(mdif,na.rm=TRUE))
			
			des <- cbind(a,b,c,mdc,mdf)
			rownames(des) <- c("N","Mean","SD","Min","Max")
			if (object@dimensions==1) {
				des <- des[,1:3]
				colnames(des) <- c("a","b","c")
			} else {
				colnames(des) <- c(paste("a",1:object@dimensions,sep=""),"d","c","MDISC","MDIF")
			}
			des <- round(des,4)
			cat(paste(object@mod.lab[object@model=="drm"],"- Item Parameter Descriptives:\n"))
			print(format(des,justify="right"),quote=FALSE)
			cat("\n")
		} else if (object@n[2]==1) { 
			cat("There is only one item. No descriptives available.\n\n")
		}
	}
	
	cat(paste("Number of Polytomous Items:",object@n[3],"\n"))
	if (object@n[3]>=1) {
		k <- 4
		for (i in 1:length(object@mod.lab)) {
			if (object@model[i]=="drm") next
			items <- object@items[[i]]
			if (length(object@model[object@model!="drm"])>1) {
				cat(paste(object@mod.lab[i],"\n"))
				cat(paste("Number of",object@mod.lab[i],"Items:",object@n[k],"\n\n"))
				j <- k
				k <- k+1
			} else {
				cat(paste("Polytomous Model:",object@mod.lab[i],"\n\n"))
				j <- 3
			}
			if (descrip==TRUE) {
				if (object@n[j]>1) {
					a1 <- as.matrix(object@a)[items,]
					b1 <- as.matrix(object@b)[items,]
					c1 <- as.matrix(object@c)[items,]
					if (length(c1[is.na(c1)])==length(c1)) c1[is.na(c1)] <- 0
					catk <- object@cat[items]
					
					a <- NULL
					if (object@model[i]=="nrm"|object@model[i]=="mcm") {
						mcatk <- max(catk)
						for (m in 1:object@dimensions) {
							tmp <- a1[,(((m-1)*mcatk)+1):(m*mcatk)]
							a <- cbind(a,c(length(tmp[!is.na(tmp)]),mean(tmp,na.rm=TRUE),.sd(tmp),min(tmp,na.rm=TRUE),max(tmp,na.rm=TRUE)))
						}
					} else {
						for (m in 1:object@dimensions) {
							tmp <- a1[,m]
							a <- cbind(a,c(length(tmp[!is.na(tmp)]),mean(tmp,na.rm=TRUE),.sd(tmp),min(tmp,na.rm=TRUE),max(tmp,na.rm=TRUE)))
						}
					}
					b <- c(length(b1[!is.na (b1)]),mean(b1,na.rm=TRUE),.sd(b1),min(b1,na.rm=TRUE),max(b1,na.rm=TRUE))
					c <- c(length(c1[!is.na (c1)]),mean(c1,na.rm=TRUE),.sd(c1),min(c1,na.rm=TRUE),max(c1,na.rm=TRUE))
					mdisc <- sqrt(apply(as.matrix(a1)^2,1,sum,na.rm=TRUE))
					mdif <- -b1/mdisc
					mdc <- c(length(mdisc[!is.na(mdisc)]),mean(mdisc,na.rm=TRUE),.sd(mdisc),min(mdisc,na.rm=TRUE),max(mdisc,na.rm=TRUE))
					mdf <- c(length(mdif[!is.na(mdif)]),mean(mdif,na.rm=TRUE),.sd(mdif),min(mdif,na.rm=TRUE),max(mdif,na.rm=TRUE))
					
					if (object@model[i]=="mcm") {
						if (object@dimensions==1) {
							des <- cbind(a,b,c)
							colnames(des) <- c("a","c","d")
						} else {
							des <- cbind(a,b,c,mdc,mdf)
							colnames(des) <- c(paste("a",1:object@dimensions,sep=""),"d","c","MDISC","MDIF")
						}
					} else {
						if (object@dimensions==1) {
							des <- cbind(a,b)
							colnames(des) <- c("a","d")
						} else {
							des <- cbind(a,b,mdc,mdf)
							colnames(des) <- c(paste("a",1:object@dimensions,sep=""),"d","MDISC","MDIF")
						}
					}
					rownames(des) <- c("N","Mean","SD","Min","Max")
					des <- round(des,4)
					if (length(object@mod.lab[object@mod.lab!="drm"])>=1) {
						cat(paste(object@mod.lab[i],"- Item Parameter Descriptives:\n"))
					} 
					print(format(des,justify="right"),quote=FALSE)
					cat("\n")
				} else if (object@n[j]==1) { 
					cat("There is only one item. No descriptives available.\n\n")
				}
			}
		}
	}
}

summary.link <- function(object, ..., descrip=FALSE) {
	nms <- matrix(c("MM","MS","HB","SL","RM","Mean/Mean","Mean/Sigma","Haebara","Stocking-Lord","Reckase-Martineau"),5,2)
	tmp.nms <- NULL
	for (j in 1:length(names(object@constants))) {
		tmp.nms <- c(tmp.nms, nms[nms[,1]==names(object@constants)[j],2])
	}
	
	cat("------- ",object@grp.names," -------\n")
	cat("Linking Constants\n\n")
	if (is.list(object@constants[[1]])) {
		for (j in 1:length(tmp.nms)) {
			cat(tmp.nms[j],"\n")
			if (length(object@constants[[j]] )==2) { #Oblique Procrustes
				print(formatC(object@constants[[j]]$A,digits=6,format="f",drop0trailing=FALSE),quote=FALSE)
				cat("\n")
				print(formatC(object@constants[[j]]$m,digits=6,format="f",drop0trailing=FALSE),quote=FALSE)
				cat("\n")
			} else {
				tmp <- object@constants[[j]]
				tmp1 <- formatC(cbind(tmp$A,NA,tmp$T,NA,tmp$K),digits=6,format="f",drop0trailing=FALSE)
				tmp1[,c(ncol(tmp$A)+1,2*ncol(tmp$A)+2)] <- NA
				print(tmp1,na.print="",quote=F)
				cat("\n")
				print(formatC(tmp$m,digits=6,format="f",drop0trailing=FALSE),quote=FALSE)
				cat("\n")
			}
		}
	} else {
		tmp <- NULL
		for (j in 1:length(tmp.nms)) {
			tmp <- rbind(tmp,object@constants[[j]])
		}
		rownames(tmp) <- tmp.nms
		print(formatC(tmp,digits=6,format="f",drop0trailing=FALSE),quote=FALSE)
	}
	cat("\n")
	
	if (descrip==TRUE) {
		cat("Common Item Descriptive Statistics\n\n")
		if (length(object@n)>3) n <- object@n[-3] else n <- object@n
		if (length(n[n>0])==2) {
			n <- n[-1][n[-1]>0]
		} else {
			n <- n[c(2:length(n),1)]
		}
		ml <- object@mod.lab
		if (length(ml)>1) ml <- c(ml, "All")
		for (j in 1:length(n[n>0])) {
			cat("Model:",ml[j],"\n")
			cat("Number of Items:",n[j],"\n\n")
			print(format(object@descriptives[[j]],justify="right"),quote=FALSE)
			cat("\n")
		}
		cat("\n")
	}
	
}