##   Summarize the item parameters in an {irt.pars} object
summary.irt.pars <- function(object, ..., descrip=FALSE) {

	object <- sep.pars(object)
	summary(object, descrip=descrip)
	
}


##   Summarize a list that may contain item parameters or
##   the output from running {plink}
summary.list <- function(object, ..., descrip=FALSE) {

	##   Check to see if the first list element is also a list.
	##   The only scenario where summary applies is when
	##   the first object contains the output from {plink}
	##   for multiple groups and rescaled item parameters
	##   and/or ability estimates were returned
	if (is.list(object[[1]])) {
	
		##   Check to see if the first object in this list is of class {link}
		if (is.link(object[[1]][[1]])) {
		
			##   Create a matrix containing the shortened name
			##   and the descrtiptive name for all the linking methods
			##   (both unidimensional and multidimensional)
			nms <- matrix(c("MM","MS","HB","SL","LS","Mean/Mean","Mean/Sigma","Haebara","Stocking-Lord","Least Squares"),5,2)
			
			##   Loop through all of the {link} objects and summarize
			##   the linking constants and (if specified) the descriptive
			##   statistics for the common items
			for (i in 1:length(object[[1]])) {
				summary(object[[1]][[i]],descrip=descrip)
			}
			
			##   Check to see if ability estimates are included in the output
			##   If yes, compute descriptive statistics for each group
			if (length(object$ability)) {
			
				##   Initialize object to store the descriptives
				tmp.ability <- NULL
				
				for (i in 1:length(object$ability)) {
					tmp <- object$ability[[i]]
					
					if (is.matrix(tmp)|is.data.frame(tmp)) {
						##   Compute the descriptives
						means <- apply(tmp,2,mean,na.rm=T)
						sds <- apply(tmp,2,sd,na.rm=T)
						tmp.min <- apply(tmp,2,min,na.rm=T)
						tmp.max <- apply(tmp,2,max,na.rm=T)
						 
						##   Compile the descriptives and create column names
						tmp1 <- rbind(means,sds,tmp.min,tmp.max)
						th <- paste("theta",1:ncol(tmp),sep="")
						colnames(tmp1) <- paste(names(object$ability)[i],":",th,sep="")
						
					} else if (is.vector(tmp)) {
						##   Compute the descriptives
						means <- mean(tmp,na.rm=T)
						sds <- sd(tmp,na.rm=T)
						tmp.min <- min(tmp,na.rm=T)
						tmp.max <- max(tmp,na.rm=T)
						
						##   Compile the descriptives and create column names
						tmp1 <- as.matrix(c(means,sds,tmp.min,tmp.max))
						colnames(tmp1) <- names(object$ability)[i]
					}
					
					tmp.ability <- cbind(tmp.ability, tmp1)
				}
				
				##   Print ability descriptives
				rownames(tmp.ability) <- c("Mean","SD","Min","Max")
				tmp.ability <- round(tmp.ability,4)
				cat("Ability Descriptive Statistics\n\n")
				print(format(as.data.frame(tmp.ability),nsmall=4),quote=FALSE)
				cat("\n")
			}
			
		} else {
			stop("The objects in the list are not of class {link}")
		}
		
	##   If the first element is not a list
	}  else {
	
		##   If the first object is of class {link} and the second
		##   element is not. This is the case where there were only
		##   two groups run through {plink} and rescaled item
		##   parameters and/or ability estimates were returned
		if (is.link(object[[1]]) & !is.link(object[[2]])) {
			summary(object[[1]],descrip=descrip)
			
			##   Check to see if ability estimates are included in the output
			##   If yes, compute descriptive statistics for each group
			if (length(object$ability)) {
			
				##   Initialize object to store the descriptives
				tmp.ability <- NULL
				
				for (i in 1:length(object$ability)) {
					tmp <- object$ability[[i]]
					
					if (is.matrix(tmp)|is.data.frame(tmp)) {
						##   Compute the descriptives
						means <- apply(tmp,2,mean,na.rm=T)
						sds <- apply(tmp,2,sd,na.rm=T)
						tmp.min <- apply(tmp,2,min,na.rm=T)
						tmp.max <- apply(tmp,2,max,na.rm=T)
						 
						##   Compile the descriptives and create column names
						tmp1 <- rbind(means,sds,tmp.min,tmp.max)
						th <- paste("theta",1:ncol(tmp),sep="")
						colnames(tmp1) <- paste(names(object$ability)[i],":",th,sep="")
						
					} else if (is.vector(tmp)) {
						##   Compute the descriptives
						means <- mean(tmp,na.rm=T)
						sds <- sd(tmp,na.rm=T)
						tmp.min <- min(tmp,na.rm=T)
						tmp.max <- max(tmp,na.rm=T)
						
						##   Compile the descriptives and create column names
						tmp1 <- as.matrix(c(means,sds,tmp.min,tmp.max))
						colnames(tmp1) <- names(object$ability)[i]
					}
					
					tmp.ability <- cbind(tmp.ability, tmp1)
				}
				
				##   Print ability descriptives
				rownames(tmp.ability) <- c("Mean","SD","Min","Max")
				tmp.ability <- round(tmp.ability,4)
				cat("Ability Descriptive Statistics\n\n")
				print(format(as.data.frame(tmp.ability),nsmall=4),quote=FALSE)
				cat("\n")
			}
			
		##   This is the case where the output from {plink}
		##   is for three or more groups, but no rescaled
		##   parameters were returned
		} else if (is.link(object[[1]]) & is.link(object[[2]])) {
			for (i in 1:length(object)) {
				summary(object[[i]],descrip=descrip)
			}
			
		##   This is the case where the list includes {sep.pars}
		##   or {irt.pars} objects. In these instances, create a
		##   dividing line for each group
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



##   Summarize the item parameters in a {sep.pars} object
summary.sep.pars <- function(object, ..., descrip=FALSE){
	
	if (object@location==TRUE) {
		pm <- as.poly.mod(object@n[1],object@model,object@items)
		object <- sep.pars(list(object@a,object@b,object@c),object@cat,pm,object@dimensions,TRUE,FALSE)
	}
	
	dimensions <- object@dimensions
	
	##   The default SD function in R uses n-1 in the denominator
	##   This function uses n in the denominator
	.sd <- function(x) {
		z <- x[!is.na(x)]
		out <- sqrt(sum((z-mean(z))^2)/length(z))
		return(out)
	}
	
	##   Print the total number of items
	cat(paste("Total Number of Items:",object@n[1],"\n\n"))
	
	##   Print the descriptive statistics for all of the item parameters
	if (descrip==TRUE) {
		if (length(object@model)>1) {
		
			##   Initialize an object to store the slope parameters
			##   A special object is created for these parameters
			##   because there may be multiple dimensions
			 a1 <- vector("list",dimensions)
			
			##   Loop through all of the response models
			for (i in 1:length(object@mod.lab)) {
			
				##   Identify the items associated with the given model
				items <- object@items[[i]]
				
				##   Identify the number of response categories for these items
				catk <- object@cat[items]
				
				if (object@model[i]=="nrm"|object@model[i]=="mcm") {
					mcatk <- max(catk)
					
					##   Extract the slopes for each dimension
					for (m in 1:dimensions) {
						a1[[m]] <- c(a1[[m]],as.vector(object@a[items,(((m-1)*mcatk)+1):(m*mcatk)]))
					}
				} else {
					##   Extract the slopes for each dimension
					for (m in 1:object@dimensions) {
						a1[[m]] <- c(a1[[m]],object@a[items,m])
					}
				}
			}
			
			##   Loop through all of the dimensions and identify the number
			##   of slope parameters, then compute the mean, SD, min, and max
			##   of each of these parameters
			a <- NULL
			for (i in 1:dimensions) {
				tmp <- a1[[i]]
				a <- cbind(a,c(length(tmp[!is.na(tmp)]),mean(tmp,na.rm=TRUE),.sd(tmp),min(tmp,na.rm=TRUE),max(tmp,na.rm=TRUE)))
			}
			
			##   Extract the difficulty/threshold/step/category parameters
			b1 <- object@b
			
			##   Extract the lower asymptote parameters
			c1 <- object@c
			
			##   Identify the number of parameters, and compute the 
			##   mean, SD, min, and max for each parameter type
			b <- c(length(b1[!is.na (b1)]),mean(b1,na.rm=TRUE),.sd(b1),min(b1,na.rm=TRUE),max(b1,na.rm=TRUE))
			c <- c(length(c1[!is.na (c1)]),mean(c1,na.rm=TRUE),.sd(c1),min(c1,na.rm=TRUE),max(c1,na.rm=TRUE))
			
			##   Compute the MDISC and MDIF (these will only
			##   be printed if there are two or more dimensions)
			mdisc <- sqrt(apply(as.matrix(object@a)^2,1,sum,na.rm=TRUE))
			mdif <- -b1/mdisc
			mdc <- c(length(mdisc[!is.na(mdisc)]),mean(mdisc,na.rm=TRUE),.sd(mdisc),min(mdisc,na.rm=TRUE),max(mdisc,na.rm=TRUE))
			mdf <- c(length(mdif[!is.na(mdif)]),mean(mdif,na.rm=TRUE),.sd(mdif),min(mdif,na.rm=TRUE),max(mdif,na.rm=TRUE))
			
			##   Combine the descriptives for each model
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
			print(format(as.data.frame(des),nsmall=4),quote=FALSE)
			cat("\n")
		}
	}
	
	##   Print the number of dichotomous items
	cat(paste("Number of Dichotomous Items:",object@n[2],"\n"))
	
	if (object@n[2]>=1) cat(paste("Dichotomous Model:",object@mod.lab[object@model=="drm"],"\n\n")) else cat("\n")
	
	##   Print the descriptives for the dichotomous items if there is at least one item
	if (descrip==TRUE) {
		if (object@n[2]>1) {
		
			a <- NULL
			for (j in 1:object@dimensions) {
				##   Extract the slopes for the dichotomous items
				a1 <- object@a[object@items$drm,j]
				
				##   Identify the number of slope parameters for each dimension
				##   and compute the mean, SD, min, and max
				a <- cbind(a,c(length(a1[!is.na (a1)]), mean(a1,na.rm=TRUE),.sd(a1),min(a1,na.rm=TRUE),max(a1,na.rm=TRUE)))
			}
			
			##   Extract the difficulty parameters
			b1 <- object@b[object@items$drm,]
			
			##   Extract the lower asymptote parameters
			c1 <- object@c[object@items$drm,]
			
			##   Identify the number of parameters, and compute the 
			##   mean, SD, min, and max for each parameter type
			b <- c(length(b1[!is.na (b1)]),mean(b1,na.rm=TRUE),.sd(b1),min(b1,na.rm=TRUE),max(b1,na.rm=TRUE))
			c <- c(length(c1[!is.na (c1)]),mean(c1,na.rm=TRUE),.sd(c1),min(c1,na.rm=TRUE),max(c1,na.rm=TRUE))
			
			##   Compute the MDISC and MDIF (these will only
			##   be printed if there are two or more dimensions)
			mdisc <- sqrt(apply(as.matrix(object@a[object@items$drm,])^2,1,sum,na.rm=TRUE))
			mdif <- -b1/mdisc
			mdc <- c(length(mdisc[!is.na(mdisc)]),mean(mdisc,na.rm=TRUE),.sd(mdisc),min(mdisc,na.rm=TRUE),max(mdisc,na.rm=TRUE))
			mdf <- c(length(mdif[!is.na(mdif)]),mean(mdif,na.rm=TRUE),.sd(mdif),min(mdif,na.rm=TRUE),max(mdif,na.rm=TRUE))
			
			##   Combine these descriptives then print them
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
			print(format(as.data.frame(des),nsmall=4),quote=FALSE)
			cat("\n")
		} else if (object@n[2]==1) { 
			cat("There is only one item. No descriptives available.\n\n")
		}
	}
	
	##   Print the total number of polytomous items
	cat(paste("Number of Polytomous Items:",object@n[3],"\n"))
	
	##   If there is at least one polytomous item
	if (object@n[3]>=1) {
		
		##   Index the values in the object {n} in the sep.pars object
		k <- 4
		
		##   Loop through all of the polytomous models
		for (i in 1:length(object@mod.lab)) {
			if (object@model[i]=="drm") next
			
			##   Identify the number of items associated with the given model
			items <- object@items[[i]]
			
			##   Print the name of the given model and the associated number of items
			if (length(object@model[object@model!="drm"])>1) {
				cat(paste(object@mod.lab[i],"\n"))
				cat(paste("Number of",object@mod.lab[i],"Items:",object@n[k],"\n\n"))
				j <- k
				k <- k+1
			} else {
				cat(paste("Polytomous Model:",object@mod.lab[i],"\n\n"))
				j <- 3
			}
			
			##   Compute and print the descriptive statistics
			if (descrip==TRUE) {
				if (object@n[j]>1) {
					
					##   Extract the item parameters for the given model
					a1 <- as.matrix(as.matrix(object@a)[items,])
					b1 <- as.matrix(object@b)[items,]
					c1 <- as.matrix(object@c)[items,]
					if (length(c1[is.na(c1)])==length(c1)) c1[is.na(c1)] <- 0
					catk <- object@cat[items]
					
					
					##   Identify the number of slope parameters for each dimension
					##   and compute the mean, SD, min, and max
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
					
					##   Identify the number of threshold/step/category and lower asymptote parameters
					##   and compute the mean, SD, min, and max for each parameter type
					
					b <- c(length(b1[!is.na (b1)]),mean(b1,na.rm=TRUE),.sd(b1),min(b1,na.rm=TRUE),max(b1,na.rm=TRUE))
					c <- c(length(c1[!is.na (c1)]),mean(c1,na.rm=TRUE),.sd(c1),min(c1,na.rm=TRUE),max(c1,na.rm=TRUE))
					
					##   Compute the MDISC and MDIF (these will only
					##   be printed if there are two or more dimensions)
					mdisc <- sqrt(apply(as.matrix(a1)^2,1,sum,na.rm=TRUE))
					mdif <- -b1/mdisc
					mdc <- c(length(mdisc[!is.na(mdisc)]),mean(mdisc,na.rm=TRUE),.sd(mdisc),min(mdisc,na.rm=TRUE),max(mdisc,na.rm=TRUE))
					mdf <- c(length(mdif[!is.na(mdif)]),mean(mdif,na.rm=TRUE),.sd(mdif),min(mdif,na.rm=TRUE),max(mdif,na.rm=TRUE))
					
					##   Compile the descriptives for each model and print them
					if (object@model[i]=="mcm") {
						if (object@dimensions==1) {
							des <- cbind(a,b,c)
							colnames(des) <- c("a","b","c")
						} else {
							des <- cbind(a,b,c,mdc,mdf)
							colnames(des) <- c(paste("a",1:object@dimensions,sep=""),"d","c","MDISC","MDIF")
						}
					} else {
						if (object@dimensions==1) {
							des <- cbind(a,b)
							colnames(des) <- c("a","b")
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
					print(format(as.data.frame(des),nsmall=4),quote=FALSE)
					cat("\n")
				} else if (object@n[j]==1) { 
					cat("There is only one item. No descriptives available.\n\n")
				}
			}
		}
	}
}



##   Summarize the linking constants and common item parameters
summary.link <- function(object, ..., descrip=FALSE) {

	##   Create a matrix containing the shortened name
	##   and the descrtiptive name for all the linking methods
	##   (both unidimensional and multidimensional)
	nms <- matrix(c("MM","MS","HB","SL","LS","Mean/Mean","Mean/Sigma","Haebara","Stocking-Lord","Least Squares"),5,2)
	
	##   Initialize an object to hold the names for all used methods
	tmp.nms <- NULL
	
	##   Loop through the object containing the linking constants
	##   and identify the names of all used methods
	for (j in 1:length(names(object@constants))) {
		tmp.nms <- c(tmp.nms, nms[nms[,1]==names(object@constants)[j],2])
	}
	
	##   Print a header for each group
	cat("------- ",object@grp.names," -------\n")
	cat("Linking Constants\n\n")
	
	##   Determine if the constants correspond to a unidmensional
	##   versus multidimensional method (the linking constants
	##   are stored as a list in the multidimensional case)
	if (is.list(object@constants[[1]])) {
	
		##   Loop through all of the methods
		for (j in 1:length(tmp.nms)) {
		
			##   Print the method name
			cat(tmp.nms[j],"\n")
			
			##   Print the linking constants for the oblique procrustes method
			if (length(object@constants[[j]] )==2) { 
				print(format(object@constants[[j]]$A,digits=6,drop0trailing=FALSE),nsmall=6,quote=FALSE)
				cat("\n")
				print(format(object@constants[[j]]$m,digits=6,drop0trailing=FALSE),nsmall=6,quote=FALSE)
				cat("\n")
			
			##   Print the linking constants for the orthogonal procrustes methods
			} else {
				tmp <- object@constants[[j]]
				tmp1 <- format(cbind(tmp$T,NA,tmp$K),digits=6,drop0trailing=FALSE)
				tmp1[,ncol(tmp$T)+1] <- NA
				print(tmp1,na.print="",quote=F)
				cat("\n")
				print(format(tmp$m,digits=6,drop0trailing=FALSE),nsmall=6,quote=FALSE)
				cat("\n")
			}
		}
		
	##   This is for the unidimensional case
	} else {
		tmp <- NULL
		
		##   Combine the constants for all methods into a matrix
		for (j in 1:length(tmp.nms)) {
			tmp <- rbind(tmp,object@constants[[j]])
		}
		rownames(tmp) <- tmp.nms
		print(format(data.frame(tmp),digits=6,nsmall=6, drop0trailing=FALSE),quote=FALSE)
	}
	cat("\n")
	
	##   Print the descriptive statistics for the common item parameters
	if (descrip==TRUE) {
		cat("Common Item Descriptive Statistics\n\n")
		
		##   Determine if there are multiple models
		if (length(object@n)>3) n <- object@n[-3] else n <- object@n
		if (length(n[n>0])==2) {
			n <- n[-1][n[-1]>0]
		} else {
			n <- n[c(2:length(n),1)]
		}
		
		##   Identify all of the item response models
		##   used for the common items
		ml <- object@mod.lab
		if (length(ml)>1) ml <- c(ml, "All")
		
		##   Loop through all of the item response models and 
		##   print the descriptive statistics for the item parameters
		for (j in 1:length(n[n>0])) {
			cat("Model:",ml[j],"\n")
			cat("Number of Items:",n[j],"\n\n")
			print(format(as.data.frame(object@descriptives[[j]]),nsmall=6),quote=FALSE)
			cat("\n")
		}
		cat("\n")
	}
}