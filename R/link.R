link.con <- function(x, method="ALL") {
	nms <- matrix(c("MM","MS","HB","SL","RM","Mean/Mean","Mean/Sigma","Haebara","Stocking-Lord","Reckase-Martineau"),5,2)
	
	link.con.loop <- function(obj, method) {
		if (is.numeric(obj[[1]])) dimensions <- 1 else dimensions <- nrow(obj[[1]]$A)
		cons <- NULL
		cons.nms <- NULL
		if (dimensions==1) {
			if (toupper(method[1])=="ALL") {
				for (j in 1:length(obj)) {
					cons <- rbind(cons,obj[[j]])
					cons.nms <- c(cons.nms, nms[nms[,1]==names(obj[j]),2])
				}
			} else {
				for (j in 1:length(method)) {
					cons <- rbind(cons,eval(parse(text=paste("obj$",method[j],sep=""))))
					cons.nms <- c(cons.nms, nms[nms[,1]==method[j],2])
				}
			}
			rownames(cons) <- cons.nms
		} else {
			cons <- obj
			cons$HB$T <- NULL
			cons$HB$K <- NULL
			cons$SL$T <- NULL
			cons$SL$K <- NULL
			if (("ALL" %in% toupper(method))==FALSE) {
				if (("RM" %in% toupper(method))==FALSE) cons$RM <- NULL
				if (("HB" %in% toupper(method))==FALSE) cons$HB <- NULL
				if (("SL" %in% toupper(method))==FALSE) cons$SL <- NULL
			}
			if (length(cons)==1) cons <- cons[[1]]
		}
		return(cons)
	}
	
	if (is.list(x)) {
		if (is.list(x[[1]])) {
			if (is.link(x[[1]][[1]])) {
				# This scenario is when there are more than two groups and rescaled item parameters are returned
				out <- vector("list",length(x[[1]]))
					for (i in 1:length(x[[1]])) {
						out[[i]] <- link.con.loop(x[[1]][[i]]@constants, method)
						names(out)[[i]] <- x[[1]][[i]]@grp.names
					}
			} else {
				stop("The objects in the list are not of class {link}")
			}
		}  else {
			if (is.link(x[[1]]) & !is.link(x[[2]])) {
				out <- link.con.loop(x[[1]]@constants, method)
			} else if (is.link(x[[1]]) & is.link(x[[2]])) {
				out <- vector("list", length(x))
				for (i in 1:length(x)) {
					out[[i]] <- link.con.loop(x[[i]]@constants, method)
					names(out)[[i]] <- x[[i]]@grp.names
				}
			} else {
				stop("The objects in the list are not of class {link}")
			}
		}
	} else if (is.link(x)) {
		out <- link.con.loop(x@constants, method)
	}
return(out)
}

link.pars <- function(x) {
	if (length(x$pars)>0) {
		return(x$pars@pars)
	} else {
		cat("There are no rescaled item parameters present. Re-run plink and specify an argument for {rescale}\n")
	}
}

link.ability <- function(x) {
	if (length(x$ability)>0) {
		return(x$ability)
	} else {
		cat("There are no rescaled ability estimates present. Re-run plink and specify an object for the argument {ability}\n")
	}
}