link.con <- function(x) {
	if (is.list(x)) {
		if (is.list(x[[1]])) {
			if (is.link(x[[1]][[1]])) {
				out <- vector("list",length(x[[1]]))
				for (i in 1:length(x[[1]])) {
					obj <- x[[1]][[i]]
					if (is.null(obj@SL)) {
						tmpa <- rbind(obj@MM,obj@MS,obj@HB)
						rownames(tmpa) <- c("MM","MS","HB")
					} else {
						tmpa <- rbind(obj@MM,obj@MS,obj@HB,obj@SL)
						rownames(tmpa) <- c("MM","MS","HB","SL")
					}
					out[[i]] <- tmpa
					names(out)[[i]] <- names(x[[1]])[i]
				}
			} else {
				stop("The objects in the list are not of class {link}")
			}
		}  else {
			if (is.link(x[[1]]) & !is.link(x[[2]])) {
				obj <- x[[1]]
				if (is.null(obj@SL)) {
					out <- rbind(obj@MM,obj@MS,obj@HB)
					rownames(out) <- c("MM","MS","HB")
				} else {
					out <- rbind(obj@MM,obj@MS,obj@HB,obj@SL)
					rownames(out) <- c("MM","MS","HB","SL")
				}
				if (x[[1]]@base.grp==1) {
					cat("Base Group: 1\n")
				} else {
					cat("Base Group: 2\n")
				}
			} else if (is.link(x[[1]]) & is.link(x[[2]])) {
				out <- vector("list", length(x))
				for (i in 1:length(x)) {
					obj <- x[[i]]
					if (is.null(obj@SL)) {
						tmpa <- rbind(obj@MM,obj@MS,obj@HB)
						rownames(tmpa) <- c("MM","MS","HB")
					} else {
						tmpa <- rbind(obj@MM,obj@MS,obj@HB,obj@SL)
						rownames(tmpa) <- c("MM","MS","HB","SL")
					}
					out[[i]] <- tmpa
					names(out)[[i]] <- names(x)[i]
				}
			} else {
				stop("The objects in the list are not of class {link}")
			}
		}
	} else if (is.link(x)) {
		if (is.null(x@SL)) {
			out <- rbind(x@MM,x@MS,x@HB)
			rownames(out) <- c("MM","MS","HB")
		} else {
			out <- rbind(x@MM,x@MS,x@HB,x@SL)
			rownames(out) <- c("MM","MS","HB","SL")
		}
		if (x@base.grp==1) {
			cat("Base Group: 1\n")
		} else {
			cat("Base Group: 2\n")
		}
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
