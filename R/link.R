##   Extract linking constants from {plink} output
link.con <- function(x, method="ALL") {
	
	##   Create a matrix containing the shortened name
	##   and the descrtiptive name for all the linking methods
	##   (both unidimensional and multidimensional)
	nms <- matrix(c("MM","MS","HB","SL","LS","Mean/Mean","Mean/Sigma","Haebara","Stocking-Lord","Least Squares"),5,2)
	
	##   Function to 
	link.con.loop <- function(obj, method) {
	
		##   Determine if the constants correspond to a unidimensional or multidimensional model
		if (is.numeric(obj[[1]])) dimensions <- 1 else dimensions <- nrow(obj[[1]][[1]])
		
		##   Initialize an object to store the linking constants
		cons <- NULL
		
		##   Initialize an object to store the names of the used linking methods
		cons.nms <- NULL
		
		if (dimensions==1) {
			if (toupper(method[1])=="ALL") {
			
				##   Combine all of the linking constants into a matrix
				for (j in 1:length(obj)) {
					cons <- rbind(cons,obj[[j]])
					cons.nms <- c(cons.nms, nms[nms[,1]==names(obj[j]),2])
				}
			} else {
			
				##   Combine the linking constants for the specified models into a matrix
				for (j in 1:length(method)) {
					cons <- rbind(cons,eval(parse(text=paste("obj$",method[j],sep=""))))
					cons.nms <- c(cons.nms, nms[nms[,1]==method[j],2])
				}
			}
			rownames(cons) <- cons.nms
		
		##   For multiple dimensions
		} else {
		
			##   Extract the list of the linking matrices/vectors
			cons <- obj
			
			##   Remove the constants for certain methods
			if (("ALL" %in% toupper(method))==FALSE) {
				if (("LS" %in% toupper(method))==FALSE) cons$LS <- NULL
				if (("HB" %in% toupper(method))==FALSE) cons$HB <- NULL
				if (("SL" %in% toupper(method))==FALSE) cons$SL <- NULL
			} 
			if (length(cons)==1) cons <- cons[[1]]
		}
		return(cons)
	}
	
	##   This will occur when there are more than two groups and/or when
	##   rescaled item and/or ability parameters are returned
	if (is.list(x)) {
		
		##   This will occur when there are more than two groups
		if (is.list(x[[1]])) {
			
			##   This scenario is when there are more than two groups and rescaled 
			##   item and/or ability parameters are returned
			if (is.link(x[[1]][[1]])) {
			
				##   Initialize an object to store the extracted linking constants
				out <- vector("list",length(x[[1]]))
				
				##   Loop through all of the sets of linking output and extract the specified constants
				for (i in 1:length(x[[1]])) {
					out[[i]] <- link.con.loop(x[[1]][[i]]@constants, method)
					names(out)[[i]] <- x[[1]][[i]]@grp.names
				}
			} else {
				stop("The objects in the list are not of class {link}")
			}
		
		##   If the first element is not a list
		}  else {
			##   This will occur when there were only two groups and 
			##   rescaled item and/or ability parameters were returned
			if (is.link(x[[1]]) & !is.link(x[[2]])) {
				out <- link.con.loop(x[[1]]@constants, method)
				
			##   This will occur when there are more than two groups and 
			##   rescaled item and/or ability parameters were not returned
			} else if (is.link(x[[1]]) & is.link(x[[2]])) {
				out <- vector("list", length(x))
				
				##   Loop through all of the sets of linking output and extract the specified constants
				for (i in 1:length(x)) {
					out[[i]] <- link.con.loop(x[[i]]@constants, method)
					names(out)[[i]] <- x[[i]]@grp.names
				}
			} else {
				stop("The objects in the list are not of class {link}")
			}
		}
		
	##   This is the scenario where there were only two groups
	##   and no rescaled parameters were returned
	} else if (is.link(x)) {
		out <- link.con.loop(x@constants, method)
	}
return(out)
}



##   Extract rescaled item parameters from {plink} output
link.pars <- function(x) {

	if (length(x$pars)>0) {
		return(x$pars@pars)
	} else {
		cat("There are no rescaled item parameters present. Re-run plink and specify an argument for {rescale}\n")
	}
}



##   Extract rescaled ability parameters from {plink} output
link.ability <- function(x) {

	if (length(x$ability)>0) {
		return(x$ability)
	} else {
		cat("There are no rescaled ability estimates present. Re-run plink and specify an object for the argument {ability}\n")
	}
}



##   Extract expected probabilities from {irt.prob} output
get.prob <- function(x) {

	if (class(x)=="irt.prob") {
		##   Extract the expected probabilities
		out <- x@prob
		
	} else if (class(x)=="list") {
		if (class(x[[1]])=="irt.prob") {
			out <- vector("list",length(x))
			for (i in 1:length(x)) {
				out[[i]] <- x[[i]]@prob
			}
			names(out) <- names(x)
		} else {
			stop("{x} does not contain an object of class irt.prob.")
		} 
	} else {
		stop("{x} does not contain an object of class irt.prob.")
	}
	return(out)
}