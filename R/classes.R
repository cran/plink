# poly.mod
setClass("poly.mod",representation(model="character", items="list"))
.Valid.poly.mod <- function(object) {
	out <- TRUE
	if (length(object@model)!=length(object@items)) out <- "Number of models does not match the number of model items."
	if (length(which(object@model %in% c("drm","grm","gpcm","nrm","mcm") ))!=length(object@model)) {
		out <- "One or more of the specified models is not {drm, gpcm, grm, nrm, mcm}"
	}
	return(out)
}
setValidity("poly.mod", .Valid.poly.mod)


as.poly.mod <- function(n, model="drm", items=NULL) {
	if (missing(items)) {
		if (length(model)==1) items <- list(1:n )
	} else if (!is.list(items)) {
		items <- list(items)
	}
	if (length(unlist(items))!=n) stop("Number of items in {items} does not match {n}")
	names(items) <- model
	new("poly.mod", model=model, items=items)
}


is.poly.mod <- function(x) {
	is(x, "poly.mod")
}


# sep.pars
setClass("sep.pars", representation(a="matrix", b="matrix", c="matrix", cat="numeric", n="numeric", mod.lab="character", loc.out="logical"), prototype(loc.out=FALSE), contains="poly.mod")
.Valid.sep.pars <- function(object) {
	if (length(object@mod.lab)==length(object@model)) TRUE else paste("Number of model labels does not match the number of models.")
}
setValidity("sep.pars", .Valid.sep.pars)

is.sep.pars <- function(x) {
	is(x, "sep.pars")
}


# irt.prob
setClass("irt.prob",representation(prob="data.frame",p.cat="numeric", mod.lab="character"), contains="poly.mod")
.Valid.irt.prob<- function(object) {
	if (length(object@mod.lab)==length(object@model)) TRUE else paste("Number of model labels does not match the number of models.")
}
setValidity("irt.prob", .Valid.irt.prob)

is.irt.prob <- function(x) {
	is(x, "irt.prob")
}

#irt.pars
setClassUnion("list.num",c("list","numeric"))
setClassUnion("list.mat",c("list","matrix","data.frame","numeric","NULL"))
setClassUnion("list.poly",c("list","poly.mod"))
setClass("irt.pars",representation(pars="list.mat", cat="list.num", poly.mod="list.poly", common="list.mat", location="logical", groups="numeric"))

setMethod("initialize", "irt.pars", function(.Object, pars, cat, poly.mod, common=NULL, location=NULL, groups=1) {
	
	# Reformat pars to be a mtrix or a list of matrices
	if (is.data.frame(pars)|is.numeric(pars)) pars <- as.matrix(pars)
	if (is.list(pars)) {
		n <- length(pars) 
		if (n==1) {
			pars <- pars[[1]]
		} else if (n>1) {
			for (i in 1:n) {
				pars[[i]] <- as.matrix(pars[[i]])
			}
		}
	} else {
		 n <- 1
	}
	
	# Peform checks on cat
	if (n==1) {
		if (is.list(cat)) {
			if (length(cat)==1) {
				cat <- cat[[1]]
			} else {
				stop("There is more than one list element for {cat}. {cat} should be a vector")
			}
		}
	} else if (n>1) {
		if (is.list(cat)) {
			if (length(cat)!=n) {
				stop(paste("The number of category elements",length(cat),"does equal the number of item parameter sets", n))
			}
			tmp <- NULL
			for (i in 1:length(cat)) {
				if (!is.numeric(cat[[i]])) tmp <- c(tmp,i)
			}
			if (length(tmp)) {
				if (length(tmp==1)) {
					stop(paste("Category",tmp,"is not an object of class 'numeric'"))
				} else {
					stop(paste("Categories",tmp,"are not objects of class 'numeric'"))
				}
			}
		} else {
			stop("There are two or more sets of item parameters. {cat} should be a list.")
		} 
	}
	
	# Check poly.mod
	if (n==1) {
		if (is.list(poly.mod)) {
			if (length(poly.mod)==1) {
				poly.mod <- poly.mod[[1]]
			} else {
				stop("There is more than one list element for {poly.mod}.")
			}
		}
	} else if (n>1) {
		if (is.list(poly.mod)) {
			if (length(poly.mod)==1) {
				stop("There are two or more sets of item parameters. {poly.mod} should be a list.")
			} else {
				if (length(poly.mod)!=n) {
					stop(paste("The number of poly.mod list elements",length(poly.mod),"does equal the number of item parameter sets", n))
				}
				tmp <- NULL
				for (i in 1:length(poly.mod)) {
					if (!is.poly.mod(poly.mod[[i]])) tmp <- c(tmp,i)
				}
				if (length(tmp)) {
					if (length(tmp==1)) {
						stop(paste("{poly.mod} list element",tmp,"is not an object of class 'poly.mod'"))
					} else {
						stop(paste("{poly.mod} list elements",tmp,"are not objects of class 'poly.mod'"))
					}
				}
			}
		}
	}
	
	for (i in 1:n) {
		if (n==1) {
			if (length(cat)!=nrow(pars)) {
				stop("The number of items does not match the number of categories")
			}
			if(length(unlist(poly.mod@items))!=nrow(pars)) {
				stop("The number of items in {poly.mod} does not match number of items in {pars}")
			}
		} else if (n>1) {
			if (length(cat[[i]])!=nrow(pars[[i]])) {
				stop(paste("The number of items in {pars} element",i,"do not match the number of categories in {cat} element",i))
			}
			if(length(unlist(poly.mod[[i]]@items))!=nrow(pars[[i]])) {
				stop(paste("The number of items in {poly.mod} element",i,"do not match number of items in {pars} element",i))
			}
		}
	}
	
	if (is.null(location)) {
		location <- rep(FALSE,n)
	} else {
		if (length(location)!=n) {
			stop("The length of {location} should equal the number of groups")
		}
	}
	
		
	if(is.null(common)) {
		if (n>1) {
			stop("There are two or more sets of item parameters. {common} should be included.")
		}
	} else if (is.list(common)) {
		tmp <- NULL
		for (i in 1:length(common)) {
			if (!is.matrix(common[[i]])) tmp <- c(tmp,i)
		}
		if (length(tmp)) {
			if (length(tmp==1)) {
				stop(paste("{common} element",tmp,"is not an object of class 'matrix'"))
			} else {
				stop(paste("{common} elements",tmp,"are not objects of class 'matrix'"))
			}
		}
		if (length(common)!=(n-1)) {
			stop(paste("The number of common item sets",length(common),"does correspond to the number of item parameter sets", n))
		} else {
			if (n==2) common <- common[[1]]
		}
	}else if (is.data.frame(common)) {
		if (n==1) {
			common <- as.matrix(common)
		} else if (n>2) {
			stop("There are more than two sets of item parameters. {common} should be a list.")
		}
	} else if (is.matrix(common)) {
		if (n>2) {
			stop("There are more than two sets of item parameters. {common} should be a list.")
		}
	} else {
		if (n==1) {
			stop("{common} must be a matrix")
		} else if (n>1) {
			stop("There are more than two sets of item parameters. {common} should be a list.")
		}
	}
		
	.Object@pars <- pars
	.Object@cat <- cat
	.Object@poly.mod <- poly.mod
	.Object@common <- common
	.Object@location <- location
	.Object@groups<- n
	.Object
})

is.irt.pars <- function(x) {
	is(x, "irt.pars")
}


#link
setClassUnion("num.null",c("matrix","numeric","NULL"))
setClassUnion("list.dat",c("list","data.frame"))
setClass("link", representation(MM="numeric", MS="numeric", HB="numeric", SL="num.null", descriptives="list.dat", iterations="numeric", convergence="character", base.grp="numeric", include.mcm.nrm="logical") )

is.link <- function(x) {
	is(x, "link")
}
