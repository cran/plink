##   {poly.mod} class
##   This is a class that characterizes the item response models associated with each item
##   It includes a validity check between the number of models in {model} and the number
##   of list elements for {items}.  It also checks to make sure only acceptable values for
##   {model} are included

setClass("poly.mod",representation(model="character", items="list"))
.Valid.poly.mod <- function(object) {
	out <- TRUE
	if (length(object@model)!=length(object@items)) {
		out <- "Number of models does not match the number list elements for {items}."
	}
	if (length(which(object@model %in% c("drm","grm","gpcm","nrm","mcm") ))!=length(object@model)) {
		out <- "One or more of the specified models is not {drm, grm, gpcm, nrm, mcm}"
	}
	return(out)
}
setValidity("poly.mod", .Valid.poly.mod)



##   as.poly.mod function
##   This function facilitates the creation of a {poly.mod} object

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



##   is.poly.mod function
##   This function checks to see if {x} is an object of class {poly.mod}

is.poly.mod <- function(x) {
	is(x, "poly.mod")
}



##   {sep.pars} class
##   This class stores a set of separated item parameters, the number of response categories,
##   the total number of items (as well as the number of items associated with each item
##   response model), descriptive labels for the models (e.g., 3PL and Graded Response Model),
##   an indicator for whether various models include a location parameter, and the number of
##   dimensions. This class also extends the {poly.mod} class. It includes validity checks to
##   ensure that the number number of rows in {a}, {b}, and {c} and the length of {cat} all
##   equal {n}. It also checks to make sure that the number of model labels is the same as the 
##   number of models (from the {poly.mod} class).

setClass("sep.pars", representation(a="matrix", b="matrix", c="matrix", cat="numeric", n="numeric", mod.lab="character", location="logical", dimensions="numeric"), prototype(loc.out=FALSE, dimensions=1), contains="poly.mod")
.Valid.sep.pars <- function(object) {
	out <- TRUE
	if (nrow(object@a)!=object@n[1])  paste("The number of a-parameters does not match the total number of items.")
	if (nrow(object@b)!=object@n[1])  paste("The number of b-parameters does not match the total number of items.")
	if (nrow(object@b)!=object@n[1])  paste("The number of c-parameters does not match the total number of items.")
	if (length(object@cat)!=object@n[1])  paste("The length of {cat} does not match the total number of items.")
	if (length(object@mod.lab)!=length(object@model)) paste("Number of model labels does not match the number of models.")
	return(out)
}
setValidity("sep.pars", .Valid.sep.pars)



##   is.sep.pars function
##   This function checks to see if {x} is an object of class {sep.pars}

is.sep.pars <- function(x) {
	is(x, "sep.pars")
}



##   {irt.prob} class
##   This class stores a T x Nk data.frame of response probabilities for T theta values 
##   (or a combination of theta values if there is more than one dimension), a vector identifying
##   the number of columns in the data.frame associated with each item, descriptive labels for
##   the modeled items, the number of dimensions, and a list of item parameters. This class also
##   extends the {poly.mod} class.

setClass("irt.prob",representation(prob="data.frame", p.cat="numeric", mod.lab="character", dimensions="numeric", D="numeric", pars="list"), contains="poly.mod")
.Valid.irt.prob<- function(object) {
	
	out <- TRUE
	if ((ncol(object@prob)-object@dimensions)!=sum(object@p.cat)) paste("The number of columns in {prob} does not align with {p.cat} and {dimensions}")
	if (length(object@mod.lab)!=length(object@model)) paste("Number of model labels does not match the number of models.")
}
setValidity("irt.prob", .Valid.irt.prob)



##   is.sep.pars function
##   This function checks to see if {x} is an object of class {sep.pars}

is.irt.prob <- function(x) {
	is(x, "irt.prob")
}


##   {irt.pars} class
##   This class stores information for one or more sets of item parameters. It uses
##   three union classes {list.num}, {list.mat}, and {list.poly}. It includes slots for
##   item parameters, numbers of response categores, {poly.mod} objects, common
##   item matrices, whether location parameters are used, the number of groups, and
##   numbers of dimensions. The item parameters are stored as a matrix when there
##   is one group or as a list of matrices when there are multiple groups. Similarly,
##   the slots {cat}, {poly.mod}, and {common} are stored as a vector, poly.mod object,
##   or matrix respectively when there is one group or as a list with these formats
##   when there is more than one group.

setClassUnion("list.num",c("list","numeric"))
setClassUnion("list.mat",c("list","matrix","data.frame","numeric","NULL"))
setClassUnion("list.poly",c("list","poly.mod"))
setClass("irt.pars",representation(pars="list.mat", cat="list.num", poly.mod="list.poly", common="list.mat", location="logical", groups="numeric", dimensions="numeric"))

setMethod("initialize", "irt.pars", function(.Object, pars, cat, poly.mod, common=NULL, location=NULL, groups=1, dimensions=1) {
	
	##   Reformat {pars} to be a matrix or a list of matrices
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
	
	##   Check the {cat} object
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
			##   Check to see if each of the elements in {cat} is a numeric vector
			##   If not numeric, add the list element number to a temporary object so that
			##   an error message can be displayed
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
	
	##   Check the {poly.mod} object
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
				##   Check to see if each of the elements in {poly.mod} is of class {poly.mod}
				##   If not numeric, add the list element number to a temporary object so that
				##   an error message can be displayed
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
	
	##   Check the correspondence between {pars}/{cat} and {pars}/{poly.mod}
	##   That is, do the number of response categories match the number of items in {pars}
	##   and do the number of items in the {poly.mod} object(s) match the number
	##   of items in pars.
	for (i in 1:n) {
		if (n==1) {
			if (length(cat)!=nrow(pars)) {
				stop("The length of {cat} does not match the number of items")
			}
			if(length(unlist(poly.mod@items))!=nrow(pars)) {
				stop("The number of items in {poly.mod} does not match number of items in {pars}")
			}
		} else if (n>1) {
			if (length(cat[[i]])!=nrow(pars[[i]])) {
				stop(paste("The length of {cat} in element", i, "does not match the number of items in {pars} element",i))
			}
			if(length(unlist(poly.mod[[i]]@items))!=nrow(pars[[i]])) {
				stop(paste("The number of items in {poly.mod} element",i,"does not match number of items in {pars} element",i))
			}
		}
	}
	
	##   Check the {location} object
	if (is.null(location)) {
		location <- rep(FALSE,n)
	} else {
		if (length(location)!=n) {
			stop("The length of {location} should equal the number of groups")
		}
	}
	
	##   Check the {common} object
	if(is.null(common)) {
		if (n>1) {
			stop("There are two or more sets of item parameters. {common} should be included.")
		}
	} else if (is.list(common)) {
		##   Check to see if each of the elements in {common} is a matrix
		##   If not numeric, add the list element number to a temporary object so that
		##   an error message can be displayed
		tmp <- NULL
		for (i in 1:length(common)) {
			if (!is.matrix(common[[i]])) {
				if (is.data.frame(common[[i]])) {
					common[[i]] <- as.matrix(common[[i]])
				} else {
					tmp <- c(tmp,i)
				}
			}
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
	
	##   Check the {dimensions} object
	if (length(dimensions)!=groups) {
		stop("The number of elements in the vector {dimensions} must be equal to the number of groups.")
	}
	
	
	##   Check the correspondence between the number of response categories in {cat}
	##   and the response model associated with each item in {poly.mod}. 
	for (i in 1:n) {
		if (n==1) {
			tmp.cat <- list(cat) 
			pm <- list(poly.mod)
		} else {
			tmp.cat <- cat
			pm <- poly.mod
		}
		##   Loop through the response models
		for (j in 1:length(pm[[i]]@model)) {
			if (pm[[i]]@model[j]=="drm") {
				dichot <- tmp.cat[[i]][pm[[i]]@items$drm]
				if (length(dichot[dichot==2])!=length(pm[[i]]@items$drm)) {
					stop(paste("One or more values in {cat} for group",i,"does not align with the DRM items in the {poly.mod} object"))
				}
			} else {
				poly <- tmp.cat[[i]][pm[[i]]@items[[j]]]
				if (length(poly[poly>2])!=length(pm[[i]]@items[[j]])) {
					stop(paste("One or more values in {cat} for group",i,"are less than 2 for the", toupper(pm[[i]]@model[[j]]),"items in the {poly.mod} object"))
				}
			}
		}
	}
	
	
	##   Check the correspondence between the number of response categories in {cat}
	##   and the number of item parameters in {pars}
	for (i in 1:n) {
		if (n==1) {
			tmp.pars <- list(pars)
			tmp.cat <- list(cat) 
			pm <- list(poly.mod)
		} else {
			tmp.pars <- pars
			tmp.cat <- cat
			pm <- poly.mod
		}
		##   Loop through the response models
		for (j in 1:length(pm[[i]]@model)) {
		
			##   Number of items associated with the given model
			np <- length(pm[[i]]@items[[j]])
			
			if (pm[[i]]@model[j]=="drm") {
				tmp <- !is.na(tmp.pars[[i]][pm[[i]]@items$drm,])
				if (np==1) tmp <- t(tmp)
				dichot.p <- apply(tmp,1,sum)
				
				if (length(dichot.p[dichot.p<=(dimensions[i]+2)])!=length(pm[[i]]@items$drm)) {
					stop(paste("One or more NRM items in group",i,"as specified by the {poly.mod} object do not have the correct number of parameters"))
				}
			} else if (pm[[i]]@model[j]=="grm") {
				poly.c <- tmp.cat[[i]][pm[[i]]@items$grm]
				
				tmp <- !is.na(tmp.pars[[i]][pm[[i]]@items$grm,])
				if (np==1) tmp <- t(tmp)
				poly.p <- apply(tmp,1,sum)
				
				if (location[i]==TRUE) {
					if (length(poly.p[poly.p==(dimensions[i]+poly.c)])!=length(pm[[i]]@items$grm)) {
						stop(paste("One or more GRM items in group",i,"as specified by the {poly.mod} object do not have the correct number of parameters"))
					}
				} else {
					if (length(poly.p[poly.p==(dimensions[i]+poly.c-1)])!=length(pm[[i]]@items$grm)) {
						stop(paste("One or more GRM items in group",i,"as specified by the {poly.mod} object do not have the correct number of parameters"))
					}
				}
			} else if (pm[[i]]@model[j]=="gpcm") {
				poly.c <- tmp.cat[[i]][pm[[i]]@items$gpcm]
				
				tmp <- !is.na(tmp.pars[[i]][pm[[i]]@items$gpcm,])
				if (np==1) tmp <- t(tmp)
				poly.p <- apply(tmp,1,sum)
	
				if (location[i]==TRUE) {
					if (length(poly.p[poly.p==(dimensions[i]+poly.c)])!=length(pm[[i]]@items$gpcm)) {
						stop(paste("One or more GPCM items in group",i,"as specified by the {poly.mod} object do not have the correct number of parameters"))
					}
				} else {
					if (length(poly.p[poly.p==(dimensions[i]+poly.c-1)])!=length(pm[[i]]@items$gpcm)) {
						stop(paste("One or more GPCM items in group",i,"as specified by the {poly.mod} object do not have the correct number of parameters"))
					}
				}
			} else if (pm[[i]]@model[j]=="nrm") {
				poly.c <- tmp.cat[[i]][pm[[i]]@items$nrm]
				
				tmp <- !is.na(tmp.pars[[i]][pm[[i]]@items$nrm,])
				if (np==1) tmp <- t(tmp)
				poly.p <- apply(tmp,1,sum)
				
				if (length(poly.p[poly.p==dimensions[i]*poly.c*2])!=length(pm[[i]]@items$nrm)) {
					stop(paste("One or more NRM items in group",i,"as specified by the {poly.mod} object do not have the correct number of parameters"))
				}
			} else if (pm[[i]]@model[j]=="mcm") {
				poly.c <- tmp.cat[[i]][pm[[i]]@items$mcm]
				
				tmp <- !is.na(tmp.pars[[i]][pm[[i]]@items$mcm,])
				if (np==1) tmp <- t(tmp)
				poly.p <- apply(tmp,1,sum)
				
				if (length(poly.p[poly.p==(dimensions[i]*poly.c*2+(poly.c-1))])!=length(pm[[i]]@items$mcm)) {
					stop(paste("One or more MCM items in group",i,"as specified by the {poly.mod} object do not have the correct number of parameters"))
				}
			}
			
		}
	}
	
	##   Check to make sure the common items in each adjacent groups have the same number of parameters.
	##   Because the correspondence between {cat} and {pars} has already been checked {cat} is used
	##   to make the comparison
	if (n>1) {
		if (n==2) com <- list(common) else com <- common
		for (i in 1:(n-1)) {
			cat1 <- cat[[i]][com[[i]][,1]]
			cat2 <- cat[[i+1]][com[[i]][,2]]
			if (sum(cat1-cat2)!=0) {
				stop(paste("One or more common items in groups ",i,"/",i+1," do not have the same number of parameters",sep="")) 
			}
		}
	}
	
	.Object@pars <- pars
	.Object@cat <- cat
	.Object@poly.mod <- poly.mod
	.Object@common <- common
	.Object@location <- location
	.Object@groups<- n
	.Object@dimensions <- dimensions
	.Object
})



##   is.poly.mod function
##   This function checks to see if {x} is an object of class {irt.pars}

is.irt.pars <- function(x) {
	is(x, "irt.pars")
}



##   {link} class
##   This class stores information from the estimation of linking constants in {plink}. 
##   It uses four union classes {num.null}, {list.dat}, {list.null}, {pars.null}. 
##   It also includes slots for rescaled item and ability parameters

setClassUnion("num.null",c("matrix","numeric","NULL"))
setClassUnion("list.dat",c("list","data.frame"))
setClassUnion("list.null",c("list","NULL"))
setClassUnion("pars.null",c("irt.pars","NULL"))

setClass("link", representation(constants="list", descriptives="list.dat", iterations="numeric", objective="numeric", convergence="character", base.grp="numeric", n="numeric", grp.names="character", mod.lab="character", dilation="character") )

##   is.link function
##   This function checks to see if {x} is an object of class {link}

is.link <- function(x) {
	is(x, "link")
}