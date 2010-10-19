##   Import item parameters or ability estimates from BILOG-MG 3
read.bilog <- function(file, ability=FALSE, pars.only=TRUE, as.irt.pars=TRUE) {
	
	##   Import item parameters
	if (ability==FALSE) {
	
		##   Read in the data
		pars <- read.fwf(file, c(8,8,rep(10,13),4,1,1),sep="\t", skip=4)
		pars <- pars[,-15]
		
		colnames(pars) <- c("item", "subtest", "intercept","int.se", "slope", "slope.se", "threshold", "thresh.se", "dispersion", "disp.se", "asymptote", "asymp.se", "drift", "drift.se","stream.loc","key","dummy.vals")
		
		##   Eliminate all columns except those for the slope, difficulty, and lower asymptote parameters
		if (pars.only==TRUE) pars <- pars[,c(5,7,11)]
		
		##   Create an irt.pars object
		if (as.irt.pars==TRUE) {
			if (pars.only==FALSE) pars <- pars[,c(5,7,11)]
			n <- nrow(pars)
			pm <- as.poly.mod(n)
			pars <- as.irt.pars(pars, cat=rep(2,n), poly.mod=pm)
		}
	
	##   Import ability estimates
	} else {
		##   Read in the data
		pars <- read.fwf(file, c(3,13,4,5,10,12,12,11,10),sep="\t", skip=2)
		
		pars[,2] <- suppressWarnings(as.numeric(as.character(pars[,2])))
		pars[,1] <- c(0,pars[,1][-length(pars[,1])])
		pars[,2] <- c(0,pars[,2][-length(pars[,2])])
		pars <- pars[seq(2,nrow(pars),2),1:7]
		pars[,5] <- pars[,5]/100
		pars[pars[,7]==999,7] <- NA
		
		colnames(pars) <- c("group","id","n.pos","n.cor","p","theta","theta.se")
	}
	
	return(pars)
}



##   Import item parameters or ability estimates from PARSCALE 4
read.parscale <- function(file, ability=FALSE,  loc.out=FALSE, pars.only=TRUE, as.irt.pars=TRUE) {
	
	##   Import item parameters
	if (ability==FALSE) {
		
		##   Read the third line of the .PAR file which identifies the number
		##   of parameter blocks, the total number of items, and the model code
		nums <- as.numeric(read.fwf(file, c(8,rep(5,5)), skip=2, n=1)[-1])
		
		##   Read the fourth line in the .PAR file which identifies
		##   the number of items in each block
		block <- as.numeric(read.fwf(file, rep(5,nums[1]), skip=3, n=1))
		
		##   Initialize an object to store the item parameters
		pars <- NULL
		
		##   Initialize a list to store the deviation threshold/step parameters for each block
		pars1 <- vector("list",length(block))
		
		##   Initialize an object to store the number of response 
		##   categories for each item
		cat <- NULL
		
		##   Identify the starting row for reading in the
		##   parameters for a given block
		skip <- 5
		
		##   Identify the starting row for reading in the
		##   category parameters for a given block
		skip1 <- skip+block[1]
		
		##   Loop through each block of parameters
		for (i in 1:length(block)) {
			if (i>1) {
				##   Update the starting rows for the given block
				skip <- 5 + sum(block[1:(i-1)]) + (i-1)*2
				skip1 <- 5 + sum(block[1:i]) + (i-1)*2
			}
			
			##   Read the item parameters for the given block
			tmp <- read.fwf(file, c(8,5,4,rep(10,6)), skip=skip, n=block[i], colClasses=c("character","numeric","character",rep("numeric",6)))
			
			##   Add this block of parameters to the full set of parameters
			pars <- rbind(pars, tmp)
			
			##   If there is more than one item in the block
			if (is.data.frame(tmp)) {
				cat <- c(cat, as.numeric(tmp[,2]))
				
			##   If there is only one item in the block
			} else {
				cat <- c(cat, as.numeric(tmp[2]))
			}
			
			##   Read in the deviation threshold/step parameters parameters for the given block
			tmp1 <- read.fwf(file, rep(10,cat[length(cat)]), skip=skip1, n=2)
			
			##   For the graded response models
			##   Eliminate the last category parameter
			if (nums[3]<5) {
				tmp1 <- tmp1[,-ncol(tmp1)] 
			
			##   For the  partial credit models
			##   Eliminate the first category parameter
			} else {
				tmp1 <- tmp1[,-1]
			}
			pars1[[i]] <- tmp1
		}
		colnames(pars) <- c("block.name","cat","item","slope","slope.se","location","loc.se","asymptote","asymp.se")
		
		##   For polytomous items
		if (max(cat)>2) {
			##   Create a matrix for the threshold/step parameters
			step <- matrix(NA, nrow(pars), 2*(max(cat)-1))
			
			##   Because multiple items can all have the same
			##   deviation threshold/step parameters, expand {pars1}
			##   so that each item has a corresponding list element
			##   containing threshold/step parameters
			pars1 <- rep(pars1,block)
			
			##   Loop through all of the list elements for {pars1}
			##   and insert the deviation threshold/step parameters
			##   into the matrix {step}
			for (i in 1:length(pars1)) {
				##   We only need deviation threshold/step parameters
				##   for polytomous items. In the list, these will be
				##   formatted as a data.frame whereas the parameters
				##   for dichotomous items will be formatted as a 
				##   vector (of zeros)
				if (is.data.frame(pars1[[i]])) {
					step[i,] <- unlist(pars1[[i]])
				}
			}
			
			if (nums[3]<5) prefix <- "thresh" else prefix <- "step"
			
			##   Create column names for the matrix of deviation threshold/step parameters
			step.names <- character()
			for (j in 1:(max(cat)-1)) {
				step.names <- c(step.names, paste(prefix,j,sep=""), paste(prefix,j,".se",sep=""))
			}
			colnames(step) <- step.names
			
			##   Combine the matrix of item parameters with the matrix 
			##   of deviation threshold/step parameters
			pars <- cbind(pars, step)
		}
		pars <- data.frame(pars)
		
		##   For polytomous items, reformulate the parameters to
		##   exclude the location parameter (if applicable)
		if (max(cat)>2) {
			if (loc.out==TRUE) {
				if (nums[3] %in% c(2,4,6,8)) {
					pars[cat>2,6] <- mean(pars[cat>2,seq(10,ncol(pars),2)], na.rm=T)
					pars[cat>2,7] <- NA
					pars[cat>2,seq(10,ncol(pars),2)] <- pars[cat>2,seq(10,ncol(pars),2)]-pars[cat>2,6]
				}
			} else {
				if (nums[3] %in% c(1,3,5,7)) {
					pars[cat>2,seq(10,ncol(pars),2)] <- pars[cat>2,seq(10,ncol(pars),2)]+ pars[cat>2,6]
					pars[cat>2,6:7] <- 0
				} 
			}
		} else {
			loc.out <- FALSE
		}
		
		##   Get rid of all the columns of containing standard errors
		##   the item names, and the number of categories
		if (pars.only==TRUE) pars <- pars[,seq(4,ncol(pars),2)]
		
		##   Create an {irt.pars} object
		if (as.irt.pars==TRUE) {
		
			##   Get rid of all the columns of containing standard errors
			##   the item names, and the number of categories
			if (pars.only==FALSE) pars <- pars[,seq(4,ncol(pars),2)]
			
			##   Make sure that the matrix of parameters is formatted properly
			##   The threshold/step parameters in the matrix need to be shifted
			##   one or two columns to the left depending on if {loc.out} equals true
			if (max(cat)>2) {
				if (loc.out==TRUE) {
					pars[cat>2,3:(ncol(pars)-1)] <- pars[cat>2,4:ncol(pars)]
					pars <- pars[,-ncol(pars)]
				} else {
					pars[cat>2,2:(ncol(pars)-2)] <- pars[cat>2,4:ncol(pars)]
					pars <- pars[,1:(ncol(pars)-2)]
				}
			} else {
				loc.out <- FALSE
			}
			pars <- as.matrix(pars)
			colnames(pars) <- NULL
			
			##   Number of items
			n <- nrow(pars)
			
			##   Create the necessary {poly.mod} object
			if (min(cat)==2) {
				if (max(cat)==2) {
					mod <- "drm"
					items <- list(1:n)
				} else {
					if (nums[3]<5) {
						mod <- c("drm", "grm")
					} else {
						mod <- c("drm", "gpcm")
					}
					ni <- 1:n
					items <- list(ni[cat==2], ni[cat>2])
				}
			} else {
				if (nums[3]<5) {
					mod <- "grm"
				} else {
					mod <- "gpcm"
				}
				items <- 1:n
			}
			pm <- as.poly.mod(n, mod, items)
			pars <- as.irt.pars(pars, cat=cat, poly.mod=pm, location=loc.out)
		}
		
	##   Import ability estimates
	} else {
		pars <- read.fwf(file, c(21,32,12,12))
		id <- as.character(pars[seq(1,nrow(pars),2),1])
		ability <- pars[seq(2,nrow(pars),2),3:4]
		ability[ability==999] <- NA
		pars <- cbind(id,ability)
		colnames(pars) <- c("id","theta","theta.se")
	}
	
	return(pars)
}




##   Import item parameters or ability estimates from MULTILOG 7
read.multilog <- function(file, cat, poly.mod, ability=FALSE, contrast="dev", drm.3PL=TRUE, loc.out=FALSE, as.irt.pars=TRUE) {
	
	##   Import item parameters
	if (ability==FALSE) {
	
		##   Read in 30 columns
		##   This should be more columns than actually exist
		##   However, individuals can specify an output format
		##   (i.e., for the number of columns) in MULTILOG
		##   This simply reads in more columns and eliminates
		##   all columns with no data
		pars <- read.fwf(file, rep(12,30))
		
		##   Get rid of columns with no data
		tmp <- apply(is.na(pars),2,sum)
		pars <- pars[,tmp<nrow(pars)]
		
		##   The last row in the .PAR file should equal
		##   -1.00000     0.00000     1.00000
		##   Remove this row
		if (sum(abs(pars[nrow(pars),1:3])-c(1,0,1))==0) pars <- pars[-nrow(pars),]
		
		##   Prepare to extract the item parameters for each item
		mod <- poly.mod@model
		items <- poly.mod@items
		p.cat <- cat
		
		##   Loop through all of the item response models
		for (i in 1:length(mod)) {
			##   Identify the number of parameters that should 
			##   be extracted for each item
			if (mod[i]=="drm") {
				if (drm.3PL==TRUE) {
					p.cat[items[[i]]] <- 4
					
					## Check for 2PL items
					p.cat[items[[i]][is.na(pars[items[[i]],3])]] <- 2
				}
				
			} else if (mod[i] %in% c("gpcm","nrm","mcm")) {
				p.cat[items[[i]]] <- (p.cat[items[[i]]]-1)*3
			}
		}
		
		##   Reformat the read-in data as a vector
		##   and eliminate all missing values
		pars <- as.vector(t(as.matrix(pars)))
		pars <- pars[!is.na(pars)]
		
		##   Initialize a list to hold the parameters for each item
		p <- vector("list", length(cat))
		
		##   Loop through all the items and extract the 
		##   appropriate number of parameters
		k <- 1
		for (i in 1:length(cat)) {
			p[[i]] <- pars[k:(k-1+p.cat[i])]
			if (p.cat[i]==2) p[[i]] <- c(pars[k:(k-1+p.cat[i])], 0)
			k <- k+p.cat[i]
		}
		
		##   Determine the maximum number of parameters across
		##   all items, for each parameter type for use in compiling 
		##   the final matrix of item parameters
		a.max <- c.max <- 1
		d.max <- 0
		
		if (!is.null(poly.mod@items$drm)) {
			if (drm.3PL==TRUE) d.max <- 1
		}
		if (!is.null(poly.mod@items$grm)) {
			tmp <- max(cat[poly.mod@items$grm])
			if (tmp>c.max) {
				c.max <- tmp
				if (loc.out==TRUE) c.max <- c.max+1
			}
		}
		if (!is.null(poly.mod@items$gpcm)) {
			tmp <- max(cat[poly.mod@items$gpcm])
			if (tmp>c.max) {
				c.max <- tmp
				if (loc.out==TRUE) c.max <- c.max+1
			}
		}
		if (!is.null(poly.mod@items$nrm)) {
			tmp <- max(cat[poly.mod@items$nrm])
			if (tmp>a.max) a.max <- c.max <- tmp
		}
		if (!is.null(poly.mod@items$mcm)) {
			tmp <- max(cat[poly.mod@items$mcm])
			if (tmp>a.max) a.max <- c.max <- tmp
			d.max <- tmp-1
		}
		
		##   Identify the number of columns that should be
		##   included in the final matrix of item parameters
		col.max <- a.max+c.max+d.max
		
		##   MULTILOG implements deviation, polynomial, and triangular contrasts
		##   when estimating the item parameters for all but the 1PL, 2PL, and GRM
		##   These contrast coefficients need to be transformed to get the actual
		##   item parameters.   All of the contrasts used in this function come from  
		##   Du Toit, M. (Ed.) IRT from SSI. (2003) pp. 570-575
		
		##   Create a list polynomial contrasts
		poly <- list(
		##  Two categories
		c(-.5,.5),
		
		##  Three categories
		matrix(c(-1,0,1,
			.58,-1.15,.58), 2, 3, byrow=TRUE),
			
		##  Four categories
		matrix(c(-1.5, -.5, .5, 1.5,
			1.12, -1.12, -1.12, 1.12,
			-.5, 1.5, -1.5, .5), 3, 4, byrow=TRUE),
			
		##  Five categories
		matrix(c(-2, -1, 0, 1, 2,
			1.69, -.85, -1.69, -.85, 1.69,
			-1, 2, 0, -2, 1,
			.38, -1.51, 2.27, -1.51, .38), 4, 5, byrow=TRUE),
		
		##  Six categories
		matrix(c(-2.5, -1.5, -.5, .5, 1.5, 2.5,
			2.28, -.46, -1.83, -1.83, -.46, 2.28,
			-1.56, 2.18, 1.25, -1.25, -2.18, 1.56,
			.79, -2.37, 1.58, 1.58, -2.37, .79,
			-.26, 1.32, -2.64, 2.64, -1.32, .26), 5, 6, byrow=TRUE),
		
		##  Seven categories
		matrix(c(-3, -2, -1, 0, 1, 2, 3,
			2.89, 0, -1.73, -2.31, -1.73, 0, 2.89,
			-2.16, 2.16, 2.16, 0, -2.16, -2.16, 2.16,
			1.28, -2.98, 0.43, 2.56, .43, -2.98, 1.28,
			-.58, 2.31, -2.89, 0, 2.89, -2.31, .58,
			.17, -1.04, 2.61, -3.48, 2.61, -1.04, .17), 6, 7, byrow=TRUE),
		
		##  Eight categories
		matrix(c(-3.5, -2.5, -1.5, -.5, .5, 1.5, 2.5, 3.5,
			3.5, .5, -1.5, -2.5, -2.5, -1.5, .5, 3.5,
			-2.79, 1.99, 2.79, 1.2, -1.2, -2.79, -1.99, 2.79,
			1.83, -3.39, -.78, 2.35, 2.35, -.78, -3.39, 1.83,
			-.97, 3.19, -2.36, -2.08, 2.08, 2.36, -3.19, .97,
			.4, -1.99, 3.59, -1.99, -1.99, 3.59, -1.99, .4,
			-.11, .77, -2.32, 3.87, -3.87, 2.32, -.77, 0.11), 7, 8, byrow=TRUE),
		
		##  Nine categories
		matrix(c(-4, -3, -2, -1, 0, 1, 2, 3, 4,
			4.12, 1.03, -1.18, -2.5, -2.94, -2.5, -1.18, 1.03, 4.12,
			-3.45, 1.72, 3.2, 2.22, 0, -2.22, -3.2, -1.72, 3.45,
			2.42, -3.64, -1.9, 1.56, 3.12, 1.56, -1.9, -3.64, 2.42,
			-1.43, 3.94, -1.43, -3.22, 0, 3.22, 1.43, -3.94, 1.43,
			.7, -2.96, 3.83, .17, -3.48, .17, 3.83, -2.96, .7,
			-.26, 1.59, -3.7, 3.7, 0, -3.7, 3.7, -1.59, .26,
			.07, -.55, 1.91, 3.82, 4.78, -3.82, 1.91, -.55, .07), 8, 9, byrow=TRUE))
		
		##   Number of items
		n <- length(cat)
		
		##   Determine the contrasts to be used with each parameter for each item
		if (is.character(contrast)) {
			##   In the first formulation of the argument {contrast}, the user
			##   can supply a character value identifying that a specific contrast
			##   should be used for all parameters for all items. Create a list
			##   object that associates all parameters for all items with this type
			##   of contrast
			con <- list(dev=NULL,poly=NULL,tri=NULL)
			contrast <- tolower(contrast)
			eval(parse(text=paste("con$",contrast,"=1:",n,sep="")))
			con <- rep(con,3)
		} else {
			if (length(contrast)!=9) stop("The object {constant} must be a list of length nine")
			
			##   Initialize an object to store the final set of contrasts
			con <- contrast
			
			##   In the second  and third formulation, a list is supplied where
			##   the first three elements correspond to slope parameters for the 
			##   deviation, polynomial, and triangular contrasts respectively. The
			##   next three elements correspond to the category parameters for 
			##   the three contrasts, and the last three elements correspond to the 
			##   lower asymptote parameters for the three contrasts
			
			##   In the second formulation, the list elements can include
			##   character values corresponding to the various response
			##   models (e.g., drm, grm, gpcm, nrm, and mcm).  This
			##   indicates that the given parameter for all items associated
			##   with a given model should use the specified contrast
			if (is.character(unlist(con))) {
				for (i in 1:length(con)) {
					tmp <- NULL
					for (j in 1:length(mod)) {
						if (mod[j] %in% con[[i]]) tmp <- c(tmp, items[[j]])
					}
					con[[i]] <- tmp
				}
				
			##   In the third formulation, the list elements include item
			##   numbers, indicating that the given parameters for the 
			##   given item should use the specified contrast. All item
			##   numbers do not need to be included. In the case where
			##   no item number is specified for a given parameter and 
			##   contrast, the deviation contrast is used
			} else { 
				if (!is.numeric(unlist(con))) stop("Under formulation three, all the values must be numeric")
				
				##   Extract the list elements associated with the slope parameters
				ak <- con[1:3]
				
				##   Check to see if all items are included for this parameter
				##   type for one of the contrast types
				tmp <- unlist(ak)
				tmp1 <- c(1:n)%in%tmp
				
				##   If not, add missing item numbers to the list of deviation contrasts
				if (sum(tmp1)!=n) ak[[1]] <- c(ak[[1]], c(1:n)[tmp1==FALSE])
				
				
				
				##   Extract the list elements associated with the category parameters
				ck <- con[4:6]
				
				##   Check to see if all items are included for this parameter
				##   type for one of the contrast types
				tmp <- unlist(ck)
				tmp1 <- c(1:n)%in%tmp
				
				##   If not, add missing item numbers to the list of deviation contrasts
				if (sum(tmp1)!=n) ck[[1]] <- c(ck[[1]], c(1:n)[tmp1==FALSE])
				
				
				
				##   Extract the list elements associated with the lower asymptote parameters
				dk <- con[7:9]
				
				##   Check to see if all items are included for this parameter
				##   type for one of the contrast types
				tmp <- unlist(dk)
				tmp1 <- c(1:n)%in%tmp
				
				##   If not, add missing item numbers to the list of deviation contrasts
				if (sum(tmp1)!=n) dk[[1]] <- c(dk[[1]], c(1:n)[tmp1==FALSE])
				
				##   Recompile the list of contrasts
				con <- c(ak,ck,dk)
			}
		}
		names(con) <- c("dev.a","poly.a","tri.a","dev.c","poly.c","tri.c","dev.d","poly.d","tri.d")
		
		##   Reformat the contrast coefficients to traditional IRT parameters
		##   Loop through each item response model
		for (i in 1:length(mod)) {
		
			##   Loop through all items for the given model
			for (j in items[[i]]) {
			
				##   Dichotomous items. Only 3PL items need to
				##   be reformatted using the contrasts
				if (mod[[i]]=="drm") {
					if (drm.3PL==TRUE) {
						if (p.cat[items[[i]]][j]!=2) {
							p[[j]] <- p[[j]][1:3]
							p[[j]][2] <- p[[j]][2]/-p[[j]][1]
							p[[j]][1] <- p[[j]][1]/1.7
							if (j %in% con$tri.c) {
								p[[j]][3] <- exp(-p[[j]][3])/(1+exp(-p[[j]][3]))
							} else {
								p[[j]][3] <- exp(p[[j]][3])/(1+exp(p[[j]][3]))
							}
						}
					}
					p[[j]] <- c(p[[j]], rep(NA,col.max-length(p[[j]])))
					
				##   GRM items. These items are not affected
				##   by the contrasts, but reformulate the items if
				##   a location parameter is desired
				} else if (mod[[i]]=="grm") {
					if (loc.out==TRUE) {
						ck <- p[[j]][-1]
						ck <- c(mean(ck),ck-mean(ck))
						p[[j]] <- c(p[[j]][1],ck)
					}
					p[[j]] <- c(p[[j]], rep(NA,col.max-length(p[[j]])))
					
				##   GPCM, NRM, and MCM items
				} else if (mod[[i]] %in% c("gpcm","nrm","mcm")) {
					tmp <- p.cat[j]/3
					ak <- p[[j]][1:tmp]
					ck <- p[[j]][(tmp+1):(tmp*2)]
					dk <- p[[j]][(2*tmp+1):(3*tmp-1)]
					
					##   Reformat the slope parameters
					if (j %in% con$dev.a) {
						C <- cbind(-1/(tmp+1), matrix(-1/(tmp+1), tmp, tmp)+diag(rep(1,tmp)))
						if (mod[[i]]=="gpcm") {
							C <- C[,1]
						}
					} else if (j %in% con$poly.a) {
						C <- poly[[tmp]]
					} else if (j %in% con$tri.a) {
						C <- matrix(0, tmp, tmp+1)
						for (k in 1:(tmp)) {
							C[k,(k+1):(tmp+1)] <- -1
						}
					}
					ak <- ak%*%C
					
					##   Reformat the step/category parameters
					if (j %in% con$dev.c) {
						C <- cbind(-1/(tmp+1), matrix(-1/(tmp+1), tmp, tmp)+diag(rep(1,tmp)))
					} else if (j %in% con$poly.c) {
						C <- poly[[tmp]]
					} else if (j %in% con$tri.c) {
						C <- matrix(0, tmp, tmp+1)
						for (k in 1:(tmp)) {
							C[k,(k+1):(tmp+1)] <- -1
						}
					}
					ck <- ck%*%C
					if (mod[[i]]=="gpcm") ck <- ck[-1]
					
					##   Reformat the lower asymptote parameters
					if (j %in% con$dev.d) {
						C <- cbind(-1/tmp, matrix(-1/tmp, tmp-1, tmp-1)+diag(rep(1,tmp-1)))
					} else if (j %in% con$poly.d) {
						C <- poly[[tmp-1]]
					} else if (j %in% con$tri.d) {
						C <- matrix(0, tmp-1, tmp)
						for (k in 1:(tmp-1)) {
							C[k,(k+1):tmp] <- -1
						}
					}
					dk <- dk%*%C
					
					num <- exp(dk)
					dk <- num/sum(num)
					
					##   For the GPCM items, reformulated them to
					##   include a location parameter (if necessary)
					if (mod[[i]]=="gpcm") {
						ak <- ak[length(ak)]
						if (loc.out==TRUE) {
							ck <- c(mean(ck),ck-mean(ck))
						}
						p[[j]] <- c(ak, ck)
						
					##   For NRM and MCM items, include missing values (as needed)
					##   for the various parameter blocks so that the parameters will
					##   be formatted properly in the final matrix of parameters
					} else if (mod[[i]]=="nrm") {
						ak <- c(ak, rep(NA,a.max-length(ak)))
						ck <- c(ck, rep(NA,c.max-length(ck)))
						p[[j]] <- c(ak, ck)
					} else if (mod[[i]]=="mcm") {
						ak <- c(ak, rep(NA,a.max-length(ak)))
						ck <- c(ck, rep(NA,c.max-length(ck)))
						dk <- c(dk, rep(NA,d.max-length(dk)))
						p[[j]] <- c(ak, ck, dk)
					}
					
					##   Add  missing values to the right of the parameters for
					##   each item (as needed) so that all of the list elements
					##   can be stacked on top of one another to create the 
					##   final matrix of item parameters
					p[[j]] <- c(p[[j]], rep(NA,col.max-length(p[[j]])))
				} 
			}
		}
		
		##   Initialize an object to store the compiled set of item parameters
		pars <- NULL
		
		##   Loop through all the items
		for (i in 1:length(p)){
			pars <- rbind(pars,p[[i]])
		}
		
		##   Create an {irt.pars} object
		if (as.irt.pars==TRUE) {
			pars <- as.irt.pars(pars, cat=cat, poly.mod=poly.mod, location=loc.out)
		}
	
	##   Import ability estimates
	} else {
		pars <- read.fwf(file, c(10,10,4))
		names(pars) <- c("theta","theta.se","freq")
	}
	return(pars)
}




##   Import item parameters or ability estimates from TESTFACT 4
read.testfact <- function(file, ability=FALSE, guessing=FALSE, bifactor=FALSE, as.irt.pars=TRUE) {

	##   Import item parameters
	if (ability==FALSE) {
	
		##   Read in the first row of item parameters
		tmp1 <- scan(file, skip=1, what="character", quiet=TRUE, nlines=1)
		
		##   Read in the second row of item parameters
		tmp2 <- scan(file, skip=2, what="character", quiet=TRUE, nlines=1)
		
		##   In some instances the parameters for a given item can
		##   wrap onto a second line. Check to see if the first element 
		##   of the second row is "2" (indicating the second item). If not,
		##   the parameters are wrapped and the number of values
		##   associated with each item needs to be adjusted; that is,
		##   set equal to the length of the elements in the first two rows
		##   instead of just the length of the elements in the first row
		if (tmp2[1]!="2") tmp1 <- c(tmp1,tmp2)
		
		##   Read in all of the data from the .PAR file
		pars <- scan(file, skip=1, what="character", quiet=TRUE)
		
		##   Number of items
		items <- length(pars)/length(tmp1)
		
		##   Initialize a matrix to store the item parameters
		pars <- matrix(pars, items, length(tmp1), byrow=TRUE)
		
		##   Eliminate the item number and name
		pars <- pars[,-c(1,2)]
		
		if (guessing==TRUE) {
			##   Number of dimensions
			dimensions <- length(tmp1)-4
			
			##   Respecify the values to be numeric
			pars <- matrix(as.numeric(pars), items, dimensions+2)
			
			##   Re-order the columns
			pars <- pars[,c(3:(dimensions+2),1,2)]
			
			colnames(pars) <- c(paste("a",1:dimensions,sep=""),"d","c")
		} else {
			##   Number of dimensions
			dimensions <- length(tmp1)-3
			
			##   Respecify the values to be numeric
			pars <- matrix(as.numeric(pars), items, dimensions+1)
			
			##   Re-order the columns
			pars <- pars[,c(2:(dimensions+1),1)]
			
			colnames(pars) <- c(paste("a",1:dimensions,sep=""),"d")
		}
		
		##   Create an {irt.pars} object
		if (as.irt.pars==TRUE) {
			n <- nrow(pars)
			pm <- as.poly.mod(n)
			pars <- as.irt.pars(pars, cat=rep(2,n), poly.mod=pm, dimensions=dimensions)
		}
		
	##   Import ability estimates
	} else {
	
		##   Read in the data
		tmp <- scan(file,what="character", quiet=TRUE)
		
		if (bifactor==FALSE) {
			
			##   Initialize an object to flag the posterior SDs
			##   (i.e., the values with asterisks)
			flag <- NULL
			
			##   Loop through all of the character strings
			for (i in 1:length(tmp)) {
			
				##   Check the length of the character string
				if (nchar(tmp[i])==6) {
				
					##   Check for the asterisk
					if (substr(tmp[i],6,6)=="*") {
						flag <- c(flag,i)
					}
				}
			}
			##   Remove all posterior SDs
			tmp <- tmp[-flag]
		}
		
		##   The first five elements are not factor scores
		##   Start with the fifth element and loop through
		##   all of the values until "2" is found, indicating
		##   the start of the data for examinee 2. Use this
		##   information to determine the number of dimensions
		for (i in 6:20) {
			if (tmp[i]=="2") dimensions <- i-6
		}
		
		##   Format the factor scores as a matrix
		pars <- matrix(tmp, ncol=dimensions+5, byrow=T)
		
		##   Eliminate the examinee information (ID numbers etc)
		pars <- pars[,-c(1:5)]
		
		##   Respecify the values to be numeric
		pars <- matrix(as.numeric(pars), nrow(pars), dimensions)
		
		colnames(pars) <- paste("theta",1:dimensions,sep="")
	}
	
	return(pars)
}




##   Reformat item parameters from the R package eRm
read.erm <- function(x, loc.out=FALSE, as.irt.pars=TRUE) {

	groups <- x$ngroups
	time <- x$mpoints
	items <- ncol(x$X)/time
	cat <- apply(x$X,2,max,na.rm=TRUE)[1:items]
	beta <- x$betapar
	
	##   Identify the item, category, time, and group ordering of the beta parameters
	it <- rep(rep(1:items,cat),groups*time)
	t <- rep(1:2,each=sum(cat)*time)
	g <- rep(rep(1:2,each=sum(cat)),time)
	
	##   Combine the parameters into a matrix
	pars <- vector("list",groups*time)
	comb <- expand.grid(list(1:groups,1:time))
	comb$grp <- 1:nrow(comb)
	for (i in 1:items) {
		for (j in 1:groups) {
			for (k in 1:time) {
				tmp <- c(beta[it==i & g==j & t==k],rep(NA,max(cat)-cat[i]))
				pars[[comb$grp[comb[,1]==j & comb[,2]==k]]] <- rbind(pars[[comb$grp[comb[,1]==j & comb[,2]==k]]],tmp)
				names(pars)[comb$grp[comb[,1]==j & comb[,2]==k]] <- paste("group",j,".time",k,sep="")
			}
		}
	}
	
	##   Update {cat} to correspond to the number of response categories
	cat <- cat+1
	names(cat) <- NULL
	
	##   Create the poly.mod object
	if (min(cat)==2 & max(cat==2)) {
		mod <- "drm"
		tmp.it <- 1:items
	} else if (min(cat)==2 & max(cat>2)) {
		mod <- c("drm","gpcm")
		tmp <- 1:items
		tmp.it <- list(tmp[cat==2],tmp[cat>2])
	} else if (min(cat)>2 & max(cat>2)) {
		mod <- "gpcm"
		tmp.it <- 1:items
	}
	pm <- as.poly.mod(items,mod,tmp.it)
	
	##   Reformat the item parameters as a {sep.pars} object
	for (i in 1:length(pars)) {
		rownames(pars[[i]]) <- NULL
		colnames(pars[[i]]) <- NULL
		pars[[i]] <- sep.pars(pars[[i]], cat=cat, poly.mod=pm, loc.out=loc.out)
	}
	
	##   Create an {irt.pars} object
	if (length(pars)==1) {
		pars <- as.irt.pars(pars[[1]])
	} else {
		common <- vector("list",(groups*time)-1)
		for (i in 1:length(common)) {
			common[[i]] <- matrix(1:items,items,2)
		}
		pars <- as.irt.pars(pars,common,grp.names=names(pars))
	}
	
	if (as.irt.pars==FALSE) {
		pars <- pars@pars
	}
	
	return(pars)
}




##   Reformat item parameters from the R package ltm
read.ltm <- function(x, loc.out=FALSE, as.irt.pars=TRUE) {

	##   Identify the class of the parameter object
	##   (different models are saved as different classes)
	cls <- class(x)
	
	dimensions <- 1
	
	##   Rasch model
	if (cls=="rasch") {
		pars <- cbind(coef(x)[,2:1])
		colnames(pars) <- c("a","b")
		cat <- rep(2,nrow(pars))
		cls <- "drm"
	
	##   3PL
	} else if (cls=="tpm") {
		pars <- coef(x)[,3:1]
		colnames(pars) <- c("a","b","c")
		cat <- rep(2,nrow(pars))
		cls <- "drm"
		
	##   Graded response model
	} else if (cls=="grm") {
		pars <- coef(x)
		
		##   All items have the same number of response categories
		if (is.matrix(pars)) {
			##   Re-order the columns
			pars<- pars[,c(ncol(pars),1:(ncol(pars)-1))]
			
			##   Identify the number of response categories
			cat <- rep(ncol(pars)-1,nrow(pars))
			cat[cat==1] <- 2
		
		##   There are different numbers of response
		##   categories across items
		} else {
			##   Identify the number of response categories for each item
			cat <- unlist(lapply(pars,length))
			names(cat) <- NULL
			
			##   Initialize a matrix to store all of the item parameters
			p <- matrix(NA,length(pars),max(cat))
			
			##   Loop through all of the items and populate {p}
			for (i in 1:length(pars)) {
				p[i,1:cat[i]] <- pars[[i]][length(pars[[i]]):1]
			}
			pars <- p
		}
		
		##   If the parameters are formatted using the slope/intercept
		##   parameterization, reformat them to use a slope/difficulty parameterization
		if (x$IRT.param==FALSE) {
			pars[,2:ncol(pars)] <- pars[,2:ncol(pars)]/matrix(pars[,1],nrow(pars),ncol(pars)-1)
		}
		colnames(pars) <- c("a",paste("b",1:(ncol(pars)-1),sep=""))
		
	##   Generalized partial credit model
	} else if (cls=="gpcm") {
		pars <- coef(x)
		
		##   All items have the same number of response categories
		if (is.matrix(pars)) {
			##   Re-order the columns
			pars<- pars[,c(ncol(pars),1:(ncol(pars)-1))]
			
			##   Identify the number of response categories
			cat <- rep(ncol(pars)-1,nrow(pars))
			cat[cat==1] <- 2
			
		##   There are different numbers of response
		##   categories across items
		} else {
			##   Identify the number of response categories for each item
			cat <- unlist(lapply(pars,length))
			names(cat) <- NULL
			
			##   Initialize a matrix to store all of the item parameters
			p <- matrix(NA,length(pars),max(cat))
			
			##   Loop through all of the items and populate {p}
			for (i in 1:length(pars)) {
				p[i,1:cat[i]] <- pars[[i]][length(pars[[i]]):1]
			}
			pars <- p
		}
		
		##   If the parameters are formatted using the slope/intercept
		##   parameterization, reformat them to use a slope/difficulty parameterization
		if (x$IRT.param==FALSE) {
			pars[,2:ncol(pars)] <- -pars[,2:ncol(pars)]/matrix(pars[,1],nrow(pars),ncol(pars)-1)
		}
		colnames(pars) <- c("a",paste("b",1:(ncol(pars)-1),sep=""))
		
	##   2PL and M2PL
	} else if (cls=="ltm") {
	
		##   slope/difficulty parameterization for the 2PL
		if (x$IRT.param==TRUE) {
			pars <- cbind(coef(x)[,2:1],0)
			colnames(pars) <- c("a","b","c")
			
		##   slope/intercept parameterization
		} else {
			##   2PL
			if (x$ltst$factors==1) {
				pars <- cbind(coef(x)[,2],-coef(x)[,1]/coef(x)[,2],0)
				colnames(pars)[1:2] <- c("a","b","c")
				
			##   M2PL (only two dimensions)
			} else {
				pars <- cbind(coef(x)[,3:1])
				colnames(pars)[1:3] <- c("a1","a2","d")
				dimensions <- 2
			}
		}
		cat <- rep(2,nrow(pars))
		cls <- "drm"
	}
	
	pm <- as.poly.mod(nrow(pars), cls)
	pars <- sep.pars(pars, cat=cat, poly.mod=pm, dimensions=dimensions, loc.out=loc.out)
	
	##   Create an {irt.pars} object
	pars <- as.irt.pars(pars)
	
	if (as.irt.pars==FALSE) {
		pars <- pars@pars
	}
	
	return(pars)
}



##   Import item parameters or ability estimates from ICL
read.icl <- function(file, poly.mod, ability=FALSE,  loc.out=FALSE, as.irt.pars=TRUE) {

	##   Import item parameters
	if (ability==FALSE) {
	
		##   Read in the data
		pars <- read.table(file,sep="\t",fill=TRUE)
		
		##   Add columns of NAs to the right to adjust for 
		##   item parameters that got wrapped to the
		##   next row during the import
		pars <-cbind(pars, matrix(NA,nrow(pars),2*ncol(pars)))
		
		##   Create an object to identify the current item number
		id <- 1
		
		##   Loop through the rows and place any item parameters
		##   that got wrapped into the same row as the rest
		##   of the parameters for the given item
		repeat{
			##   Check to see if the first value in the given row
			##   is equal to the expected item number (id)
			if (pars[id,1]==id) {
				if (id==nrow(pars)) break
				id <- id+1
				
			##   If not, item parameters have been wrapped
			} else {
				tmp <- pars[id,!is.na(pars[id,])]
				tmp1 <- length(pars[id-1,!is.na(pars[id-1,])])
				pars[id-1,(tmp1+1):(tmp1+length(tmp))] <- tmp
				pars <- pars[-id,]
				if ((id-1)==nrow(pars)) break
			}
		}
		
		##   Eliminate columns where there are no parameters
		pars <- pars[,apply(is.na(pars),2,sum)<nrow(pars)]
		
		##   Remove the column of item numbers
		pars <- pars[,-1]
		
		##   Determine the number of response categories
		cat <- apply(!is.na(pars),1,sum)
		cat[poly.mod@items$drm] <- 2
		
		pars <- sep.pars(pars, cat=cat, poly.mod=poly.mod, loc.out=loc.out)
		
		##   Create an {irt.pars} object
		pars <- as.irt.pars(pars)
		
		if (as.irt.pars==FALSE) {
			pars <- pars@pars
		}
		
	##   Import ability estimates
	} else {
		pars <- read.table(file, sep="\t")
	}
	
	return(pars)
}



##   Import item parameters or ability estimates from BMIRT
read.bmirt <- function(file, ability=FALSE, sign.adjust=TRUE, loc.out=FALSE, pars.only=TRUE, as.irt.pars=TRUE) {

	##   Import item parameters
	if (ability==FALSE) {
	
		##   Read in the data as a vector
		tmp <- scan(file, skip=1, quiet=TRUE)
		nc <- 40
		
		##   Initialize an object to store the item parameters
		pars <- NULL
		
		##   Each row in the .PAR file ends with 1.00  1.00 0.00
		##   Loop through all of the values and find this
		##   sequence to identify start and end ranges for the 
		##   parameters associated with each item
		
		##   Value indexing the location to start for 
		##   extracting the parameters
		start <- 1
		
		##   Checking for this sequence starting with the third value. This
		##   avoids problems where the loading on the first dimension is
		##   zero for the first item
		for (i in 3:(length(tmp)-2)) {
		
			##   If this sequence is found
			if (tmp[i]==1 & tmp[i+1]==1 & tmp[i+2]==0) {
				len <- i-start
				
				##   Extract the appropriate item parameters
				pars <- rbind(pars,c(tmp[start:(i-1)],rep(NA,nc-len)))
				
				##   Increment the start index for the next set of parameters
				start <- i+3
			}
		}
		
		##   Eliminate columns where there are no parameters
		tmp <- apply(is.na(pars),2,sum)
		pars <- pars[,tmp!=nrow(pars)]
		
		##   Extract the number of response categories
		cat <- pars[,2]
		
		##   Eliminate the columns of item numbers and response categories
		pars <- pars[,-c(1,2)]
		
		##   Determine the number of dimensions
		tmp <- pars[1,!is.na(pars[1,])]
		
		##   M2PL
		if (cat[1]==2) {
			dimensions <- length(tmp)-1
		##   M3PL
		} else if (cat[1]==1) {
			dimensions <- length(tmp)-2
		##   MGPCM
		} else {
			dimensions <- length(tmp)-cat[1]+1
		}
		
		##   Recode the number of response categories for  M3PL items
		cat[cat==1] <- 2
		
		##   Number of items
		n <- nrow(pars)
		
		##   Create the poly.mod object
		if (length(cat[cat==2])==n) {
			##   All dichotomous items
			pm <- as.poly.mod(n)
		} else {
			##   Mixed-format items
			items <- 1:n
			pm <- as.poly.mod(n, c("drm","gpcm"), list(items[cat==2], items[cat>2]))
		}
		pars <- sep.pars(pars, cat=cat, poly.mod=pm, dimensions=dimensions, loc.out=loc.out)
		
		##   Reformat the difficulty/step parameters to coincide with the 
		##   traditional formulation of unidimensional and multidimensional models
		if (dimensions==1) {
			pars@b[cat==2] <- pars@b[cat==2]/pars@a[cat==2]
		} else if (dimensions>1) {
			if (sign.adjust==TRUE) pars@b <- pars@b*-1
		}
		
		##   Create an {irt.pars} object
		pars <- as.irt.pars(pars)
		
		if (as.irt.pars==FALSE) {
			pars <- pars@pars
		}
		
	##   Import ability estimates
	} else {
		pars <- read.table(file, sep=" ", skip=2)
		
		##   Eliminate columns where there are no parameters
		tmp <- apply(is.na(pars),2,sum)
		pars <- pars[,tmp!=nrow(pars)]
		
		##   Eliminate the column of examinee numbers
		pars <- pars[,-1]
		
		##   Number of dimensions
		dimensions <- ncol(pars)/2
		
		nms <- c(paste("theta",rep(1:dimensions),sep=""), paste("theta",rep(1:dimensions),".se",sep=""))
		colnames(pars) <- nms
		
		##   Remove the columns with the standard errors
		if (pars.only==TRUE) pars <- pars[,1:dimensions]
	}
	
	return(pars)
}
