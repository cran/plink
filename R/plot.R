plot.irt.pars <- function(x, y, ..., type, separate, combine, items, item.names, item.nums, panels, drift, groups, grp.names, sep.mod, drift.sd, save.hist) {
	
	if (!missing(y)) {
		if (is.irt.pars(y)) stop("It looks like you are trying to compare the parameters from two {irt.pars} objects.\n In the current specification there is no way of knowing which items are common.\n Create a single {irt.pars} object using the combine.pars() function then try again.")
	}
	
	##   Delete plot history
	if (!missing(save.hist)) {
		if (save.hist==FALSE) {
			if (exists(".SavedPlots",where=1)) rm(".SavedPlots",pos=1)
		}
	}
	
	##   Check to see if another graphics device is already open
	if (names(dev.cur())=="null device") dev.flag <- FALSE else dev.flag <- TRUE
	
	##   Create an object to indicate whether the saving
	##   of plot histories has already been specified
	record.flag <- FALSE
	
	##   Create an object to indicate whether the margins for 
	##   a plot have already been initialized (if this is not 
	##   checked, there can be problems in creating multi-group plots)
	pan.flag <- FALSE
	
	##   Set the default type of MD plot
	if (missing(type)) type <- "wireframe"
	
	##   Extract/create values to be passed along to the {mixed} function
	##   for computing response probabilities (if applicable)
	dots <- list(...)
	if (length(dots$catprob)) catprob <- dots$catprob else catprob <- FALSE
	if (length(dots$incorrect)) incorrect <- dots$incorrect else incorrect <- FALSE
	if (length(dots$D)) D <- dots$D else D <- 1
	if (length(dots$D.drm)) D.drm <- dots$D.drm else D.drm <- D
	if (length(dots$D.gpcm)) D.gpcm <- dots$D.gpcm else D.gpcm <- D
	if (length(dots$D.grm)) D.grm <- dots$D.grm else D.grm <- D
	
	##   Check to see of the argument {drift} is missing
	##   if it is, create a plot other than a comparison plot
	if (missing(drift)) {
		##   Set the default number of panels to 20
		if (missing(panels)) panels <-20
		
		##   For the plots that utilize them, set item.nums to TRUE
		if (missing(item.nums)) item.nums <- TRUE
		
		##   Number of groups
		ng <- x@groups
		
		##   Separate out the item parameters for each group
		pars <- sep.pars(x)
		
		if (ng==1) {
			x@pars <- list(x@pars)
			pars <- list(pars)
		}
		
		##   Initialize an object to store the trellis objects
		pl.out <- list()
		
		##   When there are multiple groups, if values are supplied for
		##   {items}, they should be in the form of a list. If not, all
		##   items in each group will be used. Create a flag for the 
		##   group iterations below to skip over the warning and just
		##   use all items (if applicable)
		items.flag1 <- FALSE
		
		##   When there are multiple groups, if values are supplied for
		##   {item.names}, they should be in the form of a list. If not, the
		##   default item names will be used. Create a flag for the 
		##   group iterations below to skip over the warning and just
		##   use the item names (if applicable)
		items.flag2 <- FALSE
		
		
		##   Loop through all of the groups
		for (grp in 1:ng) {
		
			##   Number of dimensions for the given group
			dimensions <- x@dimensions[grp]
			
			##   Use all items
			if (missing(items)) {
				gr.items <- 1:nrow(x@pars[[grp]])
				
			##   Use some subset of items
			} else {
				if (ng>1) {
					if (items.flag1==FALSE) {
						if(!is.list(items)) {
							warning("{x} has more than one group, but {items} is not a list. All items used for each group")
							gr.items <- 1:nrow(x@pars[[grp]])
							items.flag1 <- TRUE
						} else {
							gr.items <- items[[grp]]
						}
					} else {
						##   use all items
						gr.items <- 1:nrow(x@pars[[grp]])
					}
				} else {
					gr.items <- items
				}
			}
			
			##   Create a default set of item names
			if (missing(item.names)) {
				grn.items <- paste("Item",1:nrow(x@pars[[grp]]))
				grn.items <- grn.items[gr.items]
				
			##   Use a set of user-specified item names
			} else {
				if (ng>1) {
					if (items.flag2==FALSE) {
						if(!is.list(item.names)) {
							warning("{x} has more than one group, but {item.names} is not a list. Default item names used")
							grn.items <- paste("Item",1:nrow(x@pars[[grp]]))
							grn.items <- grn.items[gr.items]
							items.flag2 <- TRUE
						} else {
							grn.items <- item.names[[grp]][gr.items]
						}
					} else {
						##   Create a default set of item names
						grn.items <- paste("Item",1:nrow(x@pars[[grp]]))
						grn.items <- grn.items[gr.items]
					}
				} else {
					grn.items <- item.names
				}
			}
			
			##   Check to see if a vectorplot should be created
			if (type %in% c("vectorplot1","vectorplot2","vectorplot3")) {
			
				if (dimensions<7) {
					if (dimensions==1) stop("Vector plots are not supported for single dimensions")
					
					##   Number of items for the given group
					ni <-length(gr.items)
					
					##   Number of permutations of two-dimensional axes to create
					##   One for each pair of dimensions
					perm <- 0
					for (i in 1:ng) {
						for (j in 1:(x@dimensions[i]-1)) {
							for (k in (j+1):x@dimensions[i]) {
								perm <- perm+1
							}
						}
					}
					
					##   Check to see if the number of permutations is less than
					##   the number of panels per page
					if (perm<panels) {
						
						##   If there are fewer permutations, set the number
						##   of panels equal to the number of permutations
						panels <- perm
						
						##   Check to see if another graphics device is already open
						if (dev.flag==FALSE) {
							if (ng>1) {
							
								##   Check to see if the saving of plot histories
								##   has already been initialized (from a previous group)
								if (record.flag==FALSE) {
								
									##   If not, record the plot history (this only works for Windows)
									if (Sys.info()["sysname"] == "Windows") {
										cat("Use PgUp and PgDn to view different plot pages\n")
										windows(record=TRUE)
										record.flag <- TRUE
									}
								}
							}
						}
					
					##   If there are more permutations than panels or the same number
					} else {
						
						##   Check to see if another graphics device is already open
						if (dev.flag==FALSE) {
						
							##   Check to see if the saving of plot histories
							##   has already been initialized (from a previous group)
							if (record.flag==FALSE) {
							
								##   If not, record the plot history (this only works for Windows)
								if (Sys.info()["sysname"] == "Windows") {
									cat("Use PgUp and PgDn to view different plot pages\n")
									windows(record=TRUE)
									record.flag <- TRUE
								}
							}
						}
					}
					
					##   Reset the number of panels if greater than 36
					if (panels>36) {
						warning("The number of panels cannot exceed 36")
						panels <- 20
					}
					
					##   Identify the dimensions of panels (i.e., how to organize the panels on the page)
					pan <- matrix(c(rep(1:6,c(2,2,5,7,9,11)),rep(1:6,c(1,5,6,8,10,6))),36,2)
					
					##   Check to see in the margins for the panel have already been initialized
					##   (Re-initializing them creates problems for multiple groups)
					if (pan.flag==FALSE) {
						
						##   If not, set the numbers of rows and columns (of plots)
						##   to be created on a single page, and set the outer margins
						par(mfrow=c(pan[panels,1],pan[panels,2]), mai=c(0.25,0.25,0.25,0.25))
						pan.flag <- TRUE
					}
					
					##   Extract the slope parameters for the subset of item for the given group
					a.all <- as.matrix(pars[[grp]]@a)[gr.items,]
					
					##   Extract the slope parameters for the subset of item for the given group
					b <- as.matrix(as.matrix(pars[[grp]]@b)[gr.items,])
					
					##   Lower asymptote parameters are not needed for the vectorplots
					
					##   Transpose  the object {a.all} if there is only one item
					##   this will keep the parameters formatted as a 1 x m matrix
					if (is.vector(a.all) & ni==1) a.all <- t(a.all)
					
					
					##   Loop through all pairs of dimensions
					for (i in 1:(dimensions-1)) {
						for (j in (i+1):dimensions) {
						
							##   Extract the slopes for each of the given dimensions
							a <- a.all[,c(i,j)]
							if (is.vector(a) & ni==1) a <- t(a)
							
							##   Compute multidimensional discrimination
							mdisc <- sqrt(apply(a^2,1,sum,na.rm=T))
							
							##   Compute directional cosine
							dcos <- a/matrix(mdisc,length(mdisc),ncol(a))
							if (is.vector(dcos) & ni==1) dcos <- t(dcos)
							
							##   Compute multidimensional difficulty
							den.con <- apply(!is.na(b),1,sum)
							mdiff <- -apply(b,1,sum,na.rm=T)/(den.con*mdisc)
							
							##   Remove items that do not load on either of the plotted dimensions
							a <- a[mdisc!=0,]
							if (is.vector(a) & ni==1) a <- t(a)
							dcos <- dcos[mdisc!=0,]
							if (is.vector(dcos) & ni==1) dcos <- t(dcos)
							mdiff <- mdiff[mdisc!=0]
							grn.items <- grn.items[mdisc!=0]
							mdisc <- mdisc[mdisc!=0]
							
							##   Create a set of X and Y coordinates for each vector
							x1 <- dcos[,1]*mdiff
							y1 <- dcos[,2]*mdiff
							x2 <- mdisc*dcos[,1]+x1
							y2 <- mdisc*dcos[,2]+y1
							
							##   Create a set of X and Y coordinates for item numbers
							md <- mdiff+0.05
							md[mdiff<0] <- mdiff[mdiff<0]-0.05
							if (type!="vectorplot3")  md <- mdiff-0.05
							x3 <- dcos[,1]*md
							y3 <- dcos[,2]*md
							
							##   Compute the reference composite for vectorplot2
							ref <- abs(eigen(t(a)%*%a)$vectors[,1])
							
							##   Vectorplot3 - This plots arrows that originate at the 
							##   origin and extend to the point of MDIF in the direction
							##   of the directional cosine
							if (type=="vectorplot3") {
							
								##   Identify the min and max values to use for the axes scales
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
							
								##   Identify the min and max values to use for the axes scales
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
							
							##   Compile the scale ranges for both axes
							xmm <- c(xmin,xmax)
							ymm <- c(ymin,ymax)
							
							##   Check to see if the user supplied values for the
							##   axes scales (change the computed values as necessary)
							if (length(dots$xlim)) xmm <- dots$xlim
							if (length(dots$ylim)) ymm <- dots$ylim
							
							##   Plot crosshairs
							plot(50,50, axes=FALSE,xlim=xmm,ylim=ymm,
							mar=rep(0.1,4),xlab="",ylab="")
							axis(1,pos=0,at=seq(xmm[1],xmm[2],0.5))
							axis(2,pos=0,at=seq(ymm[1],ymm[2],0.5))
							
							##   Plot the arrows
							if (type=="vectorplot1") {
								suppressWarnings(arrows(x1,y1,x2,y2,length=.1))
							} else if (type=="vectorplot2") {
								suppressWarnings(arrows(x1,y1,x2,y2,length=.1))
								arrows(max(par("usr")[c(1,3)])*ref[1],max(par("usr")[c(1,3)])*ref[2],min(par("usr")[c(2,4)])*ref[1],
									min(par("usr")[c(2,4)])*ref[2],lwd=3,length=.2,col="darkred")
							} else if (type=="vectorplot3") {
								tmp <- rep(0,length(mdisc))
								suppressWarnings(arrows(tmp,tmp,x1,y1,length=.1))
							}
							
							##   Plot item numbers (if applicable)
							if (item.nums==TRUE) text(x3,y3,gr.items, col=2)
							t1 <- (abs(xmm[1])/abs(xmm[1]))*0.02
							t2 <- (abs(ymm[1])/abs(ymm[1]))*0.03
							
							##   Print theta labels for the two axes for the given pair of groups
							if (ng>1) {
								text(xmm[2],t2*ymm[2], eval(parse(text=paste("expression(paste(theta[G",grp,".",i,"]))",sep=""))),cex=1.2,pos=3)
								text(t1*xmm[2],ymm[2], eval(parse(text=paste("expression(paste(theta[G",grp,".",j,"]))",sep=""))),cex=1.2,pos=4)
							} else {
								text(xmm[2],t2*ymm[2], eval(parse(text=paste("expression(paste(theta[",i,"]))",sep=""))),cex=1.2,pos=3)
								text(t1*xmm[2],ymm[2], eval(parse(text=paste("expression(paste(theta[",j,"]))",sep=""))),cex=1.2,pos=4)
							}
						}
					}
				} else {
					stop("Vector plots are not supported for more than six dimensions")
				}
				
			##   Create information plots
			} else if (type %in% c("information1","information2")) {
				##   This will be implemented in a future release
				##   The two plots will be information surfaces and clamshells
			
			##   Create item response curves/surfaces 
			##   These are based on the plot.irt.prob function
			} else {
				if (length(dots$theta)) {
					theta <- dots$theta 
				} else {
					if (dimensions==1) {
						theta <- seq(-4,4,.05) 
					} else if (dimensions %in% 2:3) {
						theta <- seq(-4,4,.5)
					} else {
						theta <- -4:4
					}
				}
				
				##   Compute response probabilities
				tmp <- mixed(pars[[grp]], theta=theta, catprob=catprob, location=pars[[grp]]@location, incorrect=incorrect, D.drm=D.drm, D.gpcm=D.gpcm, D.grm=D.grm)
				
				##   If there is more than one group, compile the trellis objects in
				##   a list then print them later
				if (ng==1) {
					plot(tmp,type=type,separate=separate,combine=combine, items=items, item.names=item.names, item.nums=item.nums,panels=panels,...)
				} else {
					if (names(x@pars)[1]=="group1") {
						nms <- paste("G",grp,": ",grn.items,sep="")
					} else {
						nms <- paste(names(x@pars)[grp],": ",grn.items,sep="")
					}
					pl.out[[grp]] <- plot(tmp,type=type,separate=separate,combine=combine, items=items, item.names=nms, item.nums=item.nums,panels=panels,mg=1,...)
				}
			}
		}
		
		##   Print compiled trellis objects
		if (length(pl.out)) {
			for (i in 1:ng) {
				print(pl.out[[i]])
			}
		}
		
	##   Create plots comparing common item parameters and item/test characteristic curves/surfaces
	} else { 
	
		##   These plots cannot be created for only one group
		if (x@groups==1) {
			stop("To examine item parameter drift, {x} must include parameters for two or more groups.")
		} else {
			if (missing(sep.mod)) sep.mod <- FALSE
			if (missing(groups)) groups <- 1:(x@groups-1)
			
			##   Make sure there are fewer group pairs than the total number of groups
			groups <- groups[groups<x@groups]
			if (length(groups)==0) warning("Invalid group number(s). Plots will be created for all groups.")
			
			##   Separate out the item parameters for each group
			pars <- sep.pars(x)
			
			##   Rename the default group names (in the irt.pars object)
			if (names(x@pars)[1]=="group1") names(x@pars) <- paste("Group",1:length(x@pars))
			
			##   Check to see what values were supplied for the {drift} argument
			##   Reformat the "pars" value
			if ("pars" %in% drift) {
				tmp <- drift[drift!="a" & drift!="b" & drift!="c" & drift!="pars"]
				drift <- c("a","b","c",tmp)
			}
			
			##   Reformat the "all" value
			if ("all" %in% drift) drift <- c("a","b","c","TCC","ICC")
			
			##   Compute response probabilities for the TCCs and ICCs
			if ("TCC" %in% drift|"ICC" %in% drift) {
				prob <- mixed(x, theta=theta, catprob=catprob, incorrect=incorrect, D.drm=D.drm, D.gpcm=D.gpcm, D.grm=D.grm)
			}
			
			##   Check to see if a value is specified for {drift.sd}
			if (missing(drift.sd)) drift.sd <- 3
			
			##   Make sure the common item object is a list
			if (is.matrix(x@common)) x@common <- list(x@common)
			
			if (missing(grp.names)) nms <- names(x@pars) else nms <- grp.names
			
			##   Check to see is a graphics device has already been initialized
			if (names(dev.cur())=="null device") {
			
				##   For Windows, save the plot history so that different plot
				##   windows can be seen using PgUp and PgDn
				if (Sys.info()["sysname"] == "Windows") {
					if (length(groups)>1|length(drift)>1) {
						cat("Use PgUp and PgDn to view different plot pages\n")
						windows(record=TRUE)
					}
				}
			}
			
			##   The bult-in SD function R uses n-1 in the denominator
			##   This function uses n in the denominator
			.sd <- function(x) {
				z <- x[!is.na(x)]
				out <- sqrt(sum((z-mean(z))^2)/length(z))
				return(out)
			}
			
			##   Create an object to store the trellis output for each group
			pl.out <- list()
			
			##   Panel function uses to specify the markers and lines 
			##   used in the parameter comparison plots
			pfunc <- function(x,y,...,sd,item.names, drift.sd) {
				if (sep.mod==TRUE) {
					panel.xyplot(x,y,cex=1.2,pch=c(1,2,0,6,8),...)
				} else {
					panel.xyplot(x,y,cex=1.2,pch=1,...)
				}
				
				##   Plot SD line
				slope <- .sd(y)/.sd(x)
				int <- mean(y)-slope*mean(x)
				panel.abline(int,slope)
				
				##   Create a confidence interval that is 3 SDs
				##   around the center line in perpendicular distance
				pd <- (drift.sd*sd)/sin(atan(1/slope))
				panel.abline(int-pd,slope,lty=2)
				panel.abline(int+pd,slope,lty=2)
				
				##   Print item numbers
				ydif <- (range(y)[2]-range(y)[1])*0.025
				if (length(item.names)) ltext(x,y-ydif,item.names,col="black",cex=.7)
			}
			
			##   Set the marker options to be used in the legend (if applicable)
			if (sep.mod==TRUE) {
				sps <- trellis.par.get("superpose.symbol")
				sps$pch <- c(1,2,0,6,8)
				sps$cex <- 1.2
				trellis.par.set("superpose.symbol", sps)
			} else {
				sps <- trellis.par.get("superpose.symbol")
				sps$pch <- 1
				sps$cex <- 1
				trellis.par.set("superpose.symbol", sps)
			}
			
			
			##   Loop through all adjacent groups
			for (i in 1:length(groups)) {
				
				grp <- groups[i]
				
				##   Identify the item numbers for the lower group
				items1 <- x@common[[grp]][,1]
				
				##   Identify the item numbers for the higher group
				items2 <- x@common[[grp]][,2]
				
				if (missing(item.nums)) {
					item.names <- NULL
				} else {
					if (item.nums==TRUE) {
						if (missing(item.names)) {
							item.names <- paste(items1, ",", items2,sep="")
						}
					} else {
						item.names <- NULL
					}
				}
				
				##   Initialize a list to store the trellis objects for each
				##   element included in {drift}
				pl.out[[i]] <- vector("list",length(drift))
				
				##   Initialize an object to increment the plot index
				##   That is, a variable that indicates which list element
				##   of {pl.out[[i]]} should hold the given trellis object
				##   "a", "b", "c", etc.
				pl <- 1
				
				##   Compare slope parameters
				if ("a" %in% drift) {
				
					##   Extract slopes for the lower group
					pa1 <- pars[[grp]]@a[items1,]
					
					##   Extract slopes for the higher group
					pa2 <- pars[[grp+1]]@a[items2,]
					
					##   Transpose these values if there is only one item
					if (is.vector(pa1)) {
						if (length(items1)==1) {
							pa1 <- t(pa1)
							pa2 <- t(pa2)
						} else {
							pa1 <- as.matrix(pa1)
							pa2 <- as.matrix(pa2)
						}
					}
					
					##   Create a title
					if ("drm" %in% pars[[grp]]@model) {
						if (length(pars[[grp]]@model)>1) {
							if (max(x@dimensions)>1) {
								a.title <- "Slope Parameters"
							} else {
								a.title <- "Discrimination/Slope Parameters"
							}
						} else {
							if (max(x@dimensions)>1) {
								a.title <- "Slope Parameters"
							} else {
								a.title <- "Discrimination Parameters"
							}
						}
					} else {
						a.title <- "Slope Parameters"
					}
					
					##   Initialize an indicator variable to identify the item
					##   response model associated with each common item
					m.type <- pa1
					
					##   Create an object that includes an incremented item
					##   number alongside the item numbers for the selected subset
					items1a=cbind(1:length(items1),items1)
					
					##   Initialize an object to store the response models for the subset of items
					mod <- NULL
					
					##   Loop through the item response models
					for (j in 1:length(pars[[grp]]@model)) {
						##   Identify the items for the given model that
						##   are part of the subset of items selected
						tmp <- pars[[grp]]@items[[j]][pars[[grp]]@items[[j]] %in% items1]
						
						##   Extract the incremented item numbers
						tmp <- items1a[items1a[,2]%in%tmp,1]
						
						##   If there are any items in the subset associated with 
						##   this model, add the model to {mod}
						if (length(tmp)) mod <- c(mod,pars[[grp]]@model[j])
						
						##   Update the indicator variable
						m.type[tmp,][!is.na(m.type[tmp,])] <- j
					}
					
					##   Compile the slopes for each group and the indicator variable
					pa <- data.frame(pa1=pa1[!is.na(pa1)],pa2=pa2[!is.na(pa2)],type=factor(m.type[!is.na(m.type)],labels=toupper(mod)))
					
					##   Compute a range of 3 SDs for creating a confidence interval
					sd.a <- .sd(pa[,1]-pa[,2])
					p.min <- min(pa[,1:2])-0.25
					p.max <- max(pa[,1:2])+0.25
					
					##   Create the plot with different markers for each response model
					if (sep.mod==TRUE) {
						pl.out[[i]][[pl]] <- xyplot(pa2~pa1,data=pa,groups=type, auto.key=list(space="bottom", columns=length(mod)),
						,xlab=nms[grp], ylab=nms[grp+1], main=a.title, panel=pfunc, xlim=c(p.min,p.max) ,ylim=c(p.min,p.max),sd=sd.a, item.names=item.names, 
						drift.sd=drift.sd)
					
					##   Create the plot with a single marker for all response models
					} else {
						pl.out[[i]][[pl]] <- xyplot(pa2~pa1,data=pa,xlab=nms[grp], ylab=nms[grp+1], main=a.title, panel=pfunc, xlim=c(p.min,p.max) ,ylim=c(p.min,p.max), sd=sd.a, item.names=item.names, drift.sd=drift.sd)
					}
					
					##   Increment the plot index
					pl <- pl+1
				}
				
				
				##   Compare difficulty/threshold/step/category parameters
				if ("b" %in% drift) {
				
					##   Extract difficulty/threshold/step/category parameters for the lower group
					pb1 <- pars[[grp]]@b[items1,]
					
					##   Extract difficulty/threshold/step/category parameters for the higher group
					pb2 <- pars[[grp+1]]@b[items2,]
					
					##   Transpose these values if there is only one item
					if (is.vector(pb1)) {
						if (length(items1)==1) {
							pb1 <- t(pb1)
							pb2 <- t(pb2)
						} else {
							pb1 <- as.matrix(pb1)
							pb2 <- as.matrix(pb2)
						}
					}
					
					##   Create a title
					b.title <- NULL
					if ("drm" %in% pars[[grp]]@model) b.title <- c(b.title,"Difficulty")
					if ("grm" %in% pars[[grp]]@model) b.title <- c(b.title,"Threshold")
					if ("gpcm" %in% pars[[grp]]@model) b.title <- c(b.title,"Step")
					if ("nrm" %in% pars[[grp]]@model|"mcm" %in% pars[[grp]]@model) b.title <- c(b.title,"Category")
					b.title <- paste(paste(b.title,collapse="/"),"Parameters")
					
					##   Initialize an indicator variable to identify the item
					##   response model associated with each common item
					m.type <- pb1
					
					##   Create an object that includes an incremented item
					##   number alongside the item numbers for the selected subset
					items1a=cbind(1:length(items1),items1)
					
					##   Initialize an object to store the response models for the subset of items
					mod <- NULL
					
					##   Loop through the item response models
					for (j in 1:length(pars[[grp]]@model)) {
						##   Identify the items for the given model that
						##   are part of the subset of items selected
						tmp <- pars[[grp]]@items[[j]][pars[[grp]]@items[[j]] %in% items1]
						
						##   Extract the incremented item numbers
						tmp <- items1a[items1a[,2]%in%tmp,1]
						
						##   If there are any items in the subset associated with 
						##   this model, add the model to {mod}
						if (length(tmp)) mod <- c(mod,pars[[grp]]@model[j])
						
						##   Update the indicator variable
						m.type[tmp,][!is.na(m.type[tmp,])] <- j
					}
					
					##   Compile the slopes for each group and the indicator variable
					pb <- data.frame(pb1=pb1[!is.na(pb1)],pb2=pb2[!is.na(pb2)],type=factor(m.type[!is.na(m.type)],labels=toupper(mod)))
					
					##   Compute a range of 3 SDs for creating a confidence interval
					sd.b <- .sd(pb[,1]-pb[,2])
					p.min <- min(pb[,1:2])-0.25
					p.max <- max(pb[,1:2])+0.25
					
					##   Create the plot with different markers for each response model
					if (sep.mod==TRUE) {
						pl.out[[i]][[pl]] <- xyplot(pb2~pb1,data=pb,groups=type,auto.key=list(space="bottom", columns=length(mod)), 
						xlab=nms[grp], ylab=nms[grp+1], main=b.title, panel=pfunc, xlim=c(p.min,p.max) ,ylim=c(p.min,p.max), sd=sd.b, item.names=item.names, drift.sd=drift.sd)
					
					##   Create the plot with a single marker for all response models
					} else {
						pl.out[[i]][[pl]] <- xyplot(pb2~pb1,data=pb,xlab=nms[grp], ylab=nms[grp+1], main=b.title, panel=pfunc, xlim=c(p.min,p.max) ,ylim=c(p.min,p.max), sd=sd.b, item.names=item.names, drift.sd=drift.sd)
					}
					
					##   Increment the plot index
					pl <- pl+1
				} 
				
				
				##   Compare lower asymptote parameters
				if ("c" %in% drift) {
				
					##   Extract asymptotes for the lower group
					pc1 <- pars[[grp]]@c[items1,]
					
					##   Extract asymptotes for the higher group
					pc2 <- pars[[grp+1]]@c[items2,]
					
					##   Check to see if there are any asymptote values 
					##   for any of the common items
					if (sum(pc1,na.rm=T)>0) {
					
						##   Transpose the extracted values if there is only one item
						if (is.vector(pc1)) {
							if (length(items1)==1) {
								pc1 <- t(pc1)
								pc2 <- t(pc2)
							} else {
								pc1 <- as.matrix(pc1)
								pc2 <- as.matrix(pc2)
							}
						}
						
						##   Create a title
						c.title <- "Lower Asymptote Parameters"
						
						##   Initialize an indicator variable to identify the item
						##   response model associated with each common item
						m.type <- pc1
						
						##   Create an object that includes an incremented item
						##   number alongside the item numbers for the selected subset
						items1a=cbind(1:length(items1),items1)
						
						##   Initialize an object to store the response models for the subset of items
						mod <- NULL
						
						##   Loop through the item response models
						for (j in 1:length(pars[[grp]]@model)) {
							##   Identify the items for the given model that
							##   are part of the subset of items selected
							tmp <- pars[[grp]]@items[[j]][pars[[grp]]@items[[j]] %in% items1]
							
							##   Extract the incremented item numbers
							tmp <- items1a[items1a[,2]%in%tmp,1]
							
							##   If there are any items in the subset associated with 
							##   this model, add the model to {mod}
							if (length(tmp)) mod <- c(mod,pars[[grp]]@model[j])
							
							##   Update the indicator variable
							m.type[tmp,][!is.na(m.type[tmp,])] <- j
						}
						
						##   Compile the slopes for each group and the indicator variable
						pc <- data.frame(pc1=pc1[!is.na(pc1)],pc2=pc2[!is.na(pc2)],type=factor(m.type[!is.na(m.type)], labels=toupper(mod[mod %in% c("drm","mcm")])))
						
						##   Compute a range of 3 SDs for creating a confidence interval
						sd.c <- .sd(pc[,1]-pc[,2])
						p.max <- max(pc[,1:2])+0.1
						
						##   Create the plot with different markers for each response model
						if (sep.mod==TRUE) {
							pl.out[[i]][[pl]] <- xyplot(pc2~pc1,data=pc,groups=type,auto.key=list(space="bottom", columns=length(mod[mod %in% c("drm","mcm")])),
							xlab=nms[grp], ylab=nms[grp+1], main=c.title, panel=pfunc, xlim=c(0,p.max) ,ylim=c(0,p.max),sd=sd.c, item.names=item.names, drift.sd=drift.sd)
							
						##   Create the plot with a single marker for all response models
						} else {
							pl.out[[i]][[pl]] <- xyplot(pc2~pc1,data=pc,xlab=nms[grp], ylab=nms[grp+1], main=c.title, panel=pfunc, xlim=c(0,p.max) ,ylim=c(0,p.max),sd=sd.c, item.names=item.names, drift.sd=drift.sd)
						}
						
						##   Increment the plot index
						pl <- pl+1
					}
				} 
				
				if ("ICC" %in% drift|"TCC" %in% drift) {
					
					##   Reorder the common items
					sort <- order(items1)
					items1 <- items1[sort]
					items2 <- items2[sort]
					
					##   Create a duplicate object of response probabilities
					##   The probability matrix for each group will be modified
					##   to hold only the columns for the common items
					prob1 <- prob
					
					##   Identify the number of columns in the matrix of
					##   probabilities associated with each item in the lower group
					tmp <- rep(1:length(prob1[[grp]]@p.cat),prob1[[grp]]@p.cat)
					
					##   Identify the columns of probabilities for the common
					##   items in the lower group
					tmp1 <- 1:length(tmp)
					tmp2 <- NULL
					for (j in items1) {
						tmp2 <- c(tmp2, tmp1[tmp==j])
					}
					
					##   Re-specify the matrix of probabilities to include 
					##   the column(s) of theta values and the appropriate 
					##   columns of probabilities for the common items
					##   in the lower group
					prob1[[grp]]@prob <- prob1[[grp]]@prob[,c(1:prob1[[grp]]@dimensions,tmp2+prob1[[grp]]@dimensions)]
					
					##   Re-specify the {p.cat} vector to identify the number of 
					##   columns in the probability matrix associated with each 
					##   common item in the lower group
					prob1[[grp]]@p.cat <- prob1[[grp]]@p.cat[items1]
				
					##   Identify the number of columns in the matrix of
					##   probabilities associated with each item in the higher group
					tmp <- rep(1:length(prob1[[grp+1]]@p.cat),prob1[[grp+1]]@p.cat)
					
					##   Identify the columns of probabilities for the common
					##   items in the higher group
					tmp1 <- 1:length(tmp)
					tmp2 <- NULL
					for (j in items2) {
						tmp2 <- c(tmp2, tmp1[tmp==j])
					}
					
					##   Re-specify the matrix of probabilities to include 
					##   the column(s) of theta values and the appropriate 
					##   columns of probabilities for the common items
					##   in the higher group
					prob1[[grp+1]]@prob <- prob1[[grp+1]]@prob[,c(1:prob1[[grp+1]]@dimensions,tmp2+prob1[[grp+1]]@dimensions)]
					
					##   Re-specify the {p.cat} vector to identify the number of 
					##   columns in the probability matrix associated with each 
					##   common item in the higher group
					prob1[[grp+1]]@p.cat <- prob1[[grp+1]]@p.cat[items2]
					
					if ("TCC" %in% drift) {
						##   Compute a vector of response weights
						##   for the lower and higher groups
						scr1 <- scr2 <- NULL
						for (j in 1:length(prob1[[grp]]@p.cat)) {
							scr1 <- c(scr1,1:prob1[[grp]]@p.cat[j])
						}
						for (j in 1:length(prob1[[grp+1]]@p.cat)) {
							scr2 <- c(scr2,1:prob1[[grp+1]]@p.cat[j])
						}
						
						p1 <- as.matrix(prob1[[grp]]@prob[,-c(1:prob1[[grp]]@dimensions)])
						p1 <- p1%*%scr1
						p2 <- as.matrix(prob1[[grp+1]]@prob[,-c(1:prob1[[grp+1]]@dimensions)])
						p2 <- p2%*%scr2
						
						if (prob1[[grp]]@dimensions==1) {
						
							th <- c(prob1[[grp]]@prob[,1],prob1[[grp+1]]@prob[,1])
							p.all <- data.frame(cbind(th,c(p1,p2),c(rep(1,length(p1)),rep(2,length(p2)))))
							names(p.all) <- c("theta","prob","gr")
							
							pl.out[[i]][[pl]] <- xyplot(prob~theta,data=p.all,type="l",
							as.table=TRUE,
							ylab="True Score",
							xlab=expression(paste(theta)),
							groups=p.all$gr,
							col=1:2, lty=1:2,
							key=list(space="bottom", text=list(c(nms[grp],nms[grp+1])),lines=list(col=1:2, lty=1:2), columns=2),...)
							
						} else {
							##   TCCs for multidimenisonal data will be implemented in a future release
						}
						
						##   Increment the plot index
						pl <- pl+1
						
					}
					
					if ("ICC" %in% drift) {
						
						dimensions <- prob1[[grp]]@dimensions
						
						if(!missing(separate)) {
							##   Re-specify {cat} so that each column in the matrix of 
							##   probabilities will be plotted in a separate panel
							if (separate==TRUE) cat <- rep(1,ncol(prob1[[grp]]@prob)-dimensions) else cat <- prob1[[grp]]@p.cat
						} else {
							separate <- FALSE
							cat <- prob1[[grp]]@p.cat
						}
						
						##   Extract the first column(s) in {prob} containing the 
						##   theta values used to compute the probabilities
						theta <- as.matrix(prob1[[grp]]@prob[,1:dimensions])
						
						##   Number of (combinations) of theta values
						nt <- nrow(theta)
						
						##   Number of "items"
						ni <- length(cat)
						
						##   If there is more than one item, create a matrix where 
						##   the probabilities for each item and category are stacked
						##   on top of each other. This is necessary
						##   to use the group argument in the various lattice plots.
						if (ncol(prob1[[grp]]@prob)>(dimensions+1)) {
							sx <- stack(prob1[[grp]]@prob[,-c(1:dimensions)]) 
						} else {
							sx <- data.frame(values=prob1[[grp]]@prob[,-c(1:dimensions)],ind=factor(rep(colnames(prob1[[grp]]@prob)[dimensions+1],nt)))
						}
						
						if (ncol(prob1[[grp+1]]@prob)>(dimensions+1)) {
							sx1 <- stack(prob1[[grp+1]]@prob[,-c(1:dimensions)]) 
						} else {
							sx1 <- data.frame(values=prob1[[grp+1]]@prob[,-c(1:dimensions)],ind=factor(rep(colnames(prob1[[grp+1]]@prob)[dimensions+1],nt)))
						}
						sx <-c(sx$values,sx1$values)
						
						##   Create indicator values for the item numbers and category numbers
						##   to identify groups in the sx matrix created above
						id <- NULL
						cid <- NULL
						
						for (j in 1:ni) {
							id <- c(id, rep(j,nt*cat[j]))
							for (k in 1:cat[j]) {
								cid <- c(cid, rep(k,nt))
							}
						}
						
						
						if (separate==TRUE) {
						
							##   Initialize an object to store the item names
							nms <- NULL
							
							##   Loop through all of the items
							for (j in 1:length(items)) {
							
								##   For dichotomous items, only include the item number
								if (prob1[[grp]]@p.cat[j]==1) {
									nms <- c(nms,paste("Item",items[j]))
									
								##   For polytomous items, include both the item
								##   number and the category number
								} else {
									for (k in 1:x@p.cat[j]) {
										nms <- c(nms, paste("Item ",items[j],".",j,sep=""))
									}
								}
							}
							id <- factor(id,seq(1:ni),nms)
						} else {
							id <- c(id,id)
						}
						
						##   Create a set of theta values to be attached to the sx
						##   matrix created above
						if (sum(cat)>1) {
							tmp <- theta
							for (j in 2:sum(cat)) {
								tmp <- rbind(tmp,theta)
							}
							theta <- tmp
						}
						
						test <- c(rep(1,nrow(theta)),rep(2,nrow(theta)))
						if (dimensions==1) theta <- c(theta,theta) else theta <- rbind(theta,theta)
						
						##   Combine the theta values, response probabilities,
						##   item indicators, and category indicators into a 
						##   single matrix
						out <- data.frame(cbind(theta,id,sx,c(cid,cid),test))
						names(out) <- c(paste("theta",1:dimensions,sep=""),"id","values","cid","test")
						
						##   Re-specify id as a factor
						out$id <- factor(out$id,seq(1:length(items1)),paste("Item ",items1,",",items2,sep=""))
						
						##   Create the formula that will be passed to the given lattice function
						if (dimensions>1) {
							tmp <- paste(names(out)[1:2],collapse="+")
							if (dimensions>2) {
								for (j in 1:(dimensions-2)) {
									out[,j+2] <- as.factor(out[,j+2])
									
								}
								tmp2 <- paste(c(names(out)[3:dimensions],"id"),collapse="+") 
							} else {
								tmp2 <- "id"
							}
							form<- as.formula(paste("values~",tmp,"|",tmp2,sep="")) 
						}
						
						##   Custom strip
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
						
						##   Determine the number of panels to print
						if (dimensions<3) np <- ni else np <- ni*length(unique(out[,1]))^(dimensions-2)
						if (missing(panels)) {
							if (np<20) panels <- np else panels <- 20
						}
						
						##   Initialize a set of colors for the strips for multiple dimensions
						cols <- c("lightpink1", "darkseagreen1", "burlywood1", "cadetblue3", "yellow", "darkorchid2", "coral", "seagreen4")
						if (dimensions<3) cols <- "lightblue" else cols <- c(cols[(dimensions-2):1],"lightblue")
					
						if (dimensions==1) {
							if (np>panels) {
							
								##   Check to see if a graphics device is already open
								if (names(dev.cur())=="null device") {
									if (Sys.info()["sysname"] == "Windows") {
										cat("Use PgUp and PgDn to view different plot pages\n")
										windows(record=TRUE)
									}
								}
								
								##   Create the multi-page plot
								pl.out[[i]][[pl]] <- xyplot(values~theta1|id,out,type="l",
								as.table=TRUE,
								ylab="Probability",
								xlab=expression(paste(theta)),
								groups=interaction(cid,test),
								par.strip.text=list(cex=0.7),
								par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1), superpose.line=list(col=c("red","blue")) ) ,
								strip=strip.irt.prob,
								layout=c(0,panels),
								key=list(space="bottom", text=list(c(nms[grp],nms[grp+1])),lines=list(col=1:2, lty=1:2), columns=2),
								ylim=c(-0.05,1.05),...)
								
							##   Create a single-page plot
							} else {
								pl.out[[i]][[pl]] <- xyplot(values~theta1|id,out,type="l",
								as.table=TRUE,
								ylab="Probability",
								xlab=expression(paste(theta)),
								groups=interaction(cid,test),
								par.strip.text=list(cex=0.7),
								par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1), superpose.line=list(col=c("red","blue")) ) ,
								strip=strip.irt.prob,
								key=list(space="bottom", text=list(c(nms[grp],nms[grp+1])),lines=list(col=1:2, lty=1:2), columns=2),
								ylim=c(-0.05,1.05),...)
							}
						} else if (dimensions>1) {
							##   ICCs for multidimenisonal data will be implemented in a future release
						}
						
						##   Increment the plot index
						pl <- pl+1
						
					}
					
				}
			}
				
			
			##   Print all of the comparison plots
			for (i in 1:length(pl.out)) {
				##   Include functionality for d.split here later
				for (j in 1:length(pl.out[[i]])) {
					print(pl.out[[i]][[j]])
				}
			}
		}
	}
}



plot.irt.prob <- function(x, y, ..., type, separate, combine, items, item.names, item.nums, panels, save.hist) {
	
	if (missing(type)) type <- "wireframe"
	dots <- list(...)
	
	##   Create vectorplots or information plots (the information plots are not currently implemented)
	if (type %in% c("vectorplot1","vectorplot2","vectorplot3","information1","information2")) {
	
		##   Determine the number of response categories for each item
		cat <- x@p.cat
		cat[cat==1] <- 2
		cat[cat>2] <- cat[cat>2]-1
		
		##   Create an {irt.pars} object
		tmp <- as.irt.pars(x@pars,cat=cat,poly.mod=as.poly.mod(length(cat),x@model,x@items),dimensions=x@dimensions)
		
		##   Create the plot
		plot(tmp,type=type,items=items, item.names=item.names, item.nums=item.nums,panels=panels,...)
		
	##   Create plots of item response curves/surfaces
	} else {
	
		##   Check to see if there is a graphics device already open
		if (dev.cur()=="null device") {
			if (!missing(save.hist)) {
				if (save.hist==FALSE) {
					##   If not, clear any saved plot history
					if (exists(".SavedPlots",where=1)) rm(".SavedPlots",pos=1)
				}
			}
		}
		
		##   Number of dimensions
		dimensions <- x@dimensions
		
		if (missing(item.nums)) item.nums <- TRUE
		
		##   Use all items
		if(missing(items)) {
			items <- 1:length(x@p.cat)
			
		##   Use a subset of items
		} else {
			tmp <- rep(1:length(x@p.cat),x@p.cat)
			tmp1 <- 1:length(tmp)
			tmp1 <- tmp1[tmp %in% items]
			
			##   Extract the response probabilities and numbers of 
			##   categories for the given subset of items
			x@prob <- x@prob[,c(1:dimensions,tmp1+dimensions)]
			x@p.cat <- x@p.cat[items]
		}
		
		
		if(!missing(separate)) {
			##   Re-specify {cat} so that each column in the matrix of 
			##   probabilities will be plotted in a separate panel
			if (separate==TRUE) cat <- rep(1,ncol(x@prob)-x@dimensions) else cat <- x@p.cat
		} else {
			separate <- FALSE
		}
		
		
		if (!missing(combine)) {
		
			##   Identify the number of categories (i.e., columns in {prob}) 
			##   to combine in each panel
			cat <- combine 
			
			
			##   Check to see if the categories specified by combine
			##   encapsulate all of the responses for all items
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
		
		##   Extract the first column(s) in {prob} containing the 
		##   theta values used to compute the probabilities
		theta <- as.matrix(x@prob[,1:dimensions])
		
		##   Number of (combinations) of theta values
		nt <- nrow(theta)
		
		##   Number of "items" or combined categories
		ni <- length(cat)
		
		##   If there is more than one item, create a matrix where 
		##   the probabilities for each item and category are stacked
		##   on top of each other. This is necessary
		##   to use the group argument in the various lattice plots.
		if (ncol(x@prob)>(dimensions+1)) {
			sx <- stack(x@prob[,-c(1:dimensions)]) 
		} else {
			sx <- data.frame(values=x@prob[,-c(1:dimensions)],ind=factor(rep(colnames(x@prob)[dimensions+1],nt)))
		}
		
		##   Create indicator values for the item numbers and category numbers
		##   to identify groups in the sx matrix created above
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
			
				##   I still need to figure out how to create
				##   item names when response categories are combined
			} else {
				if (separate==TRUE) {
				
					##   Initialize an object to store the item names
					nms <- NULL
					
					##   Loop through all of the items
					for (i in 1:length(items)) {
					
						##   For dichotomous items, only include the item number
						if (x@p.cat[i]==1) {
							nms <- c(nms,paste("Item",items[i]))
							
						##   For polytomous items, include both the item
						##   number and the category number
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
		
		##   Create a set of theta values to be attached to the sx
		##   matrix created above
		if (sum(cat)>1) {
			tmp <- theta
			for (i in 2:sum(cat)) {
				tmp <- rbind(tmp,theta)
			}
			theta <- tmp
		}
		
		##   Combine the theta values, response probabilities,
		##   item indicators, and category indicators into a 
		##   single matrix
		out <- cbind(theta,id,cid,sx)
		colnames(out)[1:dimensions] <- paste("theta",1:dimensions,sep="")
		
		##   Create the formula that will be passed to the given lattice function
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
		
		##   Custom strip
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
		
		##   Determine the number of panels to print
		if (dimensions<3) np <- ni else np <- ni*length(unique(out[,1]))^(dimensions-2)
		if (missing(panels)) {
			if (np<20) panels <- np else panels <- 20
		}
		
		##   Initialize a set of colors for the strips for multiple dimensions
		cols <- c("lightpink1", "darkseagreen1", "burlywood1", "cadetblue3", "yellow", "darkorchid2", "coral", "seagreen4")
		if (dimensions<3) cols <- "lightblue" else cols <- c(cols[(dimensions-2):1],"lightblue")
	
		if (dimensions==1) {
			if (np>panels) {
			
				##   Check to see if a graphics device is already open
				if (names(dev.cur())=="null device") {
					if (Sys.info()["sysname"] == "Windows") {
						cat("Use PgUp and PgDn to view different plot pages\n")
						windows(record=TRUE)
					}
				}
				
				##   Create the multi-page plot
				pl.out <- xyplot(values~theta1|id,out,type="l",
				as.table=TRUE,
				ylab="Probability",
				xlab=expression(paste(theta)),
				groups=cid,
				par.strip.text=list(cex=0.7),
				par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1) ) ,
				strip=strip.irt.prob,
				layout=c(0,panels),
				ylim=c(-0.05,1.05),...)
				
			##   Create a single-page plot
			} else {
				pl.out <- xyplot(values~theta1|id,out,type="l",
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
			if (np>panels) {
				
				##   Check to see if a graphics device is already open
				if (names(dev.cur())=="null device") {
					if (Sys.info()["sysname"] == "Windows") {
						cat("Use PgUp and PgDn to view different plot pages\n")
						windows(record=TRUE)
					}
				}
				
				##   Create the multi-page plot
				mlab <- "Probability"
				pl.out <- eval(parse(text=paste(type,"(form,out, as.table=TRUE, 
				zlab=list(label=mlab, rot=90),
				xlab=expression(paste(theta[1])),
				ylab=expression(paste(theta[2])),
				groups=cid,
				par.strip.text=list(cex=0.7),
				par.settings = list(strip.background = list(col = cols), layout.heights=list(strip=1.1) ) ,
				strip=strip.irt.prob,
				layout=c(0,panels),
				zlim=c(0,1),...)",sep="")))
				
			##   Create a single-page plot
			} else {
				mlab <- "Probability"
				pl.out <- eval(parse(text=paste(type,"(form,out, as.table=TRUE, 
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
		
		if (length(dots$mg)) {
			##   Return the trellis object for later printing
			return(pl.out)
		} else {
			##   Print the plot
			print(pl.out)
		}
	}
}


plot.sep.pars <- function(x, y, ..., type, separate, combine, items, item.names, item.nums, panels, save.hist) {
	
	x <- as.irt.pars(x)
	plot(x,type=type,separate=separate,combine=combine, items=items, item.names=item.names, item.nums=item.nums,panels=panels, ...)
	
}

plot.list <- function(x, y, ..., type, separate, combine, items, item.names, item.nums, panels, drift, groups, grp.names, sep.mod, drift.sd, save.hist) {

	##   Extract the rescaled item parameters from the object output by {plink}
	if (length(x$pars)) {
		tmp <- x$pars
		
		if (missing(drift.sd)) drift.sd <- 3
		plot(tmp,type=type,separate=separate,combine=combine, items=items, item.names=item.names, item.nums=item.nums,panels=panels,drift=drift,groups=groups, grp.names=grp.names,sep.mod=sep.mod, drift.sd=drift.sd, ...)
		
	##   When x is a list of irt.prob objects
	} else if (class(x[[1]])=="irt.prob") {
		
		##   Number of groups
		ng <- length(x)
		
		pl.out <- vector("list",ng)
		
		for (i in 1:ng) {
			if (ng==1) {
				plot(x[[i]],type=type,separate=separate,combine=combine, items=items, item.names=item.names, item.nums=item.nums,panels=panels,...)
			} else {
				grn.items <- unique(unlist(strsplit(names(x[[i]]@prob)[-c(1:x[[i]]@dimensions)],"\\.")))
				tmp <- suppressWarnings(as.numeric(grn.items))
				grn.items <- grn.items[is.na(tmp)]
				if (missing(grp.names)) {
					nms <- paste("G",i,": ",grn.items,sep="")
				} else {
					nms <- paste(grp.names[i],": ",grn.items,sep="")
				}
				pl.out[[i]] <- plot(x[[i]],type=type,separate=separate,combine=combine, items=items, item.names=nms, item.nums=item.nums,panels=panels,mg=1,...)
			}
		}
	
		##   Print compiled trellis objects
		if (length(pl.out)) {
			for (i in 1:ng) {
				print(pl.out[[i]])
			}
		}
		
	} else {
		stop("There were no parameters in {x}, re-run plink and specify and argument for {rescale} then try again.")
	}
}