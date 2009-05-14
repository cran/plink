read.bilog <- function(file, ability=FALSE, pars.only=TRUE, as.irt.pars=FALSE) {
	if (ability==FALSE) {
		pars <- read.fwf(file, c(8,8,rep(10,13),4,1,1),sep="\t", skip=4)
		pars <- pars[,-15]
		colnames(pars) <- c("item", "subtest", "intercept","int.se", "slope", "slope.se", "threshold", "thresh.se", "dispersion", "disp.se", "asymptote", "asymp.se", "drift", "drift.se","stream.loc","key","dummy.vals")
		if (pars.only==TRUE) pars <- pars[,c(5,7,11)]
		
		if (as.irt.pars==TRUE) {
			if (pars.only==FALSE) pars <- pars[,c(5,7,11)]
			n <- nrow(pars)
			pm <- as.poly.mod(n)
			pars <- as.irt.pars(pars, cat=rep(2,n), poly.mod=pm)
		}
	} else {
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

read.parscale <- function(file, ability=FALSE,  loc.out=FALSE, pars.only=TRUE, as.irt.pars=FALSE) {
	if (ability==FALSE) {
		nums <- as.numeric(read.fwf(file, c(8,rep(5,5)), skip=2, n=1)[-1])
		block <- as.numeric(read.fwf(file, rep(5,nums[1]), skip=3, n=1))
		pars <- NULL
		pars1 <- vector("list",length(block))
		cat <- NULL
		skip <- 5
		skip1 <- skip+block[1]
		for (i in 1:length(block)) {
			if (i>1) {
				skip <- 5 + sum(block[1:(i-1)]) + (i-1)*2
				skip1 <- 5 + sum(block[1:i]) + (i-1)*2
			}
			tmp <- read.fwf(file, c(8,5,4,rep(10,6)), skip=skip, n=block[i], colClasses=c("character","numeric","character",rep("numeric",6)))
			pars <- rbind(pars, tmp)
			if (is.data.frame(tmp)) {
				cat <- c(cat, as.numeric(tmp[,2]))
			} else {
				cat <- c(cat, as.numeric(tmp[2]))
			}
			tmp1 <- read.fwf(file, rep(10,cat[length(cat)]), skip=skip1, n=2)
			if (nums[3]<5) {
				tmp1 <- tmp1[,-ncol(tmp1)] # Graded response models
			} else {
				tmp1 <- tmp1[,-1] # Partial credit models
			}
			pars1[[i]] <- tmp1
		}
		colnames(pars) <- c("block.name","cat","item","slope","slope.se","location","loc.se","asymptote","asymp.se")
		
		if (max(cat)>2) {
			step <- matrix(NA, nrow(pars), 2*(max(cat)-1))
			pars1 <- rep(pars1,block)
			for (i in 1:length(pars1)) {
				if (is.data.frame(pars1[[i]])) {
					step[i,] <- unlist(pars1[[i]])
				}
			}
			
			if (nums[3]<5) prefix <- "thresh" else prefix <- "step"
			
			step.names <- character()
			for (j in 1:(max(cat)-1)) {
				step.names <- c(step.names, paste(prefix,j,sep=""), paste(prefix,j,".se",sep=""))
			}
			colnames(step) <- step.names
			pars <- cbind(pars, step)
		}
		pars <- data.frame(pars)
		
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
		
		if (pars.only==TRUE) pars <- pars[,seq(4,ncol(pars),2)]
		
		if (as.irt.pars==TRUE) {
			if (pars.only==FALSE) pars <- pars[,seq(4,ncol(pars),2)]
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
			
			n <- nrow(pars)
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

read.multilog <- function(file, cat, poly.mod, ability=FALSE, contrast="dev", drm.3PL=TRUE, loc.out=FALSE, as.irt.pars=FALSE) {
	if (ability==FALSE) {
		# Import and reformat item parameters as a vector
		pars <- read.fwf(file, rep(12,30))
		tmp <- apply(is.na(pars),2,sum)
		pars <- pars[,tmp<nrow(pars)]
		if (sum(abs(pars[nrow(pars),1:3])-c(1,0,1))==0) pars <- pars[-nrow(pars),]
		pars <- as.vector(t(as.matrix(pars)))
		pars <- pars[!is.na(pars)]
		
		# Extract the item parameters for each item
		mod <- poly.mod@model
		items <- poly.mod@items
		p.cat <- cat
		for (i in 1:length(mod)) {
			if (mod[i]=="drm") {
				if (drm.3PL==TRUE) p.cat[items[[i]]] <- 4
			} else if (mod[i] %in% c("gpcm","nrm","mcm")) {
				p.cat[items[[i]]] <- (p.cat[items[[i]]]-1)*3
			}
		}
		p <- vector("list", length(cat))
		k <- 1
		for (i in 1:length(cat)) {
			p[[i]] <- pars[k:(k-1+p.cat[i])]
			k <- k+p.cat[i]
		}
		
		# Determine the maximum number of response categories for NRM and MCM items
		# for use in compiling the file matrix of item parameters
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
		col.max <- a.max+c.max+d.max
		
		# Create polynomial contrasts
		poly <- list(c(-.5,.5),-1:1, -1.5:1.5, -2:2, -2.5:2.5, -3:3, -3.5:3.5, -4:4)
			
		# Create contrast lists
		n <- length(cat)
		if (is.character(contrast)) {
			# Formulation One
			con <- list(dev=NULL,poly=NULL,tri=NULL)
			contrast <- tolower(contrast)
			eval(parse(text=paste("con$",contrast,"=1:",n,sep="")))
			con <- rep(con,3)
		} else {
			if (length(con)!=9) stop("The object {constant} must be a list of length nine")
			con <- contrast
			if (is.character(unlist(con))) { # Formulation Two
				for (i in 1:length(con)) {
					tmp <- NULL
					for (j in 1:length(mod)) {
						if (mod[j] %in% con[[i]]) tmp <- c(tmp, items[[j]])
					}
					con[[i]] <- tmp
				}
			} else { # Formulation Three
				if (!is.numeric(unlist(con))) stop("Under formulation three, all the values must be numeric")
				ak <- con[1:3]
				tmp <- unlist(ak)
				tmp1 <- c(1:n)%in%tmp
				if (sum(tmp1)!=n) ak[[1]] <- c(ak[[1]], c(1:n)[tmp1==FALSE])
				
				ck <- con[4:6]
				tmp <- unlist(ck)
				tmp1 <- c(1:n)%in%tmp
				if (sum(tmp1)!=n) ck[[1]] <- c(ck[[1]], c(1:n)[tmp1==FALSE])
				
				dk <- con[7:9]
				tmp <- unlist(dk)
				tmp1 <- c(1:n)%in%tmp
				if (sum(tmp1)!=n) dk[[1]] <- c(dk[[1]], c(1:n)[tmp1==FALSE])
				
				con <- c(ak,ck,dk)
			}
		}
		names(con) <- c("dev.a","poly.a","tri.a","dev.c","poly.c","tri.c","dev.d","poly.d","tri.d")
		
		# Reformat the contrast parameters to traditional IRT parameters
		for (i in 1:length(mod)) {
			for (j in items[[i]]) {
				if (mod[[i]]=="drm") {
					if (drm.3PL==TRUE) {
						p[[j]] <- p[[j]][1:3]
						p[[j]][2] <- p[[j]][2]/-p[[j]][1]
						p[[j]][1] <- p[[j]][1]/1.7
						if (j %in% con$tri.c) {
							p[[j]][3] <- exp(-p[[j]][3])/(1+exp(-p[[j]][3]))
						} else {
							p[[j]][3] <- exp(p[[j]][3])/(1+exp(p[[j]][3]))
						}
					}
					p[[j]] <- c(p[[j]], rep(NA,col.max-length(p[[j]])))
				} else if (mod[[i]]=="grm") {
					if (loc.out==TRUE) {
						ck <- p[[j]][-1]
						ck <- c(mean(ck),ck-mean(ck))
						p[[j]] <- c(p[[j]][1],ck)
					}
					p[[j]] <- c(p[[j]], rep(NA,col.max-length(p[[j]])))
				} else if (mod[[i]] %in% c("gpcm","nrm","mcm")) {
					tmp <- p.cat[j]/3
					ak <- p[[j]][1:tmp]
					ck <- p[[j]][(tmp+1):(tmp*2)]
					dk <- p[[j]][(2*tmp+1):(3*tmp-1)]
					
					# Recode discrimination
					if (j %in% con$dev.a) {
						C <- rep(-1/(tmp+1),tmp)
						tmp.a <- ak%*%C
						ak <- c(tmp.a, tmp.a+ak)
					} else if (j %in% con$poly.a) {
						C <- poly[[tmp-1]]
						ak <- ak%*%C
					} else if (j %in% con$tri.a) {
						C <- c(0,rep(-1,tmp-1))
						ak <- ak%*%C
					}
					
					# Recode step parameters
					if (j %in% con$dev.c) {
						TC<- rep(-1/(tmp+1),tmp)
						tmp.c <- ck%*%C
						ck <- c(tmp.c, tmp.c+ck)
					} else if (j %in% con$poly.c) {
						C <- poly[[tmp-1]]
						ck <- ck%*%C
					} else if (j %in% con$tri.c) {
						C <- c(0,rep(-1,tmp-1))
						ck <- ck%*%C
					}
					
					# Recode lower asymptote parameters
					if (j %in% con$dev.d) {
						C <- rep(-1/(tmp+1),tmp-1)
						tmp.d <- dk%*%C
						dk <- c(tmp.d, tmp.d+dk)
					} else if (j %in% con$poly.d) {
						C <- poly[[tmp-2]]
						dk <- dk%*%C
					} else if (j %in% con$tri.d) {
						C <- c(0,rep(-1,tmp-2))
						dk <- dk%*%C
					}
					num <- exp(dk)
					dk <- num/sum(num)
					
					if (mod[[i]]=="gpcm") {
						ak <- ak[length(ak)]
						if (loc.out==TRUE) {
							ck <- c(mean(ck),ck-mean(ck))
						}
						p[[j]] <- c(ak, ck)
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
					p[[j]] <- c(p[[j]], rep(NA,col.max-length(p[[j]])))
				} 
			}
		}
		pars <- NULL
		for (i in 1:length(p)){
			pars <- rbind(pars,p[[i]])
		}
		if (as.irt.pars==TRUE) {
			pars <- as.irt.pars(pars, cat=cat, poly.mod=poly.mod, location=loc.out)
		}
	} else {
		# Import ability parameters
		pars <- read.fwf(file, c(10,10,4))
		names(pars) <- c("theta","theta.se","freq")
	}
	return(pars)
}

read.testfact <- function(file, ability=FALSE, guessing=FALSE, bifactor=FALSE, as.irt.pars=FALSE) {
	if (ability==FALSE) {
		tmp1 <- scan(file, skip=1, what="character", quiet=TRUE, nlines=1)
		tmp2 <- scan(file, skip=2, what="character", quiet=TRUE, nlines=1)
		if (tmp2[1]!="2") tmp1 <- c(tmp1,tmp2)
		pars <- scan(file, skip=1, what="character", quiet=TRUE)
		items <- length(pars)/length(tmp1)
		pars <- matrix(pars, items, length(tmp1), byrow=TRUE)
		pars <- pars[,-c(1,2)]
		
		if (guessing==TRUE) {
			dimensions <- length(tmp1)-4
			pars <- matrix(as.numeric(pars), items, dimensions+2)
			pars <- pars[,c(3:(dimensions+2),1,2)]
			colnames(pars) <- c(paste("a",1:dimensions,sep=""),"d","c")
		} else {
			dimensions <- length(tmp1)-3
			pars <- matrix(as.numeric(pars), items, dimensions+1)
			pars <- pars[,c(2:(dimensions+1),1)]
			colnames(pars) <- c(paste("a",1:dimensions,sep=""),"d")
		}
		
		if (as.irt.pars==TRUE) {
			n <- nrow(pars)
			pm <- as.poly.mod(n)
			pars <- as.irt.pars(pars, cat=rep(2,n), poly.mod=pm, dimensions=dimensions)
		}
	} else {
		tmp <- scan(file,what="character", quiet=TRUE)
		if (bifactor==FALSE) {
			flag <- NULL
			for (i in 1:length(tmp)) {
				if (nchar(tmp[i])==6) {
					if (substr(tmp[i],6,6)=="*") {
						flag <- c(flag,i)
					}
				}
			}
			tmp <- tmp[-flag]
		}
		# Identify the number of dimensions
		for (i in 6:20) {
			if (tmp[i]=="2") dimensions <- i-6
		}
		
		pars <- matrix(tmp, ncol=dimensions+5, byrow=T)
		pars <- pars[,-c(1:5)]
		pars <- matrix(as.numeric(pars), nrow(pars), dimensions)
		colnames(pars) <- paste("theta",1:dimensions,sep="")
	}
	
	return(pars)
}

read.ltm <- function(x, loc.out=FALSE, as.irt.pars=FALSE) {
	cls <- class(x)
	dimensions <- 1
	if (cls=="rasch") {
		pars <- cbind(coef(x)[,2:1])
		colnames(pars) <- c("a","b")
		cat <- rep(2,nrow(pars))
		cls <- "drm"
	} else if (cls=="tpm") {
		pars <- coef(x)[,3:1]
		colnames(pars) <- c("a","b","c")
		cat <- rep(2,nrow(pars))
		cls <- "drm"
	} else if (cls=="grm") {
		pars <- coef(x)
		if (is.matrix(pars)) {
			pars<- pars[,c(ncol(pars),1:(ncol(pars)-1))]
			cat <- rep(ncol(pars)-1,nrow(pars))
		} else {
			cat <- unlist(lapply(pars,length))
			names(cat) <- NULL
			p <- matrix(NA,length(pars),max(cat))
			for (i in 1:length(pars)) {
				p[i,1:cat[i]] <- pars[[i]][length(pars[[i]]):1]
			}
			pars <- p
		}
		if (x$IRT.param==FALSE) {
			pars[,2:ncol(pars)] <- pars[,2:ncol(pars)]/matrix(pars[,1],nrow(pars),ncol(pars)-1)
		}
		colnames(pars) <- c("a",paste("b",1:(ncol(pars)-1),sep=""))
	} else if (cls=="gpcm") {
		pars <- coef(x)
		if (is.matrix(pars)) {
			pars<- pars[,c(ncol(pars),1:(ncol(pars)-1))]
			cat <- rep(ncol(pars)-1,nrow(pars))
		} else {
			cat <- unlist(lapply(pars,length))
			names(cat) <- NULL
			p <- matrix(NA,length(pars),max(cat))
			for (i in 1:length(pars)) {
				p[i,1:cat[i]] <- pars[[i]][length(pars[[i]]):1]
			}
			pars <- p
		}
		if (x$IRT.param==FALSE) {
			pars[,2:ncol(pars)] <- -pars[,2:ncol(pars)]/matrix(pars[,1],nrow(pars),ncol(pars)-1)
		}
		colnames(pars) <- c("a",paste("b",1:(ncol(pars)-1),sep=""))
	} else if (cls=="ltm") {
		if (x$IRT.param==TRUE) {
			pars <- cbind(coef(x)[,2:1],0)
			colnames(pars) <- c("a","b","c")
		} else {
			if (x$ltst$factors==1) {
				pars <- cbind(coef(x)[,2:1])
				colnames(pars)[1:2] <- c("a","b")
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
	pars <- as.irt.pars(pars)
	
	if (as.irt.pars==FALSE) {
		pars <- pars@pars
	}
	
	return(pars)
}

read.icl <- function(file, poly.mod, ability=FALSE,  loc.out=FALSE, as.irt.pars=FALSE) {
	if (ability==FALSE) {
		pars <- read.table(file,sep="\t",fill=TRUE)
		pars <-cbind(pars, matrix(NA,nrow(pars),2*ncol(pars)))
		id <- 1
		repeat{
			if (pars[id,1]==id) {
				id <- id+1
			} else {
				tmp <- pars[id,!is.na(pars[id,])]
				tmp1 <- length(pars[id-1,!is.na(pars[id-1,])])
				pars[id-1,(tmp1+1):(tmp1+length(tmp))] <- tmp
				pars <- pars[-id,]
			}
			if (id==nrow(pars)) break
		}
		pars <- pars[,apply(is.na(pars),2,sum)<nrow(pars)]
		pars <- pars[,-1]
		
		# Determine the number of response categories
		cat <- apply(!is.na(pars),1,sum)
		cat[poly.mod@items$drm] <- 2
		
		pars <- sep.pars(pars, cat=cat, poly.mod=poly.mod, loc.out=loc.out)
		pars <- as.irt.pars(pars)
		
		if (as.irt.pars==FALSE) {
			pars <- pars@pars
		}
	} else {
		pars <- read.table(file, sep="\t")
	}
	return(pars)
}

read.bmirt <- function(file, ability=FALSE, loc.out=FALSE, pars.only=TRUE, as.irt.pars=FALSE) {
	if (ability==FALSE) {
		tmp <- scan(file, skip=1, quiet=TRUE)
		nc <- dimensions+20
		pars <- NULL
		start <- 1
		for (i in 1:(length(tmp)-2)) {
			if (tmp[i]==1 & tmp[i+1]==1 & tmp[i+2]==0) {
				len <- i-start
				pars <- rbind(pars,c(tmp[start:(i-1)],rep(NA,nc-len)))
				start <- i+3
			}
		}
		tmp <- apply(is.na(pars),2,sum)
		pars <- pars[,tmp!=nrow(pars)]
		cat <- pars[,2]
		cat[cat==1] <- 2
		pars <- pars[,-c(1,2)]
		
		# Determine the number of dimensions
		tmp <- pars[1,!is.na(pars[1,])]
		if (cat[1]==2) {
			dimensions <- length(tmp)-2
		} else {
			dimensions <- length(tmp)-cat[1]+1
		}
		
		n <- nrow(pars)
		
		# Create the poly.mod object
		if (length(cat[cat==2])==n) {
			# All dichotomous items
			pm <- as.poly.mod(n)
		} else {
			# Mixed-format items
			items <- 1:n
			pm <- as.poly.mod(n, c("drm","gpcm"), list(items[cat==2], items[cat>2]))
		}
		pars <- sep.pars(pars, cat=cat, poly.mod=pm, dimensions=dimensions, loc.out=loc.out)
		# Reformat the difficulty/step parameters to coincide with the traditional formulation of multidimensional models
		pars@b <- pars@b*-1
		pars <- as.irt.pars(pars)
		if (as.irt.pars==FALSE) {
			pars <- pars@pars
		}
	} else {
		pars <- read.table(file, sep=" ", skip=2)
		tmp <- apply(is.na(pars),2,sum)
		pars <- pars[,tmp!=nrow(pars)]
		pars <- pars[,-1]
		nms <- paste("theta",rep(1:dimensions,each=2),sep="")
		nms[seq(2,2*dimensions,2)] <- paste(nms[seq(2,2*dimensions,2)],"se",sep=".")
		colnames(pars) <- nms
		if (pars.only==TRUE) pars <- pars[,seq(1,2*dimensions,2)]
	}
	return(pars)
}