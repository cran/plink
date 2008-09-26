read.bilog <- function(file, ability=FALSE, pars.only=TRUE, as.irt.pars=FALSE) {

	if (ability==FALSE) {
		pars <- read.fwf(file, c(8,8,rep(10,13),4,1,1),sep="\t", skip=4)
		pars <- pars[,-15]
		colnames(pars) <- c("item", "subtest", "intercept","int.se", "slope", "slope.se", "threshold", "thresh.se", "dispersion",  "disp.se", "asymptote", "asymp.se", "drift", "drift.se","stream.loc","key","dummy.vals")
		if (pars.only==TRUE) pars <- pars[,c(5,7,11)]
	} else {
		pars <- read.fwf(file, c(3,13,4,5,10,12,12,11,10),sep="\t", skip=2)
		pars[,2] <- suppressWarnings(as.numeric(as.character(pars[,2])))
		pars[,1] <- c(0,pars[,1][-length(pars[,1])])
		pars[,2] <- c(0,pars[,2][-length(pars[,2])])
		pars <- pars[seq(2,nrow(pars),2),1:7]
		pars[,5] <- pars[,5]/100
		pars[pars[,7]==999,7] <- NA
		colnames(pars) <- c("group","id","n.pos","n.cor","p","logit","se")
	}
	
	if (as.irt.pars==TRUE) {
		if (ability==FALSE) {
			if (pars.only==FALSE) pars <- pars[,c(5,7,11)]
			n <- nrow(pars)
			pm <- as.poly.mod(n)
			pars <- as.irt.pars(pars, cat=rep(2,n), poly.mod=pm)
		}
	}
	return(pars)
}


read.parscale <- function(file, ability=FALSE, pars.only=TRUE, as.irt.pars=FALSE, loc.out=TRUE) {
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
		
		
		if (pars.only==TRUE) {
			pars <- pars[,seq(4,ncol(pars),2)]
		}
	} else {
		pars <- read.fwf(file, c(21,32,12,12))
		id <- as.character(pars[seq(1,nrow(pars),2),1])
		ability <- pars[seq(2,nrow(pars),2),3:4]
		ability[ability==999] <- NA
		pars <- cbind(id,ability)
		colnames(pars) <- c("id","ability","se")
	}
	
	if (as.irt.pars==TRUE) {
		if (ability==FALSE) {
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
	}
	return(pars)
}