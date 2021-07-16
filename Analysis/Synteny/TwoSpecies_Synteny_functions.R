



CalcMaxDist <- function(x, GeneData){
	tmp <- GeneData[which(GeneData$OG==x),]
	tmp <- tmp[order(tmp$Start),]
	maxdist <- 0
	if(length(tmp[,1])>1 & length(unique(tmp$Chr))==1){
		for(g1 in c(1:(length(tmp[,1])-1))){
			maxdist <- max(maxdist, tmp$Start[g1+1]-tmp$End[g1])
		}
	}
	return(maxdist)
}


CalcNumBetwGenes <- function(x, GeneData){
	tmp <- GeneData[which(GeneData$OG==x),]
	tmp <- tmp[order(tmp$Start),]
	numBG <- NA
	if(length(tmp[,1])>1 & length(unique(tmp$Chr))==1){
		numBG <- 0
		for(g1 in c(1:(length(tmp[,1])-1))){
			numBG <- numBG + length(GeneData[which(GeneData$End < tmp$Start[g1+1] & GeneData$Start > tmp$End[g1]),1])
		}
	}
	return(numBG)
}


RetrieveGenesInOtherSpecies <- function(x, GeneData, otherGeneData){
	xog <- unique(GeneData[which(GeneData[,"Gene"]==x), "OG"])
	otherx <- c()
	for(g in xog){
		otherx <- c(otherx, unique(otherGeneData[which(otherGeneData[,"OG"]==g),"Gene"]))
	} 
	if(length(otherx)==0){otherx <- NA}
	return(paste0(x, "_", otherx))
}


SyntenyAlongChrPlot <- function(column, main, GP.1, GP.2=NULL){
	if(column=="Index"){
		maxchrc <- "NumGenes"
		cummmaxchrc <- "CummNumGenes"
	}else if(column=="Midp"){
		maxchrc <- "Length"
		cummmaxchrc <- "CummLength"
	}else{
		quit("Invalid column name in SyntenyAlongChrPlot function")
	}
	axis1 <- ChrLen1[,cummmaxchrc]-ChrLen1[,maxchrc]/2
	axis2 <- ChrLen2[,cummmaxchrc]-ChrLen2[,maxchrc]/2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,ChrLen2[length(ChrLen2[,1]),cummmaxchrc]), xlim=c(0,ChrLen1[length(ChrLen1[,1]),cummmaxchrc]), col=NA)
	mtext(main, side = 3, line = 2, cex=2)
	mtext(SpeciesNames[which(ShortSpeciesNames==Species1)], side = 1, line = 6, cex=3)
	mtext(SpeciesNames[which(ShortSpeciesNames==Species2)], side = 2, line = 5, cex=3)
	abline(v=ChrLen1[,cummmaxchrc], lty=2, lwd=.5, col="grey60")
	abline(h=ChrLen2[,cummmaxchrc], lty=2, lwd=.5, col="grey60")
	Pos2.1 <- GP.2[,paste0(column,"1")]+ChrLen1[match(GP.2$Chr1, ChrLen1$Chr),cummmaxchrc]-ChrLen1[match(GP.2$Chr1, ChrLen1$Chr),maxchrc]
	Pos2.2 <- GP.2[,paste0(column,"2")]+ChrLen2[match(GP.2$Chr2, ChrLen2$Chr),cummmaxchrc]-ChrLen2[match(GP.2$Chr2, ChrLen2$Chr),maxchrc]
	points(Pos2.1, Pos2.2, col="grey60", pch=16, cex=.7)
	Pos1.1 <- GP.1[,paste0(column,"1")]+ChrLen1[match(GP.1$Chr1, ChrLen1$Chr),cummmaxchrc]-ChrLen1[match(GP.1$Chr1, ChrLen1$Chr),maxchrc]
	Pos1.2 <- GP.1[,paste0(column,"2")]+ChrLen2[match(GP.1$Chr2, ChrLen2$Chr),cummmaxchrc]-ChrLen2[match(GP.1$Chr2, ChrLen2$Chr),maxchrc]
	points(Pos1.1, Pos1.2, col="black", pch=16, cex=.7)
	print(head(Pos1.1))
	axis(1, at = c(0,ChrLen1[,cummmaxchrc]), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = c(0,ChrLen2[,cummmaxchrc]), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = axis1[seq(1, length(axis1), 2)], labels=ChrLen1$Chr[seq(1, length(axis1), 2)], lwd=NA, las=1, cex.axis=2.7, line=1.2)
	axis(2, at = axis2[seq(1, length(axis2), 2)], labels=ChrLen2$Chr[seq(1, length(axis2), 2)], lwd=NA, las=1, cex.axis=2.7)
	axis(1, at = axis1[seq(2, length(axis1), 2)], labels=ChrLen1$Chr[seq(2, length(axis1), 2)], lwd=NA, las=1, cex.axis=2.7, line=1.2)
	axis(2, at = axis2[seq(2, length(axis2), 2)], labels=ChrLen2$Chr[seq(2, length(axis2), 2)], lwd=NA, las=1, cex.axis=2.7)
	box()	
}



modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

AddHistlogYaxis <- function(h, col){
	for(i in c(1:length(h$counts))){
		polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,log(h$counts[i]),log(h$counts[i])), col=modif_alpha(col,.2), border=col)
	}
}



HypergeometricTest <- function(Overlap, group1, group2, Total, lab1, lab2, threshold){
	# Enrichment
	ep <- phyper(Overlap-1, group2, Total-group2, group1, lower.tail= FALSE)
	# Depletion
	dp <- phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE)
	if(ep <= threshold & dp <= threshold){
		cat("Error: both enriched and depleted!! \n")
		quit()
	}
	return(list(lab1, lab2, Overlap/(group1*group2/Total), Overlap, (group1*group2/Total), dp, ep))
}




