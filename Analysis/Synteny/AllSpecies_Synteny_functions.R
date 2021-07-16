



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


SyntenyAlongChrPlot <- function(GP.df, column, main){
	if(column=="Index"){
		maxchrc <- "NumGenes"
		cummmaxchrc <- "CummNumGenes"
	}else if(column=="Midp"){
		maxchrc <- "Length"
		cummmaxchrc <- "CummLength"
	}else{
		quit("Invalid column name in SyntenyAlongChrPlot function")
	}
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,ChrLen2[length(ChrLen2[,1]),cummmaxchrc]), xlim=c(0,ChrLen1[length(ChrLen1[,1]),cummmaxchrc]), col=NA)
	mtext(main, side = 3, line = 2, cex=2)
	mtext(SpeciesNames[which(ShortSpeciesNames==Species1)], side = 1, line = 5, cex=1.5)
	mtext(SpeciesNames[which(ShortSpeciesNames==Species2)], side = 2, line = 5, cex=1.5)
	abline(v=ChrLen1[,cummmaxchrc], lty=2, lwd=.5, col="grey60")
	abline(h=ChrLen2[,cummmaxchrc], lty=2, lwd=.5, col="grey60")
	Pos1 <- GP.df[,paste0(column,"1")]+ChrLen1[match(GP.df$Chr1, ChrLen1$Chr),cummmaxchrc]-ChrLen1[match(GP.df$Chr1, ChrLen1$Chr),maxchrc]
	Pos2 <- GP.df[,paste0(column,"2")]+ChrLen2[match(GP.df$Chr2, ChrLen2$Chr),cummmaxchrc]-ChrLen2[match(GP.df$Chr2, ChrLen2$Chr),maxchrc]
	points(Pos1, Pos2, col="darkred", pch=16, cex=.5)
	axis(1, at = c(0,ChrLen1[,cummmaxchrc]), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = c(0,ChrLen2[,cummmaxchrc]), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(1, at = ChrLen1[,cummmaxchrc]-ChrLen1[,maxchrc]/2, labels=ChrLen1$Chr, lwd.ticks=NA, las=1, cex.axis=1.5)
	axis(2, at = ChrLen2[,cummmaxchrc]-ChrLen2[,maxchrc]/2, labels=ChrLen2$Chr, lwd.ticks=NA, las=1, cex.axis=1.5)
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

AddDens <- function(d, col){
	polygon(c(0,d$x), c(0,d$y), col=modif_alpha(col,.2), border=col)
}

plotHistAllSpecies <- function(vec, type, breaks, ylim, xlim, xlab, main){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(main, side = 3, line = 2, cex=1.5)
	mtext("Frequency", side = 2, line = 3, cex=1.5)
	mtext(xlab, side = 1, line = 3, cex=1.5)
	for(Species in ShortSpeciesNames){
		h <- hist(vec[which(names(vec)==paste0(type, Species))], breaks=breaks, plot=FALSE)
		AddHistlogYaxis(h, SpeciesColors[which(ShortSpeciesNames==Species)])
	}
	axis(1, at = seq(0, xlim[2], xlim[2]/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = c(0,log(c(1,10,100,1000,10000))), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plotDensAllSpecies <- function(vec, type, adj, ylim, xlim, xlab, main, sps, colors){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(main, side = 3, line = 2, cex=1.5)
	mtext("Density", side = 2, line = 3, cex=1.5)
	mtext(xlab, side = 1, line = 3, cex=1.5)
	for(s in sps){
		d <- density(vec[which(names(vec)==paste0(type, s))], adjust=adj)
		AddDens(d, colors[which(sps==s)])
	}
	axis(1, at = seq(0, xlim[2], xlim[2]/5), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(0, ylim[2], ylim[2]/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plotOGTypePerSpecies <- function(ogT){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,40), xlim=c(.5,length(ogT[,1])+.5), col=NA)
	mtext("% of total orthologous groups", side = 2, line = 5, cex=1)
	for(sp in c(1:length(ogT[,1]))){
		points(rep(sp,4), ogT[sp,c("SingleCopy", "Inter","Intra","Tandem")]/ogT$Total*100, pch=16, col=c("royalblue4", "orangered1", "orangered2", "orangered3"), cex=2)
	}
	axis(1, at = c(1:length(ogT[,1])), labels=ogT$Species, tick=FALSE, line=2, las=1, cex.axis=1)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1)
}

BoxPlot <- function(values, pos, col, cextext=1, w=.8, den=NULL, text=FALSE){
	s <- boxplot(values, plot=FALSE)
	points(rep(pos,length(s$out)), s$out, cex=.3, col="grey60")
	lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	if(!is.null(den)){
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=modifColor(col, .3), border=col, lwd=4)		
	}
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border=col, lwd=4, density=den)
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=NA, border=col, lwd=4)		
	#points(pos, m, pch=21, bg=col, col=col, cex=4)
	if(text){
		par(xpd=TRUE) 
		text(pos, 0,  labels =length(values), pos=1, cex=cextext)
		par(xpd=FALSE) 		
	}
}

plotVerticalDensAllSpecies<- function(vec, type, adj, ylim, ylab, main, sps, colors){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5, length(sps)+.5), col=NA)
	mtext(main, side = 3, line = 2, cex=1)
	mtext(ylab, side = 2, line = 3, cex=1)
	width <- .5
	for(s in c(1:length(sps))){
		d <- density(vec[which(names(vec)==paste0(type, sps[s]))], adjust=adj)
		polygon(c(s+d$y/max(d$y)*width/2, rev(s-d$y/max(d$y)*width/2)), c(d$x,rev(d$x)), col=modifColor(colors[which(sps==sps[s])], .2), border=colors[which(sps==sps[s])], lwd=2)
	}
	axis(1, at = c(1:length(sps)), labels=sps, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(0, ylim[2], ylim[2]/6), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

plotVerticalDensAllSpeciesSCD<- function(vec, adj, ylim, ylab, main, sps, colors){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(.5, length(sps)+.5), col=NA)
	mtext(main, side = 3, line = 2, cex=1)
	mtext(ylab, side = 2, line = 3, cex=1)
	width <- .5
	distSCD <- .2
	for(s in c(1:length(sps))){
		pos <- s-distSCD/2
		d <- density(vec[which(names(vec)==paste0("SC", sps[s]))], adjust=adj)
		polygon(c(pos+d$y/max(d$y)*width/2, rev(pos-d$y/max(d$y)*width/2)), c(d$x,rev(d$x)), col=modifColor(colors[which(sps==sps[s])], .2), border=colors[which(sps==sps[s])], lwd=2)
		pos <- s+distSCD/2
		d <- density(vec[which(names(vec)==paste0("D", sps[s]))], adjust=adj)
		polygon(c(pos+d$y/max(d$y)*width/2, rev(pos-d$y/max(d$y)*width/2)), c(d$x,rev(d$x)), col=modifColor(colors[which(sps==sps[s])], .05), border=colors[which(sps==sps[s])], lwd=2)
		polygon(c(pos+d$y/max(d$y)*width/2, rev(pos-d$y/max(d$y)*width/2)), c(d$x,rev(d$x)), col=modifColor(colors[which(sps==sps[s])], .2), border=colors[which(sps==sps[s])], lwd=2, density=20)
	}
	axis(1, at = c(1:length(sps)), labels=sps, lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(0, ylim[2], ylim[2]/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}

modifColor <- function(col, change){
	apply(sapply(col, col2rgb)/255, 2,
	function(x)
	return(rgb(max(0,min(x[1]+change,1)), max(0,min(x[2]+change,1)), max(0,min(x[3]+change,1)))))
}
