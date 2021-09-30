

GetNumberOfGOtermGenes <- function(go, species, supcounts, og2gene, folder, SizeThresh){
	fileinfo = file.info(paste(folder, "/", go, ".txt", sep=""))
	empty = rownames(fileinfo[fileinfo$size == 0, ])
	if(length(empty)==0){
		HGenelist <-  read.table(paste(folder, "/", go, ".txt", sep=""), h=F)[,1]
	}else{
		HGenelist <- c()
	}
	if(length(HGenelist) >= SizeThresh){
		GOList <- unique(og2gene$OG[which(og2gene$OG %in% rownames(supcounts) & og2gene$Gene %in% HGenelist)])
		return(length(unique(og2gene$Gene[which(og2gene$OG %in% GOList & og2gene$Species==species)])))
	}else{
		return(0)
	}
}

ScatterPercentagePlot <- function(vec1, vec2, col, lab1, lab2, xlim, ylim){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-3, ylim[2]+3), xlim=c(xlim[1]-3, xlim[2]+3), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.5)
	mtext(lab2, side = 2, line = 5, cex=1.5)
	points(vec1, vec2, col=col, bg=modif_alpha(col,.5), pch=21, cex=1.5)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
}


modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

modifColor <- function(col, change){
	apply(sapply(col, col2rgb)/255, 2,
	function(x)
	return(rgb(max(0,min(x[1]+change,1)), max(0,min(x[2]+change,1)), max(0,min(x[3]+change,1)))))
}


PlotColumnVertebrateType <- function(df, pos, vtypes, col, cextext, alpha=0, den=NULL){
	len <- length(df[,1])
	base <- 0
	w <- .8
	for(vt in c(1:length(vtypes))){
		value <- length(df[which(df$VertType==vtypes[vt]),1])/len*100
		if(value>0){
			if(!is.null(den)){
				polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[vt], .1), border=col[vt], lwd=4)		
			}
			polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[vt], alpha), border=col[vt], lwd=4, density=den)
			text(pos, base+value/2, labels=length(df[which(df$VertType==vtypes[vt]),1]), cex=cextext)
			base <- base+value
		}
	}
	par(xpd=TRUE) 
	text(pos, 100,  labels =length(df[,1]), pos=3, cex=cextext)
	par(xpd=FALSE) 
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
	return(list(lab1, lab2, log(Overlap/(group1*group2/Total)), Overlap, (group1*group2/Total), dp, ep))
}

PlotHypergeomTest_VertBlanTypes <- function(OG.df, VTypes, ATypes, vcols, thresh){
	# Testing significance coexistence of amphioxus and vertebrate categories
	HyperTests <- data.frame()
	for(vtype in VTypes){
		for(atype in ATypes){
			vnum <- length(OG.df[which(OG.df$VertType == vtype),1])
			anum <- length(OG.df[which(OG.df$BlanType == atype),1])
			onum <- length(OG.df[which(OG.df$VertType == vtype & OG.df$BlanType == atype),1])
			HyperTests <- rbind(HyperTests, unlist(HypergeometricTest(onum, anum, vnum, length(OG.df[,1]), atype, vtype, thresh)))
		}
	}
	colnames(HyperTests) <- c("BlanType", "VertType", "FoldChange", "Observed", "Expected", "HpvalDepleted", "HpvalEnriched")
	HyperTests$HBonfDepleted <- as.numeric(HyperTests$HpvalDepleted)*(length(HyperTests$HpvalDepleted)+length(HyperTests$HpvalEnriched))
	HyperTests$HBonfEnriched <- as.numeric(HyperTests$HpvalEnriched)*(length(HyperTests$HpvalDepleted)+length(HyperTests$HpvalEnriched))
	HyperTests$HResult <- rep("NA", length(HyperTests[,1]))
	HyperTests$HResult[which(HyperTests$HBonfDepleted <= PQvalThreshold)] <- rep("D", length(HyperTests$HResult[which(HyperTests$HBonfDepleted <= PQvalThreshold)]))
	HyperTests$HResult[which(HyperTests$HBonfEnriched <= PQvalThreshold)] <- rep("E", length(HyperTests$HResult[which(HyperTests$HBonfEnriched <= PQvalThreshold)]))
	print(HyperTests[,c(1,2,3,8,9,10)])

	colfunc <- colorRampPalette(c("indianred4", "indianred3", "indianred1", "white", "gold1", "gold3", "gold4"))
	gradientcolors <- colfunc(50)
	gradientlims <- c(-2,2)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,10), xlim=c(0,6), col=NA)
	mtext("B.lanceolatum", side = 2, line = 5, cex=1)
	mtext("Vertebrates", side = 1, line = 5, cex=1)
	for(v in c(1:length(VTypes))){
		for(a in c(1:length(ATypes))){
			color <- findColorInGRadient(HyperTests$FoldChange[which(HyperTests$BlanType==ATypes[a] & HyperTests$VertType==VTypes[v])], gradientlims, gradientcolors)
			if(color=="grey30"){dens <- 10}else{dens <- NULL}
			polygon(c(v-1,v,v,v-1), c(a-1,a-1,a,a), col=color, border="white", density=dens)
			if(color!="grey30" & HyperTests$HResult[which(HyperTests$BlanType==ATypes[a] & HyperTests$VertType==VTypes[v])]!="NA"){
					text(v-.5, a-.5, label="*", cex=1.5)
			}
		}
	}
	axis(2, at = seq(0.5,length(ATypes),1), labels=ATypes, lwd=NA, col = NA, las=1, cex.axis=1)

	text(x = seq(0.5,length(VTypes),1),
		y = par("usr")[3] - 0.45,
		labels = VTypes,
		xpd = NA,
		srt = 35,
		cex = 1,
		adj = 1)
	printgradientlegend(c(5.1,1), .03, .2, "log(FC)", gradientlims, gradientcolors)
}

BarPlotVertTypesInBlanGenes <- function(OG.df, ATypes, VTypes, vcols){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, length(ATypes)+.5), col=NA)
	mtext("% in each category", side = 2, line = 5, cex=1.5)
	PlotColumnVertebrateType(OG.df[which(OG.df$Blan==1),], 1, VTypes, vcols, 1.5)
	PlotColumnVertebrateType(OG.df[which(OG.df$Blan>=2),], 2, VTypes, vcols, 1.5, den=10)
	if(length(OG.df[which(OG.df$Blan==0),1])>0){
		PlotColumnVertebrateType(OG.df[which(OG.df$Blan==0),], 3, VTypes, vcols, 1.5, alpha=.2)
	}
	axis(1, at = c(1:length(ATypes)), labels=paste(ATypes, "\nin B.lan"), tick=FALSE, line=2, las=1, cex.axis=1.5)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
}



findColorInGRadient <- function(value, lims, colors){
	step <- (lims[2]-lims[1])/length(colors)
	bins <- seq(lims[1], lims[2]-step, step)
	value <- as.numeric(value)
	if(is.infinite(value) | value=="NaN"){
		return("grey30")
	}
	if(value < lims[1] | value > lims[2]){
		print(paste("Error: value out of limits in findColorInGRadient", value, lims[1], lims[2]))
		quit()
	}
	for(bin in c(1:length(bins))){
		if(value>=bins[bin] & value<bins[bin]+step){
			return(colors[bin])
		}
	}
}

printgradientlegend <- function(coord, height, width, lab, lims, colors){
	x <- coord[1]
	y <- coord[2]
	for(c in c(1:length(colors))){
		polygon(c(x, x+width, x+width, x),c(y+height*(c-1), y+height*(c-1), y+height*(c-1)+height, y+height*(c-1)+height), col=colors[c], border=NA)
	}
	text(x+width*1.5, y, label=lims[1], pos=4)
	text(x+width*1.5, y+height*length(colors), label=lims[2], pos=4)
	text(x+width*1.5, y+height*length(colors)/2, label=0, pos=4)
	text(x+width*1.5, y+height*length(colors)*1.1, label=lab, pos=3)
}
