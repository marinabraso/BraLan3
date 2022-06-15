

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

GetNumberOfGOtermGenes <- function(go, species, supcounts, geneinf, folder){
	fileinfo = file.info(paste(folder, "/", go, ".txt", sep=""))
	empty = rownames(fileinfo[fileinfo$size == 0, ])
	if(length(empty)==0){
		HGenelist <-  read.table(paste(folder, "/", go, ".txt", sep=""), h=F)[,1]
	}else{
		HGenelist <- c()
	}
	OGList <- unique(geneinf$OG[which(geneinf$OG %in% supcounts$OG & geneinf$Gene %in% HGenelist)])
	return(length(unique(geneinf$Gene[which(geneinf$OG %in% OGList & geneinf$Species==species)])))
}

CalMeanExpGO <- function(go, ExpData, expcol, OGlist, og2gene, gofolder){
	fileinfo = file.info(paste(gofolder, "/", go, ".txt", sep=""))
	empty = rownames(fileinfo[fileinfo$size == 0, ])
	if(length(empty)==0){
		HGenelist <-  read.table(paste(gofolder, "/", go, ".txt", sep=""), h=F)[,1]
	}else{
		HGenelist <- c()
	}
	FinalOGList <- unique(og2gene$OG[which(og2gene$OG %in% OGlist & og2gene$Gene %in% HGenelist)])
	return(mean(ExpData[which(ExpData$Gene %in% og2gene$Gene[which(og2gene$OG %in% FinalOGList & og2gene$Species=="Blan")]),expcol]))	
}

ScatterPlotSmooth <- function(vec1, vec2, col, lab1, lab2, xlim, ylim){
	colfunc <- colorRampPalette(c("white", "gold", "darkorange", "firebrick"))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	smoothScatter(vec1, vec2, nbin = 128, bandwidth=1, colramp = colfunc, nrpoints = 0, add=TRUE, xaxs=FALSE, yaxs=FALSE, box=FALSE)
	points(vec1, vec2, col=col, pch=16, cex=.4, xpd = NA)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
}

ScatterPlotPointAlpha <- function(vec1, vec2, alpha, col, lab1, lab2, xlim, ylim){
	alpha <- alpha/max(alpha)*0.7
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	points(vec1, vec2, col=NA, bg=modif_alpha(col,alpha), pch=21, cex=2, xpd = NA)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
}

ScatterPlotPointSize <- function(vec1, vec2, sizes, col, lab1, lab2, xlim, ylim){
	sizes <- sizes/max(sizes)*5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	points(vec1, vec2, col=NA, bg=modif_alpha(col,.4), pch=21, cex=sizes, xpd = NA)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
}

ScatterPlotContour <- function(vec1, vec2, col, lab1, lab2, xlim, ylim){
	f <- kde2d(vec1, vec2, h=12, n = 50, lims = c(xlim, ylim))
	levels <- 16
	colfunc <- colorRampPalette(c("white", "gold", "darkorange", "firebrick"))
	colors <- colfunc(levels)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	fill.contour(f, nlevels = levels, col=colfunc(levels+2), axes=FALSE, frame.plot=FALSE)
	points(vec1, vec2, col=NA, bg=modif_alpha(col,.5), pch=21, cex=.5, xpd = NA)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
}

fill.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = function(n) hcl.colors(n, 
        "YlOrRd", rev = TRUE), col = color.palette(length(levels) - 
        1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
    xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, 
    ...) {
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    .filled.contour(x, y, z, levels, col)
}

ScatterPlot_GOExp <- function(vec1, vec2, sizes, col, lab1, lab2, lim){
	sizes <- sizes/max(sizes)*5
	df <- as.data.frame(cbind(vec1, vec2, sizes))
	df <- df[which(vec1<lim[2] & vec2<lim[2]),]
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(lim[1]-3, lim[2]+3), xlim=c(lim[1]-3, lim[2]+3), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	points(df$vec1, df$vec2, col=NA, bg=modif_alpha(col,.5), pch=21, cex=df$sizes)
	abline(0,1, col="darkred")
	axis(1, at = seq(lim[1],lim[2],(lim[2]-lim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(lim[1],lim[2],(lim[2]-lim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
}

ScatterPlot <- function(vec1, vec2, col, lab1, lab2, xlim, ylim){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-3, ylim[2]+3), xlim=c(xlim[1]-3, xlim[2]+3), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	points(vec1, vec2, col=NA, bg=modif_alpha(col,.5), pch=21, cex=1)
	abline(0,1, col="darkred")
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4), lwd.ticks=1, las=1, cex.axis=1.5)
}

VerticalDensitiesPlot <- function(vec1, vec2, values1, col, lab1, lab2, xlim, ylim){
	colfunc <- colorRampPalette(c("gold", "seagreen"))
	colors <- colfunc(length(values1))
	colors.alpha <- apply(sapply(colors, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=0.2))  
	colors.alpha2 <- apply(sapply(colors, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=0.1))  
	width <- .8

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.5)
	mtext(lab2, side = 2, line = 5, cex=1.5)
	for(r in c(1:length(values1))){
		interv <- as.numeric(unlist(strsplit(as.character(values1[r]), "-")))
		if(length(interv)>1){
			vec <- vec2[which(vec1>=interv[1] & vec1<=interv[2])]
		}else{
			vec <- vec2[which(vec1==values1[r])]
		}
		if(length(vec)>1){
			d  <- density(vec, adjust = 1)
			polygon(c(r+d$y/max(d$y)*width/2, rev(r-d$y/max(d$y)*width/2)), c(d$x,rev(d$x)), col=colors.alpha[r], border=colors[r], lwd=1)
		}
		points(jitter(rep(r, length(vec)), amount=width/2), vec, pch=16, col=colors.alpha2[r], cex=.5)
	}
	axis(1, at = c(1:length(values1)), labels=values1, lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
}

Scatter2HeatmapPlot <- function(vecx, vecy, col, labx, laby, main, xlim, ylim){
	xstep <- 1
	ystep <- .5
	colfunc <- colorRampPalette(c(col[2],col[3]))
	colors <- c(col[1], colfunc(50))
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-ystep/2,ylim[2]+ystep/2), xlim=c(xlim[1]-xstep/2,xlim[2]+xstep/2), col=NA)
	mtext(labx, side = 1, line = 4, cex=1.2)
	mtext(laby, side = 2, line = 3, cex=1.2)
	mtext(main, side = 3, line = 2, cex=1.2)
	df <- as.data.frame(cbind(vecx, vecy))
	props <- c()
	for(x in seq(xlim[1]-xstep/2, xlim[2]+xstep/2,xstep)){
		for(y in seq(ylim[1]-ystep/2, ylim[2]+ystep/2,ystep)){
			proppoints <- length(df[which(df$vecx>=x & df$vecx<x+xstep & df$vecy>=y & df$vecy<y+ystep),1])/length(df[,1])
			props <- c(props, proppoints)
			cellcolor <- findColorInGRadient(proppoints, c(0,max(props)), colors)
			polygon(c(x,x+xstep,x+xstep,x), c(y,y,y+ystep,y+ystep), col=cellcolor, border=NA)
		}
	}
	#lines(c(1:length(values1)), medians, col="darkred")
	axis(1, at = seq(xlim[1],xlim[2],xstep), lwd.ticks=1, las=1, cex.axis=1.2)
	axis(2, at = seq(ylim[1],ylim[2],1), lwd.ticks=1, las=1, cex.axis=1.2)
	box()
}

findColorInGRadient <- function(value, lims, colors){
	step <- (lims[2]-lims[1])/length(colors)
	bins <- seq(lims[1], lims[2]-step, step)
	value <- as.numeric(value)
	if(is.infinite(value)){
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
	if(value==lims[2]){
		return(colors[length(colors)])
	}
}

ContingencyTableCN <- function(vec1, vec2, outfile){
	spcor <- cor.test(vec1, vec2, method="spearman")
	write.table("Spearman correlation: ", file = outfile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	write.table(paste("rho = ", spcor$estimate), file = outfile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE, append=TRUE)
	write.table(paste("corr p.value = ", spcor$p.value), file = outfile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE, append=TRUE)
	breaks <- c(-Inf, 2.5, Inf)
	ContingencyTable(vec1, vec2, breaks, outfile)
	breaks <- c(seq(1.5, 5.5, 1), Inf)
	ContingencyTable(vec1, vec2, breaks, outfile)
}

ContingencyTable <- function(vec1, vec2, breaks, outfile=NA){
	categ1 <- cut(vec1, breaks = breaks)
	categ2 <- cut(vec2, breaks = breaks)
	table <- table(categ1, categ2)
	chisq <- chisq.test(table)
	print(chisq$observed/chisq$expected)
	print(chisq$p.value)
	write.table("Contingency table: ", file = outfile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE, append=TRUE)
	write.table(chisq$observed/chisq$expected, file = outfile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE, append=TRUE)
	write.table(paste("chisq p.value = ", chisq$p.value), file = outfile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE, append=TRUE)
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

PlotHypergeomTest_VertBlanTypes <- function(OG.df, VTypes, ATypes, vcols, thresh, rfolder){
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
	write.table(HyperTests[,c(1,2,3,8,9,10)], file = paste0(rfolder, "/HypergeometricTestTable_", paste0(substr(VTypes, 1, 1), collapse=""), "_", paste0(substr(ATypes, 1, 1), collapse=""),".txt"), quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

	colfunc <- colorRampPalette(c("indianred4", "indianred3", "indianred1", "white", "gold1", "gold3", "gold4"))
	gradientcolors <- colfunc(50)
	gradientlims <- c(-2.5,2.5)
	colmatrix <- matrix(rep(NA, length(ATypes)*length(VTypes)),nrow=length(VTypes),ncol=length(ATypes),byrow=T)
	FCmatrix <- matrix(rep(NA, length(ATypes)*length(VTypes)),nrow=length(VTypes),ncol=length(ATypes),byrow=T)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,10), xlim=c(0,6), col=NA)
	mtext("B.lanceolatum", side = 2, line = 5, cex=1)
	mtext("Vertebrates", side = 1, line = 5, cex=1)
	for(v in c(1:length(VTypes))){
		for(a in c(1:length(ATypes))){
			color <- findColorInGRadient(log2(as.numeric(HyperTests$FoldChange[which(HyperTests$BlanType==ATypes[a] & HyperTests$VertType==VTypes[v])])), gradientlims, gradientcolors)
			FCmatrix[v,a] <- as.numeric(HyperTests$FoldChange[which(HyperTests$BlanType==ATypes[a] & HyperTests$VertType==VTypes[v])])
			colmatrix[v,a] <- color
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
	return(list(colmatrix=colmatrix, FCmatrix=FCmatrix))
}

BarPlotVertTypesInBlanGenes <- function(OG.df, ATypes, VTypes, vcols){
	barDensities <- c(NA, 10, NA)
	barAlphas <- c(0, 0, .2)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, length(ATypes)+.5), col=NA)
	mtext("% in each category", side = 2, line = 4, cex=1.2)
	mtext("B. lanceolatum", side = 1, line = 4, cex=1.2)
	for(atype in c(1:length(ATypes))){
		PlotColumnOtherBranchType_percent(OG.df[which(OG.df$BlanType==ATypes[atype]),], atype, "VertType", VTypes, vcols, 1.5, den=barDensities[atype], alpha=barAlphas[atype])
	}
	axis(1, at = c(1:length(ATypes)), labels=ATypes, tick=FALSE, line=1, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)
}

BarPlotVertTypesBlanTypes <- function(df, ATypes, VTypes, vcols, acols){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, 6.5), col=NA)
	mtext("% in each category", side = 2, line = 4, cex=1.2)
	mtext("B. lanceolatum", side = 1, line = 4, cex=1.2)
	PlotColumnOtherBranchType(df[which(df$BlanType==ATypes[1]),], length(df[,1]), 1, "VertType", VTypes, vcols, 1)
	PlotColumnOtherBranchType(df[which(df$BlanType==ATypes[2]),], length(df[,1]), 2, "VertType", VTypes, vcols, 1)
	PlotColumnOtherBranchType(df[which(df$VertType==VTypes[1]),], length(df[,1]), 4, "BlanType", ATypes, acols, 1)
	PlotColumnOtherBranchType(df[which(df$VertType==VTypes[2]),], length(df[,1]), 5, "BlanType", ATypes, acols, 1)
	PlotColumnOtherBranchType(df[which(df$VertType==VTypes[3]),], length(df[,1]), 6, "BlanType", ATypes, acols, 1)
	axis(1, at = c(1,2,4,5,6), labels=c(ATypes, VTypes), tick=FALSE, line=1, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)
}

BarPlotSpeciesNumGenesOG <- function(df, vspecies){
	w <- .8
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, length(vspecies)+1+.5), col=NA)
	mtext("% of duplicated genes in each species", side = 2, line = 4, cex=1.2)
	pos <- 1
	value <- df$dupGenes[which(df$Species=="Blan")]/df$totalGenes[which(df$Species=="Blan")]*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(0,0,value,value), col="grey60", border=NA)
	text(pos, value/2, labels=paste(format(round(value, 1), nsmall = 1),"%"), cex=1)
	for(v in c(1:length(vspecies))){ 
		pos <- v+1
		value <- df$dupGenes[which(df$Species==vspecies[v])]/df$totalGenes[which(df$Species==vspecies[v])]*100
		polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(0,0,value,value), col="grey60", border=NA)
		text(pos, value/2, labels=paste(format(round(value, 1), nsmall = 1),"%"), cex=1)
	}
	axis(1, at = c(1:(length(vspecies)+1)), labels=c("Blan", vspecies), tick=FALSE, line=1, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, length(vspecies)+1+.5), col=NA)
	mtext("% of duplicated orthogroups in each species", side = 2, line = 4, cex=1.2)
	pos <- 1
	value <- df$dupOG[which(df$Species=="Blan")]/df$totalOG[which(df$Species=="Blan")]*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(0,0,value,value), col="grey60", border=NA)
	text(pos, value/2, labels=paste(format(round(value, 1), nsmall = 1),"%"), cex=1)
	for(v in c(1:length(vspecies))){ 
		pos <- v+1
		value <- df$dupOG[which(df$Species==vspecies[v])]/df$totalOG[which(df$Species==vspecies[v])]*100
		polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(0,0,value,value), col="grey60", border=NA)
		text(pos, value/2, labels=paste(format(round(value, 1), nsmall = 1),"%"), cex=1)
	}
	axis(1, at = c(1:(length(vspecies)+1)), labels=c("Blan", vspecies), tick=FALSE, line=1, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)
}

BarPlotVertTypesInBlanGenes_wexpected <- function(OG.df, ATypes, VTypes, vcols, colmatrix){
	tBlanType <- table(OG.df$BlanType)
	tBlanType <- tBlanType[match(ATypes, names(tBlanType))]
	tVertType <- table(OG.df$VertType)
	tVertType <- tVertType[match(VTypes, names(tVertType))]
	obsBlanVertTypes <- t(table(OG.df[,c("BlanType", "VertType")]))
	obsBlanVertTypes <- obsBlanVertTypes[match(VTypes, rownames(obsBlanVertTypes)),match(ATypes, colnames(obsBlanVertTypes))]
	print(colmatrix)

	expBlanVertTypes <- matrix(rep(NA, length(tBlanType)*length(tVertType)),nrow=length(tVertType),ncol=length(tBlanType),byrow=T)
	for(i in c(1:length(tBlanType))){
		for(j in c(1:length(tVertType))){
			expBlanVertTypes[j,i] <- tBlanType[i]*tVertType[j]/sum(tBlanType)
		}
	}
	tBlanType <- tBlanType/sum(tBlanType)*100
	tVertType <- tVertType/sum(tVertType)*100
	obsBlanVertTypes <- c(obsBlanVertTypes)/sum(obsBlanVertTypes)*100
	expBlanVertTypes <- c(expBlanVertTypes)/sum(expBlanVertTypes)*100

	den <- 45
	alpha <- 0.1
	width <- .8
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, 4), col=NA)
	mtext("% in each category", side = 2, line = 4, cex=1.2)
	mtext("B. lanceolatum", side = 1, line = 4, cex=1.2)
	PlotColumnFromArray(tBlanType, 1, rep("white",length(tBlanType)), rep("white",length(tBlanType)), c(NULL,den), alpha, width)
	PlotBetweenColumnFrom2Arrays(expBlanVertTypes, obsBlanVertTypes, 2, 3.5, c(colmatrix), width)

	PlotColumnFromArray(expBlanVertTypes, 2, c(vcols, vcols), rep("white", length(vcols)*2), rep(den, length(tVertType)*2), alpha, width)
	PlotColumnFromArray(obsBlanVertTypes, 3.5, c(vcols, vcols), c(vcols, vcols), rep(NULL, length(tVertType)*2), alpha, width)

	abline(h=tBlanType[1]/sum(tBlanType)*100, lty=2, lwd=2, col="black")
	axis(1, at = c(1,2,3.5), labels=c("", "Expected", "Observed"), tick=FALSE, line=1, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)
}

BarPlotVertTypesInBlanGenes_wexpected2 <- function(OG.df, ATypes, VTypes, vcols, FCmatrix){
	tBlanType <- table(OG.df$BlanType)
	tBlanType <- tBlanType[match(ATypes, names(tBlanType))]
	tVertType <- table(OG.df$VertType)
	tVertType <- tVertType[match(VTypes, names(tVertType))]
	obsBlanVertTypes <- t(table(OG.df[,c("BlanType", "VertType")]))
	obsBlanVertTypes <- obsBlanVertTypes[match(VTypes, rownames(obsBlanVertTypes)),match(ATypes, colnames(obsBlanVertTypes))]
	print(FCmatrix)

	expBlanVertTypes <- matrix(rep(NA, length(tBlanType)*length(tVertType)),nrow=length(tVertType),ncol=length(tBlanType),byrow=T)
	for(i in c(1:length(tBlanType))){
		for(j in c(1:length(tVertType))){
			expBlanVertTypes[j,i] <- tBlanType[i]*tVertType[j]/sum(tBlanType)
		}
	}
	tBlanType <- tBlanType/sum(tBlanType)*100
	tVertType <- tVertType/sum(tVertType)*100
	obsBlanVertTypes <- c(obsBlanVertTypes)/sum(obsBlanVertTypes)*100
	expBlanVertTypes <- c(expBlanVertTypes)/sum(expBlanVertTypes)*100

	den <- 45
	alpha <- 0.1
	width <- .8
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, 3), col=NA)
	mtext("% of B. lanceolatum orthogroups", side = 2, line = 4, cex=1.2)
	PlotColumnFromArray(tBlanType, 1, vcols[which(VTypes %in% ATypes)], vcols[which(VTypes %in% ATypes)], NULL, alpha, width)
	PlotColumnFromArray(tVertType, 2, vcols, vcols, NULL, alpha, width)
	axis(1, at = c(1,2), labels=c("In B. lanceolatum", "In vertebrates"), tick=FALSE, line=1, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)
	
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, 3.5), col=NA)
	mtext("% of B. lanceolatum orthogroups", side = 2, line = 4, cex=1.2)
	PlotBetweenColumnFrom2Arrays2(expBlanVertTypes, obsBlanVertTypes, 1, 2.5, c(FCmatrix), width)
	#PlotColumnFromArray(expBlanVertTypes, 1, c(vcols, vcols), rep("white", length(vcols)*2), rep(den, length(tVertType)*2), alpha, width)
	PlotColumnFromArray(expBlanVertTypes, 1, c(vcols, vcols), c(vcols, vcols), rep(NULL, length(tVertType)*2), alpha, width)
	PlotColumnFromArray(obsBlanVertTypes, 2.5, c(vcols, vcols), c(vcols, vcols), rep(NULL, length(tVertType)*2), alpha, width)

	abline(h=tBlanType[1]/sum(tBlanType)*100, lty=2, lwd=2, col="black")
	axis(1, at = c(1,2.5), labels=c("Expected", "Observed"), tick=FALSE, line=1, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)
}

PlotBetweenColumnFrom2Arrays <- function(vec1, vec2, pos1, pos2, cols, w){
	base1 <- 0
	base2 <- 0
	for(v in c(1:length(vec1))){
		if(vec1[v]>0 | vec2[v]){
			polygon(c(pos1+w/2,pos2-w/2,pos2-w/2,pos1+w/2), c(base1,base2,base2+vec2[v],base1+vec1[v]), col=cols[v], border=NA, lwd=2)
			base1 <- base1+vec1[v]
			base2 <- base2+vec2[v]
		}
	}
}

PlotBetweenColumnFrom2Arrays2 <- function(vec1, vec2, pos1, pos2, fc, w){
	w <- w*1.2
	base1 <- 0
	base2 <- 0
	for(v in c(1:length(vec1))){
		if(vec1[v]>0 | vec2[v]){
			if(vec2[v]>vec1[v]){
				arrows(pos1+w/2, base1+vec1[v]/2, pos2-w/2, base2+vec2[v]/2, length = 0.1, angle = 20, code = 2, col="black", lwd=3)
				text((pos1+pos2)/2, (base1+vec1[v]/2+base2+vec2[v]/2)/2, labels=format(round(fc[v], 1), nsmall = 1), cex=1)
			}
			base1 <- base1+vec1[v]
			base2 <- base2+vec2[v]
		}
	}
}

PlotColumnFromArray <- function(vec, pos, cols, bgcols, den=c(rep(NULL,length(vec))), alpha=0.1, w){
	base <- 0
	for(v in c(1:length(vec))){
		if(vec[v]>0){
			if(!is.null(den[v])){
				polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+vec[v],base+vec[v]), col=bgcols[v], border=NA, lwd=0.5)	
			}
			polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+vec[v],base+vec[v]), col=cols[v], border=NA, lwd=0.5, density=den[v])
			text(pos, base+vec[v]/2, labels=paste(format(round(vec[v], 1), nsmall = 1),"%"), cex=1)
			base <- base+vec[v]
		}
	}
	par(xpd=TRUE) 
	#text(pos, 100, labels=sum(vec), pos=3, cex=1.5)
	par(xpd=FALSE) 
}

PlotColumnOtherBranchType <- function(df, total, pos, ColName, types, col, cextext, alpha=0, den=NULL){
	base <- 0
	w <- .8
	for(t in c(1:length(types))){
		value <- length(df[which(df[,ColName]==types[t]),1])/total*100
		if(value>0){
			if(!is.null(den)){
				polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[t], .1), border=col[t], lwd=2)		
			}
			polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[t], alpha), border=col[t], lwd=2, density=den)
			text(pos, base+value/2, labels=paste(format(round(value, 1), nsmall = 1),"%"), cex=cextext)
			base <- base+value
		}
	}
	par(xpd=TRUE) 
	text(pos, 100,  labels =length(df[,1]), pos=3, cex=cextext)
	par(xpd=FALSE) 
}

PlotColumnOtherBranchType_percent <- function(df, pos, ColName, types, col, cextext, alpha=0, den=NULL){
	len <- length(df[,1])
	base <- 0
	w <- .8
	for(t in c(1:length(types))){
		value <- length(df[which(df[,ColName]==types[t]),1])/len*100
		if(value>0){
			if(!is.null(den)){
				polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[t], .1), border=col[t], lwd=2)		
			}
			polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[t], alpha), border=col[t], lwd=2, density=den)
			text(pos, base+value/2, labels=paste(format(round(value, 1), nsmall = 1),"%"), cex=cextext)
			base <- base+value
		}
	}
	par(xpd=TRUE) 
	text(pos, 100,  labels =length(df[,1]), pos=3, cex=cextext)
	par(xpd=FALSE) 
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

CalcMaxDist <- function(x, df){
	tmp <- df[which(df$OG==x),]
	tmp <- tmp[order(tmp$Start),]
	maxdist <- 0
	if(length(tmp[,1])>1 & length(unique(tmp$Chr))==1){
		for(g1 in c(1:(length(tmp[,1])-1))){
			maxdist <- max(maxdist, tmp$Start[g1+1]-tmp$End[g1])
		}
	}
	return(maxdist)
}

CalcIfConsecutive <- function(x, df){
	tmp <- df[which(df$OG==x),]
	tmp <- tmp[order(tmp$Start),]
	consec <- FALSE
	if(length(tmp[,1])>1){
		if(length(unique(tmp$Chr))==1 & max(diff(tmp$Position))==1){
			consec <- TRUE
		}
	}
	return(consec)
}

TandemIntraInterPerSpecies <- function(df, spec, stype){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, length(spec)+8), col=NA)
	mtext("% in each category", side = 2, line = 4, cex=1.2)
	w <- .6
	col <- c("gold", "forestgreen","darkgreen")
	for(sp in c(1:length(spec))){
		supdf <- df[which(df[,spec[sp]]>1 & df[,stype[sp]]=="Small-scale\nduplicates"),]
		total <- length(supdf[,1])
		inter <- length(supdf$OG[which(supdf[,paste0("NumChrs",spec[sp])]>1)])/total*100
		intra <- length(supdf$OG[which(supdf[,paste0("NumChrs",spec[sp])]==1 & supdf[,paste0("Consecutive",spec[sp])]==FALSE)])/total*100
		tandem <- length(supdf$OG[which(supdf[,paste0("NumChrs",spec[sp])]==1 & supdf[,paste0("Consecutive",spec[sp])]==TRUE)])/total*100
		polygon(c(sp-w/2,sp+w/2,sp+w/2,sp-w/2), c(0,0,inter,inter), col=col[1], border=NA, lwd=4)
		polygon(c(sp-w/2,sp+w/2,sp+w/2,sp-w/2), c(inter,inter,inter+intra,inter+intra), col=col[2], border=NA, lwd=4)
		polygon(c(sp-w/2,sp+w/2,sp+w/2,sp-w/2), c(inter+intra,inter+intra,inter+intra+tandem,inter+intra+tandem), col=col[3], border=NA, lwd=4)
	}
	axis(1, at = c(1:length(spec)), labels=spec, tick=FALSE, line=0, las=1, cex.axis=1.2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.2)
	legend("topright", c("Multichromosomal","Distant monochromosomal","Tandem"), pch=15, col=col , bty = "n", pt.cex=1.5, cex=1, xjust = 0, yjust = 0)
}

prepare_DrerGeneData <- function(file, bgeefile, tissues, rfolder){
	pfolder <- system("pwd")
	if(file.exists(file)){
		Data <- read.delim(file, h=T, stringsAsFactors=F)
	}else{
		library(BgeeDB)
		setwd(paste0("./", rfolder))
		DrerBgee <- Bgee$new(species="Danio_rerio", dataType="rna_seq")
		DrerBgeeData <- getData(DrerBgee)
		write.table(DrerBgeeData, file = bgeefile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
		Data <- as.data.frame(cbind(unique(DrerBgeeData$Gene.ID)))
		colnames(Data) <- c("Gene")
		for(i in c(1:length(tissues[,1]))){
			system_out <- system(paste0("cat ", bgeefile, " | sed 's/\"//g' | awk -F'\\t' -v tissue=\"", tissues$Drer[i], "\" '{if($7==tissue){print $5\"\\t\"$13\"\\t\"$16}}' | sort -k1,1 | awk '{if(g!=$1){if(NR!=1){mTPM=mTPM/num;} print g\"\\t\"mTPM\"\\t\"pres; g=$1; num=1;mTPM=$2;pres=$3}else{num=num+1;mTPM=mTPM+$2;if($3==\"present\"){pres=$3}}}END{mTPM=mTPM/num; print g\"\\t\"mTPM\"\\t\"pres;}' | tail -n +2"), intern=T)
			tData <- read.table(text=system_out, h=F, sep = "\t")
			Data[,paste0(tissues$Name[i],"TPM")] <- tData[match(Data$Gene, tData[,1]),2]
			Data[,paste0(tissues$Name[i],"Presence")] <- tData[match(Data$Gene, tData[,1]),3]
		}
		write.table(Data, file = file, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
		setwd(pfolder)
	}
	return(Data)
}

prepare_BlanGeneData <- function(file, MetaInfoFile, tpmfile){
	if(file.exists(file)){
		Data <- read.delim(file, h=T, stringsAsFactors=F)
	}else{
		system_out <- system(paste("cat ", MetaInfoFile," | cut -f1,13,24", sep=""), intern=T)
		MetaRNA <- read.table(text=system_out, h=F, sep = "\t")
		colnames(MetaRNA) <- c("Sample","Age","Tissue")
		MetaRNA$Tissue <- sub(" ", ".", MetaRNA$Tissue)
		MetaRNA$Age <- sub(" ", ".", MetaRNA$Age)
		MetaRNA$Age <- sub("-", ".", MetaRNA$Age)

		TPM <- read.table(tpmfile, h=T)
		TPM <- TPM[,MetaRNA$Sample[which(MetaRNA$Tissue!="egg" & MetaRNA$Sample%in%colnames(TPM))]]
		MetaRNA <- MetaRNA[which(MetaRNA$Sample %in% colnames(TPM)),]

		Tissues <- unique(MetaRNA$Tissue[which(MetaRNA$Age=="adult")])
		Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads")
		EmbAges <- unique(MetaRNA$Age[which(MetaRNA$Age!="adult")])
		EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre.metamorphic.larvae")

		# TAU calculation per adult tissue
		GeneData <- data.frame("Gene" = rownames(TPM))
		for(t in c(1:length(Tissues))){
			samples <- MetaRNA$Sample[which(MetaRNA$Tissue==Tissues[t])]
			if(length(samples) == 1){
				GeneData <- cbind(GeneData, as.numeric(TPM[, samples]))
			}else{
				GeneData <- cbind(GeneData, as.numeric(rowMeans(TPM[, samples])))	
			}
		}
		for(a in c(1:length(EmbAges))){
			samples <- MetaRNA$Sample[which(MetaRNA$Age==EmbAges[a])]
			if(length(samples) == 1){
				GeneData <- cbind(GeneData, as.numeric(TPM[, samples]))
			}else{
				GeneData <- cbind(GeneData, as.numeric(rowMeans(TPM[, samples])))	
			}
		}
		colnames(GeneData) <- c("Gene", Tissues, EmbAges)
		m.GeneData <- as.matrix(GeneData[,Tissues])
		Tau <- apply(log2(m.GeneData+1), 1, compute.tau)
		MaxGroup <- apply(log2(m.GeneData+1), 1, max.col)
		GeneData$TauTissues <- Tau
		GeneData$MaxTissue <- MaxGroup
		m.GeneData <- as.matrix(GeneData[,EmbAges])
		Tau <- apply(log2(m.GeneData+1), 1, compute.tau)
		MaxGroup <- apply(log2(m.GeneData+1), 1, max.col)
		GeneData$TauEmbAge <- Tau
		GeneData$MaxEmbAge <- MaxGroup

		GeneData$MeanAdult <- rowMeans(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age=="adult")]])
		GeneData$MaxAdult <- apply(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age=="adult")]],1,max)
		GeneData$MeanEmbr <- rowMeans(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age!="adult")]])
		GeneData$MaxEmbr <- apply(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age!="adult")]],1,max)
		return(GeneData)
		write.table(Data, file = file, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	}
}

compute.tau <- function(exp){
	if(max(exp)==0){
	  return(NA)
	}
	n=length(exp)
	newexp=exp/max(exp)
	tau=sum(1-newexp)/(n-1)
	return(tau)
}

max.col <- function(exp){
	if(max(exp)==0){
	  return(NA)
	}
	maxcollist <- paste(names(exp[which(exp==max(exp))]), sep=":")
	return(maxcollist)
}

BoxPlot_BlanTypes_VertTypes <- function(Genes, Ginf, column, ylab, ylim, vtypes, vcolors){
	dist <- .15

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0.5, length(vtypes)+.5), col=NA)
	mtext(ylab, side = 2, line = 5, cex=1.2)
	for(t in c(1:length(vtypes))){
		BoxPlot(Genes[which(Genes$Gene %in% Ginf$Gene[which(Ginf$BlanType=="Single-copy" & Ginf$VertType==vtypes[t])]), column], t-dist, vcolors[t], .7, .4)
		BoxPlot(Genes[which(Genes$Gene %in% Ginf$Gene[which(Ginf$BlanType=="Small-scale\nduplicates" & Ginf$VertType==vtypes[t])]), column], t+dist, vcolors[t], .7, .4, 10)
	}
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	text(x = c(1:4), y = par("usr")[3] - (ylim[2]-ylim[1])/10, labels = vtypes, xpd = NA, srt = 40, cex = 1.3, adj = .9)
}

BoxPlot <- function(values, pos, col, cextext=1, w=.8, den=NULL, text=FALSE){
	s <- boxplot(values, plot=FALSE)
	lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	if(!is.null(den)){
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=modifColor(col, .3), border=col, lwd=3)		
	}
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border=col, lwd=3, density=den)
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=NA, border=col, lwd=3)		
	if(text){
		par(xpd=TRUE) 
		text(pos, 0,  labels =length(values), pos=1, cex=cextext)
		par(xpd=FALSE) 		
	}
}

Hist_ExpressDomanis <- function(gp, og, t1, t2, tissues, lab1, lab2, main){
	breaks <- seq(-length(tissues[,1])-.5,length(tissues[,1])+.5,1)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(-length(tissues[,1]),length(tissues[,1])), col=NA)
	mtext(main, side = 3, line = 1, cex=1.2)
	mtext(lab1, side = 1, line = 4, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	add_relative_hist(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], breaks, "cornflowerblue")
	if(!is.na(t1)){
		add_relative_hist(gp$DiffDom[which(gp$BlanType==t1 & gp$VertType==t2)], breaks, "tomato3")
		add_relative_hist(og$DiffDom[which(og$BlanType==t1 & og$VertType==t2)], breaks, "tomato3", 20)
	}
	axis(1, at = seq(-length(tissues[,1])+1,length(tissues[,1]),2), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.5)
}

add_relative_hist <- function(vec, b, col, dens=NULL){
	h <- hist(vec, breaks=b, plot=FALSE)
	h$counts <- h$counts/length(vec)
	for(i in c(1:length(h$counts))){
		polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,h$counts[i],h$counts[i]), col=modif_alpha(col,.3), border=col, density=dens, lwd=2)
	}
	return(h$counts)
}

Dist_ExpressDomanis <- function(vec, col, tissues, lab1, lab2, main, adj){
	breaks <- seq(-length(tissues[,1])-.5,length(tissues[,1])+.5,1)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(-length(tissues[,1]),length(tissues[,1])), col=NA)
	mtext(main, side = 3, line = 1, cex=1.2)
	d  <- density(vec, adjust = adj)
	lines(d$x, d$y, col=col, lwd=3)
	axis(1, at = seq(-length(tissues[,1])+1,length(tissues[,1]),2), lwd.ticks=1, las=1, cex.axis=1.5)
}

DifferenceWithSC <- function(gp, og, tissues, lab1, lab2){
	breaks <- seq(-length(tissues[,1])-.5,length(tissues[,1])+.5,1)
	width <- .7
	color <- "tomato3"
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0-.5, 8+.5), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1)
	mtext(lab2, side = 2, line = 5, cex=1)
	# D - SC
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], gp$DiffDom[which(gp$BlanType=="Small-scale\nduplicates" & gp$VertType=="Single-copy")], 1, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], og$DiffDom[which(og$BlanType=="Small-scale\nduplicates" & og$VertType=="Single-copy")], 2, width, breaks, color, 20)
	# SC - D
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Small-scale\nduplicates")], 4, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], og$DiffDom[which(og$BlanType=="Single-copy" & og$VertType=="Small-scale\nduplicates")], 5, width, breaks, color, 20)
	# SC - O
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Ohnologs")], 7, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], og$DiffDom[which(og$BlanType=="Single-copy" & og$VertType=="Ohnologs")], 8, width, breaks, color, 20)
	axis(1, at = c(1.5,4.5,7.5), labels= c("Small-scale\nB. lanceolatum","Small-scale\nD. rerio","Ohnologs\nD. rerio"), line=2, tick=FALSE, las=1, cex.axis=1.2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.2)

	polygon(c(4,4.5,4.5,4),c(0.9,0.9,0.95,0.95), col=color, border=color, lwd=2)
	polygon(c(4,4.5,4.5,4),c(0.8,0.8,0.85,0.85), col=modif_alpha(color,.3), border=color, lwd=2)
	polygon(c(4,4.5,4.5,4),c(0.8,0.8,0.85,0.85), col=color, border=color, density=20, lwd=2)
	text(4.5, 0.92,  labels ="Independent genes", pos=4, cex=1.2)
	text(4.5, 0.82,  labels ="Union of duplicates", pos=4, cex=1.2)
}

DifferenceWithSC_TandemIntraInter <- function(gp, og, tissues, lab1, lab2){
	breaks <- seq(-length(tissues[,1])-.5,length(tissues[,1])+.5,1)
	width <- .7
	color <- "tomato3"
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0-.5, 8+.5), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1)
	mtext(lab2, side = 2, line = 5, cex=1)
	# D - SC Inter
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], gp$DiffDom[which(gp$BlanType=="Small-scale\nduplicates" & gp$VertType=="Single-copy" & gp$BlanLType=="Inter")], 1, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], og$DiffDom[which(og$BlanType=="Small-scale\nduplicates" & og$VertType=="Single-copy" & og$BlanLType=="Inter")], 2, width, breaks, color, 20)
	# D - SC Intra
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], gp$DiffDom[which(gp$BlanType=="Small-scale\nduplicates" & gp$VertType=="Single-copy" & gp$BlanLType=="Intra")], 4, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], og$DiffDom[which(og$BlanType=="Small-scale\nduplicates" & og$VertType=="Single-copy" & og$BlanLType=="Intra")], 5, width, breaks, color, 20)
	# D - SC Tandem
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], gp$DiffDom[which(gp$BlanType=="Small-scale\nduplicates" & gp$VertType=="Single-copy" & gp$BlanLType=="Tandem")], 7, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$BlanType=="Single-copy" & gp$VertType=="Single-copy")], og$DiffDom[which(og$BlanType=="Small-scale\nduplicates" & og$VertType=="Single-copy" & og$BlanLType=="Tandem")], 8, width, breaks, color, 20)
	axis(1, at = c(1.5,4.5,7.5), labels= c("Multichr.","Distant\nmonochr.","Tandem"), line=2, tick=FALSE, las=1, cex.axis=1.2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.2)

	polygon(c(4,4.5,4.5,4),c(0.9,0.9,0.95,0.95), col=color, border=color, lwd=2)
	polygon(c(4,4.5,4.5,4),c(0.8,0.8,0.85,0.85), col=modif_alpha(color,.3), border=color, lwd=2)
	polygon(c(4,4.5,4.5,4),c(0.8,0.8,0.85,0.85), col=color, border=color, density=20, lwd=2)
	text(4.5, 0.92,  labels ="Independent genes", pos=4, cex=1.2)
	text(4.5, 0.82,  labels ="Union of duplicates", pos=4, cex=1.2)
}

Plot_Difference <- function(vecSC, vec, pos, w, b, col, dens=NA){
	hSC <- hist(vecSC, breaks=b, plot=FALSE)
	h <- hist(vec, breaks=b, plot=FALSE)
	difference <- sum(abs(h$counts/length(vec)-hSC$counts/length(vecSC)))
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2),c(0,0,difference, difference), col=modif_alpha(col,.3), border=col, lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2),c(0,0,difference, difference), col=col, border=col, density=dens, lwd=2)
}

CalcVertType <- function(x, vertlist, OG3R){
	og <- x[which(names(x)=="OG")]
	x <- as.numeric(x[vertlist])
	num <- max(x)
	numwof <- max(x[which(vertlist!="Drer")])
	if(num>1){
		if(numwof<=1 & og %in% OG3R){
			num=1
		}else{
			num=2
		}
	}
	return(num)
}







