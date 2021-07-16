


Amphioxus_Vertebrate_categories <- function(OG.df){
	# Print numbers
	Shared <- length(OG.df[which(OG.df$AmphiType != "Missing" & OG.df$VertebType != "Missing"),1])
	VertebTotal <- length(OG.df[which(OG.df$VertebType != "Missing"),1])
	AmphiTotal <- length(OG.df[which(OG.df$AmphiType != "Missing"),1])
	Control <- length(OG.df[which(OG.df$AmphiType == "Missing" & OG.df$VertebType == "Missing"),1])
	print("Shared VertebTotal AmphiTotal Control")
	print(paste(Shared, VertebTotal, AmphiTotal, Control))
	print("Shared/VertebTotal*100")
	print(Shared/VertebTotal*100)
	print("Shared/AmphiTotal*100")
	print(Shared/AmphiTotal*100)

	for(sp in c(Amphi, Verteb)){
		total <- sum(OG.df[,sp]>0)
		dup <- sum(OG.df[,sp]>1)
		totalG <- sum(OG.df[which(OG.df[,sp]>0),sp])
		dupG <- sum(OG.df[which(OG.df[,sp]>1),sp])
		print(paste(sp, total, dup, dup/total*100, totalG, dupG, dupG/totalG*100))
	}

	# Testing significance coexistence of amphioxus and vertebrate categories
	OG.woM <- OG.df[which(OG.df$BlanType!="Missing" & OG.df$VertebType!="Missing"),]
	VTypes.woM <- VTypes[which(VTypes!="Missing")]
	VTypeColors.woM <- VTypeColors[which(VTypes!="Missing")]
	OG.blan <- OG.df[which(OG.df$BlanType!="Missing"),]

	HyperTests <- data.frame()
	for(vtype in c("Missing","SingleCopy","Ohnolog","Duplicated")){
		for(atype in c("Missing","SingleCopy","Duplicated")){
			vnum <- length(OG.df[which(OG.df$VertebType == vtype),1])
			anum <- length(OG.df[which(OG.df$BlanType == atype),1])
			onum <- length(OG.df[which(OG.df$VertebType == vtype & OG.df$BlanType == atype),1])
			HyperTests <- rbind(HyperTests, unlist(HypergeometricTest(onum, anum, vnum, length(OG.df[,1]), atype, vtype, PQvalThreshold)))
		}
	}
	colnames(HyperTests) <- c("BlanType", "VertebType", "FoldChange", "Observed", "Expected", "HpvalDepleted", "HpvalEnriched")
	HyperTests$HBonfDepleted <- as.numeric(HyperTests$HpvalDepleted)*(length(HyperTests$HpvalDepleted)+length(HyperTests$HpvalEnriched))
	HyperTests$HBonfEnriched <- as.numeric(HyperTests$HpvalEnriched)*(length(HyperTests$HpvalDepleted)+length(HyperTests$HpvalEnriched))
	HyperTests$HResult <- rep("NA", length(HyperTests[,1]))
	HyperTests$HResult[which(HyperTests$HBonfDepleted <= PQvalThreshold)] <- rep("D", length(HyperTests$HResult[which(HyperTests$HBonfDepleted <= PQvalThreshold)]))
	HyperTests$HResult[which(HyperTests$HBonfEnriched <= PQvalThreshold)] <- rep("E", length(HyperTests$HResult[which(HyperTests$HBonfEnriched <= PQvalThreshold)]))
	print(HyperTests[,c(1,2,8,9,10)])

	HyperTests.woM <- data.frame()
	for(vtype in c("SingleCopy","Ohnolog","Duplicated")){
		for(atype in c("SingleCopy","Duplicated")){
			vnum <- length(OG.woM[which(OG.woM$VertebType == vtype),1])
			anum <- length(OG.woM[which(OG.woM$BlanType == atype),1])
			onum <- length(OG.woM[which(OG.woM$VertebType == vtype & OG.woM$BlanType == atype),1])
			HyperTests.woM <- rbind(HyperTests.woM, unlist(HypergeometricTest(onum, anum, vnum, length(OG.woM[,1]), atype, vtype, PQvalThreshold)))
		}
	}
	colnames(HyperTests.woM) <- c("BlanType", "VertebType", "FoldChange", "Observed", "Expected", "HpvalDepleted", "HpvalEnriched")
	HyperTests.woM$HBonfDepleted <- as.numeric(HyperTests.woM$HpvalDepleted)*(length(HyperTests.woM$HpvalDepleted)+length(HyperTests.woM$HpvalEnriched))
	HyperTests.woM$HBonfEnriched <- as.numeric(HyperTests.woM$HpvalEnriched)*(length(HyperTests.woM$HpvalDepleted)+length(HyperTests.woM$HpvalEnriched))
	HyperTests.woM$HResult <- rep("NA", length(HyperTests.woM[,1]))
	HyperTests.woM$HResult[which(HyperTests.woM$HBonfDepleted <= PQvalThreshold)] <- rep("D", length(HyperTests.woM$HResult[which(HyperTests.woM$HBonfDepleted <= PQvalThreshold)]))
	HyperTests.woM$HResult[which(HyperTests.woM$HBonfEnriched <= PQvalThreshold)] <- rep("E", length(HyperTests.woM$HResult[which(HyperTests.woM$HBonfEnriched <= PQvalThreshold)]))
	print(HyperTests.woM[,c(1,2,8,9,10)])

	# Plotting
	pdf(paste(ResultsFolder, "/Amphioxus_Vertebrate_categories.pdf", sep=""), width=20, height=10)
	par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')

	## Bar plot of Blan/Verteb categories
	layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,ncol=4,byrow=T), widths=c(1), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, 3.5), col=NA)
	mtext("% in each category", side = 2, line = 5, cex=1.5)
	PlotColumnVertebrateType(OG.df[which(OG.df$Blan==1),], 1, VTypeColors, 1.5)
	PlotColumnVertebrateType(OG.df[which(OG.df$Blan>=2),], 2, VTypeColors, 1.5, den=10)
	PlotColumnVertebrateType(OG.df[which(OG.df$Blan==0),], 3, VTypeColors, 1.5, alpha=.2)
	axis(1, at = c(1:3), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3", "Missing\nin Bralan3"), tick=FALSE, line=2, las=1, cex.axis=1.5)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)

	plot.new()
	legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VTypeColors), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)

	## 2nd Bar plot of Blan/Verteb categories
	#layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,65), xlim=c(0.5, 4.5), col=NA)
	mtext("% in each category", side = 2, line = 5, cex=2)
	w <- .25
	for(vt in c(1:length(VTypes))){
		scvalue <- length(OG.df[which(OG.df$Blan==1 & OG.df$VertebType==VTypes[vt]),1])/length(OG.df[which(OG.df$Blan==1),1])*100
		dvalue <- length(OG.df[which(OG.df$Blan>=2 & OG.df$VertebType==VTypes[vt]),1])/length(OG.df[which(OG.df$Blan>=2),1])*100
		mvalue <- length(OG.df[which(OG.df$Blan==0 & OG.df$VertebType==VTypes[vt]),1])/length(OG.df[which(OG.df$Blan==0),1])*100
		polygon(c(vt-w*1.5, vt-w*.5, vt-w*.5, vt-w*1.5),c(0,0,scvalue, scvalue), col=VTypeColors[vt], border=VTypeColors[vt], lwd=4)		
		polygon(c(vt-w*.5, vt+w*.5, vt+w*.5, vt-w*.5),c(0,0,dvalue, dvalue), col=modifColor(VTypeColors[vt], .1), border=VTypeColors[vt], lwd=4)		
		polygon(c(vt-w*.5, vt+w*.5, vt+w*.5, vt-w*.5),c(0,0,dvalue, dvalue), col=VTypeColors[vt], border=VTypeColors[vt], lwd=4, den=10)		
		polygon(c(vt+w*.5, vt+w*1.5, vt+w*1.5, vt+w*.5),c(0,0,mvalue, mvalue), col=modifColor(VTypeColors[vt], .2), border=VTypeColors[vt], lwd=4)		
	}
	axis(2, at = seq(0,60,10), lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = seq(1,4,1), labels=VTypes, lwd=NA, las=1, cex.axis=1.5)

	colfunc <- colorRampPalette(c("indianred4", "indianred3", "indianred1", "white", "gold1", "gold3", "gold4"))
	gradientcolors <- colfunc(50)
	gradientlims <- c(-2,2)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,10), xlim=c(0,6), col=NA)
	mtext("B.lanceolatum", side = 2, line = 5, cex=1.5)
	mtext("Vertebrates", side = 1, line = 5, cex=1.5)
	for(v in c(1:length(unique(HyperTests$VertebType)))){
		for(a in c(1:length(unique(HyperTests$BlanType)))){
			color <- findColorInGRadient(HyperTests$FoldChange[which(HyperTests$BlanType==unique(HyperTests$BlanType)[a] & HyperTests$VertebType==unique(HyperTests$VertebType)[v])], gradientlims, gradientcolors)
			if(color=="grey30"){dens <- 10}else{dens <- NULL}
			polygon(c(v-1,v,v,v-1), c(a-1,a-1,a,a), col=color, border="white", density=dens)
			if(color!="grey30" & HyperTests$HResult[which(HyperTests$BlanType==unique(HyperTests$BlanType)[a] & HyperTests$VertebType==unique(HyperTests$VertebType)[v])]!="NA"){
					text(v-.5, a-.5, label="*", cex=1.5)
			}
		}
	}
	axis(1, at = seq(0.5,length(unique(HyperTests$VertebType)),1), labels=unique(HyperTests$VertebType), lwd=NA, col = NA, las=1, cex.axis=1)
	axis(2, at = seq(0.5,length(unique(HyperTests$BlanType)),1), labels=unique(HyperTests$BlanType), lwd=NA, col = NA, las=1, cex.axis=1)
	printgradientlegend(c(5.1,1), .03, .2, "log(FC)", gradientlims, gradientcolors)

####################################

	## Bar plot of Blan/Verteb categories
	layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,ncol=4,byrow=T), widths=c(1), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-0.5,100.5), xlim=c(0.5, 3.5), col=NA)
	mtext("% in each category", side = 2, line = 5, cex=1.5)
	mtext("B.lanceolatum", side = 1, line = 5, cex=1.5)
	PlotColumnVertebrateType(OG.woM[which(OG.woM$Blan==1),], 1, VTypeColors, 1.5)
	PlotColumnVertebrateType(OG.woM[which(OG.woM$Blan>=2),], 2, VTypeColors, 1.5, den=10)
	axis(1, at = c(1:2), labels=c("Single copy", "Duplicated"), tick=FALSE, line=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)

	plot.new()
	legend("bottomright", c("Vertebrates", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VTypeColors), bty = "n", pt.cex=2, cex=2, xjust = 0, yjust = 0)

	## 2nd Bar plot of Blan/Verteb categories
	#layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,65), xlim=c(0.5, 4.5), col=NA)
	mtext("% in each category", side = 2, line = 5, cex=1.5)
	w <- .35
	for(vt in c(1:length(VTypes.woM))){
		scvalue <- length(OG.woM[which(OG.woM$Blan==1 & OG.woM$VertebType==VTypes.woM[vt]),1])/length(OG.woM[which(OG.woM$Blan==1),1])*100
		dvalue <- length(OG.woM[which(OG.woM$Blan>=2 & OG.woM$VertebType==VTypes.woM[vt]),1])/length(OG.woM[which(OG.woM$Blan>=2),1])*100
		polygon(c(vt-w, vt, vt, vt-w),c(0,0,scvalue, scvalue), col=VTypeColors.woM[vt], border=VTypeColors.woM[vt], lwd=4)		
		polygon(c(vt, vt+w, vt+w, vt),c(0,0,dvalue, dvalue), col=modifColor(VTypeColors.woM[vt], .1), border=VTypeColors.woM[vt], lwd=4)		
		polygon(c(vt, vt+w, vt+w, vt),c(0,0,dvalue, dvalue), col=VTypeColors.woM[vt], border=VTypeColors.woM[vt], lwd=4, den=10)		
	}
	axis(2, at = seq(0,60,10), lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = seq(1,length(VTypes.woM),1), labels=VTypes.woM, lwd=NA, las=1, cex.axis=1.5)

	colfunc <- colorRampPalette(c("indianred4", "indianred3", "indianred1", "white", "gold1", "gold3", "gold4"))
	gradientcolors <- colfunc(50)
	gradientlims <- c(-2,2)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,7), xlim=c(0,4), col=NA)
	mtext("B.lanceolatum", side = 2, line = 5, cex=1.5)
	mtext("Vertebrates", side = 1, line = 5, cex=1.5)
	for(v in c(1:length(unique(HyperTests.woM$VertebType)))){
		for(a in c(1:length(unique(HyperTests.woM$BlanType)))){
			color <- findColorInGRadient(HyperTests.woM$FoldChange[which(HyperTests.woM$BlanType==unique(HyperTests.woM$BlanType)[a] & HyperTests.woM$VertebType==unique(HyperTests.woM$VertebType)[v])], gradientlims, gradientcolors)
			if(color=="grey30"){dens <- 10}else{dens <- NULL}
			polygon(c(v-1,v,v,v-1), c(a-1,a-1,a,a), col=color, border="white", density=dens)
			if(color!="grey30" & HyperTests.woM$HResult[which(HyperTests.woM$BlanType==unique(HyperTests.woM$BlanType)[a] & HyperTests.woM$VertebType==unique(HyperTests.woM$VertebType)[v])]!="NA"){
					text(v-.5, a-.5, label="*", cex=1.5)
			}
		}
	}
	axis(1, at = seq(0.5,length(unique(HyperTests.woM$VertebType)),1), labels=unique(HyperTests.woM$VertebType), lwd=NA, col = NA, las=1, cex.axis=1)
	axis(2, at = seq(0.5,length(unique(HyperTests.woM$BlanType)),1), labels=unique(HyperTests.woM$BlanType), lwd=NA, col = NA, las=1, cex.axis=1)
	printgradientlegend(c(2,4), .03, .2, "log(FC)", gradientlims, gradientcolors)


	dev.off()
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








