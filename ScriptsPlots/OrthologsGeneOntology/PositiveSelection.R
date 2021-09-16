



PositiveSelection <- function(Gene.df, OG.df){

	OG.blan <- OG.df[which(OG.df$BlanType!="Missing"),]

	pdf(paste(ResultsFolder, "/PositiveSelection.pdf", sep=""), width=20, height=10)
	par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
	layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
	mtext("M8 D", side = 2, line = 5, cex=2)
	PlotJitterPoints(Gene.df$M8.D[which(Gene.df$BlanType.AV=="SingleCopy")], 1, VTypeColors[2])
	PlotJitterPoints(Gene.df$M8.D[which(Gene.df$BlanType.AV=="Duplicated")], 2, VTypeColors[4])
	axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
	mtext("M8 q-value", side = 2, line = 5, cex=2)
	PlotJitterPoints(Gene.df$M8.Qval[which(Gene.df$BlanType.AV=="SingleCopy")], 1, VTypeColors[2])
	PlotJitterPoints(Gene.df$M8.Qval[which(Gene.df$BlanType.AV=="Duplicated")], 2, VTypeColors[4])
	axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)


	par(mar=c(15,10,10,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
	layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
	xlim <- 0.1
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, length(VTypes)+1), xlim=c(-xlim, xlim), col=NA)
	mtext("Per gene", side = 3, line = 2, cex=2)
	ypos <- 0.5
	PlotRandomizeDifferenceDupSC(Gene.df, ypos, BaseColor, 5, "BlanType.AV")
	for(v in c(1:length(VTypes))){
		ypos <- ypos+1
		PlotRandomizeDifferenceDupSC(Gene.df[which(Gene.df$VertebType.AV==VTypes[v]),], ypos, VTypeColors[v], 5, "BlanType.AV")
	}
	abline(v=0)
	axis(1, at = seq(-xlim, xlim,0.05), lwd.ticks=1, line=0, las=1, cex.axis=1.5)
	abline(h=1, lty=2, col="grey70")
	par(xpd=TRUE) 
	arrows(x0=xlim/4*3, y0=-1, x1=xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
	text(xlim/6, -1.5, pos=4, labels="+ in duplicated genes", cex=1.5)
	arrows(x0=-xlim/4*3, y0=-1, x1=-xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
	text(-xlim/6, -1.5, pos=2, labels="+ in single copy genes", cex=1.5)
	text(rep(-xlim,4), c(0:length(VTypes))+.5, labels=c("Total", VTypes), pos=4, cex=1.5)
	par(xpd=FALSE) 

	xlim <- 0.1
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, length(VTypes)+1), xlim=c(-xlim, xlim), col=NA)
	mtext("Per orthologous group", side = 3, line = 2, cex=2)
	ypos <- 0.5
	PlotRandomizeDifferenceDupSC(OG.blan, ypos, BaseColor, 5, "BlanType")
	for(v in c(1:length(VTypes))){
		ypos <- ypos+1
		PlotRandomizeDifferenceDupSC(OG.blan[which(OG.blan$VertebType==VTypes[v]),], ypos, VTypeColors[v], 5, "BlanType")
	}
	abline(v=0)
	axis(1, at = seq(-xlim, xlim,0.05), lwd.ticks=1, line=0, las=1, cex.axis=1.5)
	abline(h=1, lty=2, col="grey70")
	par(xpd=TRUE) 
	arrows(x0=xlim/4*3, y0=-1, x1=xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
	text(xlim/6, -1.5, pos=4, labels="+ in duplicated genes", cex=1.5)
	arrows(x0=-xlim/4*3, y0=-1, x1=-xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
	text(-xlim/6, -1.5, pos=2, labels="+ in single copy genes", cex=1.5)
	text(rep(-xlim,4), c(0:length(VTypes))+.5, labels=c("Total", VTypes), pos=4, cex=1.5)
	par(xpd=FALSE) 

	xlim <- 0.1
	step <- 5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 40), xlim=c(-xlim, xlim), col=NA)
	mtext("Mean adult expression", side = 2, line = 5, cex=2)
	for(i in seq(0, 40, step)){
		tmpGene.df <- Gene.df[which(Gene.df$MeanAdult>=i & Gene.df$MeanAdult<(i+step)),]
		PlotRandomizeDifferenceDupSC(tmpGene.df, i+step/2, BaseColor, 5, "BlanType.AV")
	}
	abline(v=0)
	axis(1, at = seq(-xlim, xlim,0.05), lwd.ticks=1, line=0, las=1, cex.axis=1.5)
	par(xpd=TRUE) 
	arrows(x0=xlim/4*3, y0=-5, x1=xlim/4, y1=-5, code=1, length=0.1, col="black", lwd=3)
	text(xlim/6, -7, pos=4, labels="+ in duplicated genes", cex=1.5)
	arrows(x0=-xlim/4*3, y0=-5, x1=-xlim/4, y1=-5, code=1, length=0.1, col="black", lwd=3)
	text(-xlim/6, -7, pos=2, labels="+ in single copy genes", cex=1.5)
	axis(2, at = seq(0,40,step), lwd.ticks=1, las=1, cex.axis=1.5)
	par(xpd=FALSE) 

	xlim <- 0.1
	step <- 5
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(1.5, 5.5), xlim=c(-xlim, xlim), col=NA)
	mtext("Blan copy number", side = 2, line = 5, cex=2)
	for(i in c(2:4)){
		tmpOG.blan <- OG.blan[which(OG.blan$Blan==i | OG.blan$BlanType=="SingleCopy"),]
		PlotRandomizeDifferenceDupSC(tmpOG.blan, i, BaseColor, 5, "BlanType")
	}
	tmpOG.blan <- OG.blan[which(OG.blan$Blan>=5 | OG.blan$BlanType=="SingleCopy"),]
	PlotRandomizeDifferenceDupSC(tmpOG.blan, 5, BaseColor, 5, "BlanType")
	abline(v=0)
	axis(1, at = seq(-xlim, xlim,0.05), lwd.ticks=1, line=0, las=1, cex.axis=1.5)
	par(xpd=TRUE) 
	arrows(x0=xlim/4*3, y0=0.8, x1=xlim/4, y1=0.8, code=1, length=0.1, col="black", lwd=3)
	text(xlim/6, 0.5, pos=4, labels="+ in duplicated genes", cex=1.5)
	arrows(x0=-xlim/4*3, y0=0.8, x1=-xlim/4, y1=0.8, code=1, length=0.1, col="black", lwd=3)
	text(-xlim/6, 0.5, pos=2, labels="+ in single copy genes", cex=1.5)
	axis(2, at = c(2:5), labels=c(c(2:4), "5+"), lwd.ticks=1, line=1, las=1, cex.axis=1.5)
	par(xpd=FALSE) 

	dev.off()
}


