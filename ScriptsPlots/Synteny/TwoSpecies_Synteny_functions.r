
SyntenyScatterPlotColors <- function(S1name, S2name, chrlen1, chrlen2, GP){
	colfunc <- colorRampPalette(c("forestgreen", "gold2", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
	chrlen1$Color <- colfunc(length(chrlen1$Chr))
	GP$ChrColor <- chrlen1$Color[match(GP$Chr1, chrlen1$Chr)]
	axis1 <- chrlen1$CummLength-chrlen1$Length/2
	axis2 <- chrlen2$CummLength-chrlen2$Length/2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,chrlen2$CummLength[length(chrlen2[,1])]), xlim=c(0,chrlen1$CummLength[length(chrlen1[,1])]), col=NA)
	mtext(S1name, side = 1, line = 4, cex=2)
	mtext(S2name, side = 2, line = 3, cex=2)
	abline(v=chrlen1$CummLength, lty=2, lwd=.5, col="grey60")
	abline(h=chrlen2$CummLength, lty=2, lwd=.5, col="grey60")
	points(GP$CummMidp1, GP$CummMidp2, col=GP$ChrColor, pch=16, cex=.4)
	axis(1, at = c(0,chrlen1$CummLength), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = c(0,chrlen2$CummLength), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = axis1[seq(1, length(axis1), 2)], labels=chrlen1$Chr[seq(1, length(axis1), 2)], lwd=NA, las=1, cex.axis=1)
	axis(2, at = axis2[seq(1, length(axis2), 2)], labels=chrlen2$Chr[seq(1, length(axis2), 2)], lwd=NA, las=1, cex.axis=1)
	axis(1, at = axis1[seq(2, length(axis1), 2)], labels=chrlen1$Chr[seq(2, length(axis1), 2)], lwd=NA, las=1, cex.axis=1)
	axis(2, at = axis2[seq(2, length(axis2), 2)], labels=chrlen2$Chr[seq(2, length(axis2), 2)], lwd=NA, las=1, cex.axis=1)
	box()	
}

modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

