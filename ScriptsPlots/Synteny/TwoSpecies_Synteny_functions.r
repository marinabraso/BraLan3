
SyntenyScatterPlot <- function(S1name, S2name, chrlen1, chrlen2, GP.1, GP.2=NULL, highlightedgenes=NULL){
	axis1 <- chrlen1$CummLength-chrlen1$Length/2
	axis2 <- chrlen2$CummLength-chrlen2$Length/2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,chrlen2$CummLength[length(chrlen2[,1])]), xlim=c(0,chrlen1$CummLength[length(chrlen1[,1])]), col=NA)
	mtext("Gene mid point coordinates", side = 3, line = 2, cex=2)
	mtext(S1name, side = 1, line = 6, cex=3)
	mtext(S2name, side = 2, line = 5, cex=3)
	abline(v=chrlen1$CummLength, lty=2, lwd=.5, col="grey60")
	abline(h=chrlen2$CummLength, lty=2, lwd=.5, col="grey60")
	points(GP.2$CummMidp1, GP.2$CummMidp2, col="grey60", pch=16, cex=.7)
	points(GP.1$CummMidp1, GP.1$CummMidp2, col="black", pch=16, cex=.7)
	points(GP.2$CummMidp1[which(GP.2$Gene1 %in% highlightedgenes)], GP.2$CummMidp2[which(GP.2$Gene1 %in% highlightedgenes)], col="firebrick2", pch=16, cex=.7)
	points(GP.1$CummMidp1[which(GP.1$Gene1 %in% highlightedgenes)], GP.1$CummMidp2[which(GP.1$Gene1 %in% highlightedgenes)], col="firebrick4", pch=16, cex=.7)
	axis(1, at = c(0,chrlen1$CummLength), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = c(0,chrlen2$CummLength), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = axis1[seq(1, length(axis1), 2)], labels=chrlen1$Chr[seq(1, length(axis1), 2)], lwd=NA, las=1, cex.axis=2.7, line=1.2)
	axis(2, at = axis2[seq(1, length(axis2), 2)], labels=chrlen2$Chr[seq(1, length(axis2), 2)], lwd=NA, las=1, cex.axis=2.7)
	axis(1, at = axis1[seq(2, length(axis1), 2)], labels=chrlen1$Chr[seq(2, length(axis1), 2)], lwd=NA, las=1, cex.axis=2.7, line=1.2)
	axis(2, at = axis2[seq(2, length(axis2), 2)], labels=chrlen2$Chr[seq(2, length(axis2), 2)], lwd=NA, las=1, cex.axis=2.7)
	box()	
}

SyntenyLinearPlot <- function(S1name, S2name, chrlen1, chrlen2, GP.1, GP.2=NULL){
	axis1 <- chrlen1$CummLengthMargin-chrlen1$Length/2-Margin/2
	axis2 <- chrlen2$CummLengthMargin-chrlen2$Length/2-Margin/2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0,max(chrlen1$CummLengthMargin[length(chrlen1[,1])], chrlen2$CummLengthMargin[length(chrlen2[,1])])), col=NA)
	#mtext(S1name, side = 1, line = 6, cex=3)
	#mtext(S2name, side = 2, line = 5, cex=3)
	tmp <- lapply(c(1:length(GP.1[,1])), function(x){lines(c(GP.1$CummMidp1Margin[x], GP.1$CummMidp2Margin[x]), c(0,1), col=modif_alpha("black", .3))})
	text(rep(1, length(axis1)), axis1, labels=chrlen1$Chr, las=1, cex.axis=2.7, pos=0)
	axis(1, at = axis2, labels=chrlen2$Chr, lwd=NA, las=1, cex.axis=2.7, pos=1)
}


modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

