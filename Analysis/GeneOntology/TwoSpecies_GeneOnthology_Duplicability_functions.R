


ScatterPercentagePlot <- function(vec1, vec2, col, lab1, lab2){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,103), xlim=c(0,103), col=NA)
	mtext(lab1, side = 1, line = 5, cex=1.5)
	mtext(lab2, side = 2, line = 5, cex=1.5)
	points(vec1, vec2, col=col, pch=16, cex=2)
	text(50,5,labels=bquote("Spearman correlation = " ~ .(format(round(cor(na.omit(vec1), na.omit(vec2), method="spearman")[1], 2), nsmall = 3))), cex=2)
	axis(1, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
}
