
Scatter_plot <- function(vec1, vec2, xlim, ylim, numaxis, names, col, lab1, lab2){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1],ylim[2]+(ylim[2]-ylim[1])*.1), xlim=c(xlim[1],xlim[2]+(xlim[2]-xlim[1])*.1), col=NA)
	mtext(lab1, side = 1, line = 3, cex=1)
	mtext(lab2, side = 2, line = 3, cex=1)

	points(vec1, vec2, pch=19, col=col, cex=1)
	text(vec1, vec2, labels=names, pos=3, font=2, cex=.7)

	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/numaxis),  lwd.ticks=1, las=1, cex.axis=.7)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/numaxis), lwd.ticks=1, las=1, cex.axis=.7)
}



