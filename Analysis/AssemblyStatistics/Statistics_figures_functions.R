














modif_alpha <- function(col, alpha=.5){
	if(missing(col))
	stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}


Scatter_plot <- function(vec1, vec2, xlim, ylim, names, col, lab1, lab2){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1],ylim[2]*1.1), xlim=c(xlim[1],xlim[2]*1.1), col=NA)
	mtext(lab1, side = 1, line = 3, cex=1.5)
	mtext(lab2, side = 2, line = 3, cex=1.5)

	points(vec1, vec2, pch=19, col=col, cex=2)
	#points(vec1, vec2, pch=21, bg=modif_alpha(col,0.1), col=col, cex=3)
	text(vec1, vec2, labels=names, pos=3, font=2, cex=1)

	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5),  lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1)
	box()
}



