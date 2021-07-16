


ScatterPercentagePlot <- function(vec1, vec2, col, lab1, lab2, xlim, ylim){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(ylim[1]-3, ylim[2]+3), xlim=c(xlim[1]-3, xlim[2]+3), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.5)
	mtext(lab2, side = 2, line = 5, cex=1.5)
	points(vec1, vec2, col=col, bg=modif_alpha(col,.5), pch=21, cex=3)
	#text(50,5,labels=bquote("Spearman correlation = " ~ .(format(round(cor(na.omit(vec1), na.omit(vec2), method="spearman")[1], 2), nsmall = 3))), cex=2)
	axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=1.5)
}


modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}
