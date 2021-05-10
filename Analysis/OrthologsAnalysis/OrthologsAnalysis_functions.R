


VerticalDensity <- function(values, pos, col){
	w <- .8
	d  <- density(values, adjust = 0.05)
	polygon(c(pos+d$y/max(d$y)*w/2, rev(pos-d$y/max(d$y)*w/2)), c(d$x,rev(d$x)), col=col, border=col, lwd=2)
	#points(jitter(rep(pos, length(values)), amount=w/2), values, pch=21, bg=col, col=col)
}

PlotaJitterPoints <- function(values, pos, col){
	w <- .8
	points(jitter(rep(pos, length(values)), amount=w/2), values, pch=21, bg=col, col=col)
}

BoxPlot <- function(values, pos, col){
	w <- .8
	s <- boxplot(values, plot=FALSE)
	lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border=col, lwd=2)
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2)
	#points(pos, m, pch=21, bg=col, col=col, cex=4)
}

PlotColumnVertebrateType <- function(df, pos, col){
	len <- length(df[,1])
	base <- 0
	w <- .8
	value <- length(df[which(df$VertebType=="Missing"),1])/len*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=col[1], border=NA)
	base <- base+value
	value <- length(df[which(df$VertebType=="SingleCopy"),1])/len*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=col[2], border=NA)
	base <- base+value
	value <- length(df[which(df$VertebType=="Ohnolog"),1])/len*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=col[3], border=NA)
	base <- base+value
	value <- length(df[which(df$VertebType=="Duplicated"),1])/len*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=col[4], border=NA)
	base <- base+value
}



PlotColumnBlanCN <- function(df, pos, col){
	len <- length(df[,1])
	base <- 0
	w <- .8
	value <- length(df[which(df$BlanType=="Missing"),1])/len*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=col[1], border=NA)
	base <- base+value
	value <- length(df[which(df$BlanType=="SingleCopy"),1])/len*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=col[2], border=NA)
	base <- base+value
	value <- length(df[which(df$BlanType=="Duplicated"),1])/len*100
	polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=col[3], border=NA)
	base <- base+value
}



modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}









