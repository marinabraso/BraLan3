


VerticalDensity <- function(values, pos, col){
	w <- .8
	d  <- density(values, adjust = 0.05)
	polygon(c(pos+d$y/max(d$y)*w/2, rev(pos-d$y/max(d$y)*w/2)), c(d$x,rev(d$x)), col=col, border=col, lwd=2)
	#points(jitter(rep(pos, length(values)), amount=w/2), values, pch=21, bg=col, col=col)
}

PlotJitterPoints <- function(values, pos, col){
	w <- .8
	points(jitter(rep(pos, length(values)), amount=w/2), values, pch=21, bg=col, col=col)
}

BoxPlot <- function(values, pos, col, cextext=1, w=.8, den=NULL, text=FALSE){
	s <- boxplot(values, plot=FALSE)
	lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	if(!is.null(den)){
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=modifColor(col, .3), border=col, lwd=4)		
	}
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border=col, lwd=4, density=den)
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=NA, border=col, lwd=4)		
	#points(pos, m, pch=21, bg=col, col=col, cex=4)
	if(text){
		par(xpd=TRUE) 
		text(pos, 0,  labels =length(values), pos=1, cex=cextext)
		par(xpd=FALSE) 		
	}
}

PlotMedianQuartiles <- function(values, pos, col, pch=19, cex=2){
	if(length(values)>0){
		s <- summary(values)
		#arrows(x0=pos, y0=s[2], x1=pos, y1=s[5], code=3, angle=90, length=0.02, col="grey70", lty=2)
		points(pos, s[3], col=col, pch=pch, cex=cex)	
	}
}

PlotPropOfAccelerated <- function(df, pos, totaldf, col, pch=16, cex=3, plus=0.001, QvalT=0.1){
	expected <- as.numeric(unlist(RandomizationTestProportion(length(df[which(df$M8.Qval<=QvalT),1]), length(df[,1]), length(totaldf[which(totaldf$M8.Qval<=QvalT),1]), length(totaldf[,1]), "", "")))
	arrows(x0=pos, y0=expected[6], x1=pos, y1=expected[7], code=3, angle=90, length=0.05, col="black")
	points(pos, expected[5], col="black", bg="white", pch=pch+5, cex=cex-1)
	points(pos, length(df[which(df$M8.Qval<=QvalT),1])/length(df[,1]), col=col, pch=pch, cex=cex)
	#text(pos, 0, labels=length(df[which(df$M8.Qval<=QvalT),1]), pos=3, cex=1)
}

PlotColumnVertebrateType <- function(df, pos, col, cextext, alpha=0, den=NULL){
	len <- length(df[,1])
	base <- 0
	w <- .8
	for(vt in c(1:length(unique(df$VertebType)))){
		value <- length(df[which(df$VertebType==unique(df$VertebType)[vt]),1])/len*100
		if(value>0){
			if(!is.null(den)){
				polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[vt], .1), border=col[vt], lwd=4)		
			}
			polygon(c(pos-w/2,pos+w/2,pos+w/2,pos-w/2), c(base,base,base+value,base+value), col=modifColor(col[vt], alpha), border=col[vt], lwd=4, density=den)
			text(pos, base+value/2, labels=length(df[which(df$VertebType==unique(df$VertebType)[vt]),1]), cex=cextext)
			base <- base+value
		}
	}
	par(xpd=TRUE) 
	text(pos, 100,  labels =length(df[,1]), pos=3, cex=cextext)
	par(xpd=FALSE) 
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

modifColor <- function(col, change){
	apply(sapply(col, col2rgb)/255, 2,
	function(x)
	return(rgb(max(0,min(x[1]+change,1)), max(0,min(x[2]+change,1)), max(0,min(x[3]+change,1)))))
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
	return(list(lab1, lab2, log(Overlap/(group1*group2/Total)), Overlap, (group1*group2/Total), dp, ep))
}

RandomizationTest <- function(Overlap, group1, group2, Total, lab1, lab2, threshold, iterations=1000){
	expDist <- c()
	for(i in c(1:iterations)){
		r1 <- sample(1:Total, group1)
		r2 <- sample(1:Total, group2)
		expDist <- c(expDist, sum(c(1:Total) %in% r1 & c(1:Total) %in% r2))
	}
	ep <- sum(expDist >= Overlap)/iterations
	dp <- sum(expDist <= Overlap)/iterations
	return(list(lab1, lab2, dp, ep, mean(expDist), sd(expDist)))
}

RandomizationTestProportion <- function(Overlap, group1, group2, Total, lab1, lab2, iterations=1000){
	expDist <- c()
	for(i in c(1:iterations)){
		r1 <- sample(1:Total, group1)
		r2 <- sample(1:Total, group2)
		expDist <- c(expDist, sum(c(1:Total) %in% r1 & c(1:Total) %in% r2))
	}
	expDist <- expDist/Total
	ep <- sum(expDist >= Overlap/Total)/iterations
	dp <- sum(expDist <= Overlap/Total)/iterations
	numtail <- 1/100*iterations
	expDist <- expDist[order(expDist)]
	lowertail <- expDist[numtail]
	uppertail <- expDist[iterations-numtail+1] # because R is 1 based
	return(list(lab1, lab2, dp, ep, mean(expDist), lowertail, uppertail))
}



RandomizeDifferenceDupSC <- function(df, thersh, colblan, iterations=1000){
	total <- length(df[,1])
	numD <- length(df[which(df[,colblan]=="Duplicated"),1])
	numSC <- length(df[which(df[,colblan]=="SingleCopy"),1])
	numAE <- length(df[which(df$M8.Qval<=QvalT),1])
	propD <- length(df[which(df[,colblan]=="Duplicated" & df$M8.Qval<=QvalT),1])/numD
	propSC <- length(df[which(df[,colblan]=="SingleCopy" & df$M8.Qval<=QvalT),1])/numSC
	diff <- propD-propSC
	exp.diff <- c()
	for(i in c(1:iterations)){
		rD <- sample(c(1:total), numD)
		rAE <- sample(c(1:total), numAE)
		rpropD <- sum(c(1:total) %in% rD & c(1:total) %in% rAE)/numD
		rpropSC <- sum(!(c(1:total) %in% rD) & c(1:total) %in% rAE)/numSC
		exp.diff <- c(exp.diff, rpropD-rpropSC)
	}
	tails <- quantile(exp.diff, c(thersh, 1-thersh)) 
	return(list(diff, mean(exp.diff), tails))
}


PlotRandomizeDifferenceDupSC <- function(df, ypos, col, cex, colblan){
	rdata <- unlist(RandomizeDifferenceDupSC(df, 0.05, colblan))
	points(rdata[1], ypos, col=col, bg=col, pch=16, cex=cex-1)
	arrows(x0=rdata[3], y0=ypos, x1=rdata[4], y1=ypos, code=3, angle=90, length=0.05, col="black")
	points(rdata[2], ypos, col="black", bg="white", pch=16, cex=max(.5,cex-3))

}





