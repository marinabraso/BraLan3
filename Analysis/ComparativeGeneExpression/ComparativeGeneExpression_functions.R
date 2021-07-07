



AddHistlogYaxis <- function(h, col){
	#for(i in c(1:length(h$counts))){
	#	polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,log(h$counts[i]),log(h$counts[i])), col=modif_alpha(col,.2), border=col)
	#}
	lines(h$mids, log(h$counts), col=col, lwd=4)
	points(h$mids, log(h$counts), col=col, pch=16, cex=3)
}






GetSpearmanCorr <- function(pairs, tissue){
	vec1 <- pairs[,paste0(tissue,1)]
	vec2 <- pairs[,paste0(tissue,2)]
	tmp <- cbind(vec1[which(!is.na(vec1) & !is.na(vec2))], vec2[which(!is.na(vec1) & !is.na(vec2))])
	return(cor(tmp[,1], tmp[,2], method="spearman"))
}





modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}





