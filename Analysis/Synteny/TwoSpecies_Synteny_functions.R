



CalcMaxDist <- function(x, GeneData){
	tmp <- GeneData[which(GeneData$OG==x),]
	tmp <- tmp[order(tmp$Start),]
	maxdist <- 0
	if(length(tmp[,1])>1 & length(unique(tmp$Chr))==1){
		for(g1 in c(1:(length(tmp[,1])-1))){
			maxdist <- max(maxdist, tmp$Start[g1+1]-tmp$End[g1])
		}
	}
	return(maxdist)
}
