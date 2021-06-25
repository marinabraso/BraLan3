



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


CalcNumBetwGenes <- function(x, GeneData){
	tmp <- GeneData[which(GeneData$OG==x),]
	tmp <- tmp[order(tmp$Start),]
	numBG <- NA
	if(length(tmp[,1])>1 & length(unique(tmp$Chr))==1){
		numBG <- 0
		for(g1 in c(1:(length(tmp[,1])-1))){
			numBG <- numBG + length(GeneData[which(GeneData$End < tmp$Start[g1+1] & GeneData$Start > tmp$End[g1]),1])
		}
	}
	return(numBG)
}


RetrieveGenesInOtherSpecies <- function(x, GeneData, otherGeneData){
	xog <- unique(GeneData[which(GeneData[,"Gene"]==x), "OG"])
	otherx <- c()
	for(g in xog){
		otherx <- c(otherx, unique(otherGeneData[which(otherGeneData[,"OG"]==g),"Gene"]))
	} 
	if(length(otherx)==0){otherx <- NA}		
	return(paste0(x, "_", otherx))
}















