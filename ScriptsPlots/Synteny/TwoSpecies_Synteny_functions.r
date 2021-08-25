(S1, S2, ohnologsfile, resultsfolder, countsfile, NamesFileAV, gtffolder, ChrLenFolder, SpeciesGRefs[which(ShortSpeciesNames==S1)], SpeciesGRefs[which(ShortSpeciesNames==S2)]){


ProduceFilesS1S2 <- function(S1, S2, ohnologsfile, resultsfolder, countsfile, namesfile, gtffolder, chrlenfolder, splongname1, splongname2){
	# List of onhologs
	ohnologs <- read.table(ohnologsfile, h=F)[,1]
	# Orthologous groups counts (amphioxus + vertebrates)
	counts <- read.table(countsfile, h=F, sep = "\t", row.names=NULL)
	system_out <- system(paste("head -1 ", countsfile, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]*_\\([a-z][a-z][a-z]\\)[A-Za-z\\.0-9_]*/\\1\\2/g'"), intern=T)
	Header <- read.table(text=system_out, h=F, sep = "\t")
	colnames(counts) <- c("OG", as.character(unlist(lapply(Header[1,], as.character))))
	counts$Ohnologs <- rep(FALSE, length(counts[,1]))
	counts$Ohnologs[which(counts$OG %in% ohnologs)] <- rep(TRUE, length(counts[which(counts$OG %in% ohnologs),1]))

	ProduceSpeciesfiles(S1, splongname1, colnames(counts), namesfile, gtffolder, chrlenfolder, resultsfolder)
	ProduceSpeciesfiles(S2, splongname2, colnames(counts), namesfile, gtffolder, chrlenfolder, resultsfolder)

	## Calc num chromosomes, max distance between pairs of copies and num genes between copies
	counts$NumChrs1 <- unlist(lapply(counts$OG, function(x){length(unique(GeneData1$Chr[which(GeneData1$OG==x)]))}))
	counts$MaxDist1 <- unlist(lapply(counts$OG, CalcMaxDist, GeneData=GeneData1))
	counts$BetwGenes1 <- unlist(lapply(counts$OG, CalcNumBetwGenes, GeneData=GeneData1))
	counts$DupType1 <- rep("Missing", length(counts[,1]))
	counts$DupType1[which(counts[,Species1]==1)] <- rep("SingleCopy",sum(counts[,Species1]==1))
	counts$DupType1[which(counts$NumChrs1==1 & counts[,Species1]>1)] <- rep("Intra",sum(counts$NumChrs1==1 & counts[,Species1]>1))
	counts$DupType1[which(counts$NumChrs1>1 & counts[,Species1]>1)] <- rep("Inter",sum(counts$NumChrs1>1 & counts[,Species1]>1))
	counts$DupType1[which(counts$NumChrs1==1 & counts$BetwGenes1==0 & counts$MaxDist1<=TandemMaxDist & counts[,Species1]>1)] <- rep("Tandem",sum(counts$NumChrs1==1 & counts$BetwGenes1==0 & counts$MaxDist1<=TandemMaxDist & counts[,Species1]>1))

	counts$NumChrs2 <- unlist(lapply(counts$OG, function(x){length(unique(GeneData2$Chr[which(GeneData2$OG==x)]))}))
	counts$MaxDist2 <- unlist(lapply(counts$OG, CalcMaxDist, GeneData=GeneData2))
	counts$BetwGenes2 <- unlist(lapply(counts$OG, CalcNumBetwGenes, GeneData=GeneData2))
	counts$DupType2 <- rep("Missing", length(counts[,1]))
	counts$DupType2[which(counts[,Species2]==1)] <- rep("SingleCopy",sum(counts[,Species2]==1))
	counts$DupType2[which(counts$NumChrs2==1 & counts[,Species2]>1)] <- rep("Intra",sum(counts$NumChrs2==1 & counts[,Species2]>1))
	counts$DupType2[which(counts$NumChrs2>1 & counts[,Species2]>1)] <- rep("Inter",sum(counts$NumChrs2>1 & counts[,Species2]>1))
	counts$DupType2[which(counts$NumChrs2==1 & counts$BetwGenes2==0 & counts$MaxDist2<=TandemMaxDist & counts[,Species2]>1)] <- rep("Tandem",sum(counts$NumChrs2==1 & counts$BetwGenes2==0 & counts$MaxDist2<=TandemMaxDist & counts[,Species2]>1))
	write.table(counts, file = paste0(ResultsFolder, "/counts_", Species1, "_", Species2, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)

}

ProduceSpeciesfiles <- function(Sp, splongname, headcounts, namesfile, gtffolder, chrlenfolder, resultsfolder){
	# Read gene location information Species 1
	system_out <- system(paste0("zcat ", gtffolder, "/", splongname, ".gtf.gz | sed 's/\\(.*\\)\\t.*\\t.*\\t\\(.*\\)\\t\\(.*\\)\\t.*\\t.*\\t.*\\t.*gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1\\t\\2\\t\\3\\t\\4/g' | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"), intern=T)
	genecoord1 <- read.table(text=system_out, h=F, sep = "\t")
	colnames(genecoord1) <- c("Gene","Chr","Start","End","Index")
	genecoord1$Midp <- genecoord1$Start + (genecoord1$End - genecoord1$Start)/2

	# OG to Sp gene
	system_out <- system(paste("tail -n +2 ", namesfile, " | awk -F '\\t' '{split($", which(headcounts==Sp), ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
	og2Spgene <- read.table(text=system_out, h=F, sep = "\t")
	colnames(og2Spgene) <- c("OG", "Gene")
	print(head(og2Spgene))
	genedata <- merge(og2Spgene, genecoord1, by="Gene")
	 
	# Read chomosomal length information Species 1
	system_out <- system(paste0("cat ", chrlenfolder, "/", splongname, "/", splongname, "_lengths.txt | sed 's/>//g' | sed 's/chr//g' | grep -P '^[0-9XY]+'"), intern=T)
	chrlen1 <- read.table(text=system_out, h=F, sep = "\t")
	colnames(chrlen1) <- c("Chr","Length","Num","CummLength")
	chrlen1 <- chrlen1[order(chrlen1$Chr),]
	numericchr <- as.numeric(chrlen1$Chr[grep("[0-9]", chrlen1$Chr)])
	XYchr <- chrlen1$Chr[grep("[0-9]", chrlen1$Chr, invert=TRUE)]
	chrlen1 <- chrlen1[match(c(numericchr[order(numericchr)],XYchr), chrlen1$Chr),]
	chrlen1$NumGenes <- unlist(lapply(chrlen1$Chr, function(x){max(c(0,genecoord1$Index[which(genecoord1$Chr==x)]))}))
	chrlen1$CummNumGenes <- cumsum(as.numeric(chrlen1$NumGenes))
	chrlen1$CummLength <- cumsum(as.numeric(chrlen1$Length))

	# Write files
	cat(paste0("Write files ", Sp, "\n"))
	write.table(chrlen1, file = paste0(resultsfolder, "/ChrLen_", Sp, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	write.table(genedata1, file = paste0(resultsfolder, "/GeneData_", Sp, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	return(genedata1)
}

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


SyntenyAlongChrPlot <- function(column, main, GP.1, GP.2=NULL){
	if(column=="Index"){
		maxchrc <- "NumGenes"
		cummmaxchrc <- "CummNumGenes"
	}else if(column=="Midp"){
		maxchrc <- "Length"
		cummmaxchrc <- "CummLength"
	}else{
		quit("Invalid column name in SyntenyAlongChrPlot function")
	}
	axis1 <- chrlen1[,cummmaxchrc]-chrlen1[,maxchrc]/2
	axis2 <- ChrLen2[,cummmaxchrc]-ChrLen2[,maxchrc]/2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,ChrLen2[length(ChrLen2[,1]),cummmaxchrc]), xlim=c(0,chrlen1[length(chrlen1[,1]),cummmaxchrc]), col=NA)
	mtext(main, side = 3, line = 2, cex=2)
	mtext(SpeciesNames[which(ShortSpeciesNames==S1)], side = 1, line = 6, cex=3)
	mtext(SpeciesNames[which(ShortSpeciesNames==Species2)], side = 2, line = 5, cex=3)
	abline(v=chrlen1[,cummmaxchrc], lty=2, lwd=.5, col="grey60")
	abline(h=ChrLen2[,cummmaxchrc], lty=2, lwd=.5, col="grey60")
	Pos2.1 <- GP.2[,paste0(column,"1")]+chrlen1[match(GP.2$Chr1, chrlen1$Chr),cummmaxchrc]-chrlen1[match(GP.2$Chr1, chrlen1$Chr),maxchrc]
	Pos2.2 <- GP.2[,paste0(column,"2")]+ChrLen2[match(GP.2$Chr2, ChrLen2$Chr),cummmaxchrc]-ChrLen2[match(GP.2$Chr2, ChrLen2$Chr),maxchrc]
	points(Pos2.1, Pos2.2, col="grey60", pch=16, cex=.7)
	Pos1.1 <- GP.1[,paste0(column,"1")]+chrlen1[match(GP.1$Chr1, chrlen1$Chr),cummmaxchrc]-chrlen1[match(GP.1$Chr1, chrlen1$Chr),maxchrc]
	Pos1.2 <- GP.1[,paste0(column,"2")]+ChrLen2[match(GP.1$Chr2, ChrLen2$Chr),cummmaxchrc]-ChrLen2[match(GP.1$Chr2, ChrLen2$Chr),maxchrc]
	points(Pos1.1, Pos1.2, col="black", pch=16, cex=.7)
	print(head(Pos1.1))
	axis(1, at = c(0,chrlen1[,cummmaxchrc]), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = c(0,ChrLen2[,cummmaxchrc]), labels=NA, lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = axis1[seq(1, length(axis1), 2)], labels=chrlen1$Chr[seq(1, length(axis1), 2)], lwd=NA, las=1, cex.axis=2.7, line=1.2)
	axis(2, at = axis2[seq(1, length(axis2), 2)], labels=ChrLen2$Chr[seq(1, length(axis2), 2)], lwd=NA, las=1, cex.axis=2.7)
	axis(1, at = axis1[seq(2, length(axis1), 2)], labels=chrlen1$Chr[seq(2, length(axis1), 2)], lwd=NA, las=1, cex.axis=2.7, line=1.2)
	axis(2, at = axis2[seq(2, length(axis2), 2)], labels=ChrLen2$Chr[seq(2, length(axis2), 2)], lwd=NA, las=1, cex.axis=2.7)
	box()	
}



modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

AddHistlogYaxis <- function(h, col){
	for(i in c(1:length(h$counts))){
		polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,log(h$counts[i]),log(h$counts[i])), col=modif_alpha(col,.2), border=col)
	}
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
	return(list(lab1, lab2, Overlap/(group1*group2/Total), Overlap, (group1*group2/Total), dp, ep))
}




