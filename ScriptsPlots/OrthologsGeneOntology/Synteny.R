


Synteny <- function(Gene.df, OG.df, chrlens){
	pdf(paste(ResultsFolder, "/Synteny.pdf", sep=""), width=20, height=10)
	par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
	layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)
	chrs <- unique(Gene.df$Chr)[grep("chr", unique(Gene.df$Chr))]

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(.5,length(chrs)+.5), col=NA)
	mtext("Gene density", side = 2, line = 5, cex=2)
	for(chr in c(1:length(chrs))){
		tmp.df <- Gene.df[which(Gene.df$Chr==chrs[chr]),]
		mergedlength <- system("sed 's/\\./\t/g' | awk '{if($1<$2){print $1\"\t\"$2}else{print $2\"\t\"$1}}' | sort -k1,1V | awk '{if(NR==1){start=$1; end=$2; sum=0; next}if($1 > end){sum=sum+end-start; start=$1; end=$2}else{if($2 > end){end=$2}}}END{sum=sum+end-start; print sum}'", intern=TRUE, input=paste(tmp.df[, "Start"], tmp.df[, "End"], sep="."))
		proplength <- as.numeric(mergedlength)/chrlens$Length[which(rownames(chrlens)==chrs[chr])]
		points(chr, proplength, pch=16, col="black", cex=2)
	}
	axis(1, at = c(1:length(chrs)), tick=FALSE, line=2, las=1, cex.axis=1)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,.5), xlim=c(.5,length(chrs)+.5), col=NA)
	mtext("Prop of genes in chr", side = 2, line = 5, cex=2)
	for(chr in c(1:length(chrs))){
		for(blant in c(1:length(BlanTypes))){
			for(vt in c(1:length(VTypes))){
				proplength <- length(Gene.df[which(Gene.df$Chr==chrs[chr] & Gene.df$BlanType.AV==BlanTypes[blant] & Gene.df$VertebType.AV==VTypes[vt]),1])/length(Gene.df[which(Gene.df$Chr==chrs[chr]),1])
				points(chr, proplength, pch=BlanType.pch[blant], col=VTypeColors[vt], cex=2)
			}
		}
	}
	axis(1, at = c(1:length(chrs)), tick=FALSE, line=2, las=1, cex.axis=1)
	axis(2, at = seq(0,1,.1), lwd.ticks=1, las=1, cex.axis=1)
	legend("topright", c(BlanTypes, "", VTypes), pch=c(BlanType.pch, NA, rep(15,length(VTypes))), col=c("black", "black", NA, VTypeColors), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)

	Gene.df$Length <- abs(Gene.df$End-Gene.df$Start)
	dist <- .15
	layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 25000), xlim=c(0.5, 4.5), col=NA)
	mtext("Gene length (kbp)", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
	BoxPlot(Gene.df$Length[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
	axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,25000,5000), labels=seq(0,25000,5000)/1000, lwd.ticks=1, las=1, cex.axis=2)


	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 50), xlim=c(0, 10000), col=NA)
	breaks <- seq(0, 40000000, 50)
	hist(OG.df$MaxDist[which(OGData.AV$NumChrs==1 & OGData.AV$BlanType=="Duplicated")], breaks=breaks, add=T, col="black")
	abline(v=TandemMaxDist, lty=2, col="darkred")
	axis(1, at = seq(0,10000,5000), lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = seq(0,50,5), lwd.ticks=1, las=1, cex.axis=2)

	dist <- .1
	wid <- .2
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 4000), xlim=c(0.5, 4.5), col=NA)
	mtext("Number of orthologous groups", side = 2, line = 7, cex=2)
	for(vt in c(1:length(VTypes))){
		t <- length(OG.df[which(OG.df$DupType=="Tandem" & OG.df==VTypes[vt]),1])
		polygon(c(vt,vt+wid,vt+wid,vt), c(0,0,t,t), col=VTypeColors[vt], border="black")
		a <- length(OG.df[which(OG.df$DupType=="Intra" & OG.df==VTypes[vt]),1])
		polygon(c(vt+dist,vt+wid+dist,vt+wid+dist,vt+dist), c(0,0,a,a), col=VTypeColors[vt], border="black")
		e <- length(OG.df[which(OG.df$DupType=="Inter" & OG.df==VTypes[vt]),1])
		polygon(c(vt+dist*2,vt+wid+dist*2,vt+wid+dist*2,vt+dist*2), c(0,0,e,e), col=VTypeColors[vt], border="black")
	}
	axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,10000,500), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 4.5), col=NA)
	mtext("Proportion of orthologous groups", side = 2, line = 7, cex=2)
	for(vt in c(1:length(VTypes))){
		t <- length(OG.df[which(OG.df$DupType=="Tandem" & OG.df==VTypes[vt]),1])/length(OG.df[which(OG.df$BlanType=="Duplicated" & OG.df==VTypes[vt]),1])
		polygon(c(vt,vt+wid,vt+wid,vt), c(0,0,t,t), col=VTypeColors[vt], border="black")
		a <- length(OG.df[which(OG.df$DupType=="Intra" & OG.df==VTypes[vt]),1])/length(OG.df[which(OG.df$BlanType=="Duplicated" & OG.df==VTypes[vt]),1])
		polygon(c(vt+dist,vt+wid+dist,vt+wid+dist,vt+dist), c(0,0,a,a), col=VTypeColors[vt], border="black")
		e <- length(OG.df[which(OG.df$DupType=="Inter" & OG.df==VTypes[vt]),1])/length(OG.df[which(OG.df$BlanType=="Duplicated" & OG.df==VTypes[vt]),1])
		polygon(c(vt+dist*2,vt+wid+dist*2,vt+wid+dist*2,vt+dist*2), c(0,0,e,e), col=VTypeColors[vt], border="black")
	}
	axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

	OG.blanD <- OG.df[which(OG.df$BlanType=="Duplicated"),]
	HyperTests <- data.frame()
	for(vtype in c("Missing","SingleCopy","Ohnolog","Duplicated")){
		for(dtype in c("Tandem","Intra", "Inter")){
			vnum <- length(OG.blanD[which(OG.blanD$VertebType == vtype),1])
			anum <- length(OG.blanD[which(OG.blanD$DupType == dtype),1])
			onum <- length(OG.blanD[which(OG.blanD$VertebType == vtype & OG.blanD$DupType == dtype),1])
			HyperTests <- rbind(HyperTests, unlist(HypergeometricTest(onum, anum, vnum, length(OG.blanD[,1]), dtype, vtype, PQvalThreshold)))
		}
	}
	colnames(HyperTests) <- c("VertebType", "DupType", "FoldChange", "Observed", "Expected", "HpvalDepleted", "HpvalEnriched")
	HyperTests$HqvalDepleted <- qvalue(as.numeric(HyperTests$HpvalDepleted))$qvalues
	HyperTests$HqvalEnriched <- qvalue(as.numeric(HyperTests$HpvalEnriched))$qvalues
	HyperTests$HResult <- rep("NA", length(HyperTests[,1]))
	HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)] <- rep("D", length(HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)]))
	HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)] <- rep("E", length(HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)]))
	print(HyperTests)


}





