



prepare_DrerGeneData <- function(file, bgeefile, sname, tissues, rfolder, pfolder){
	if(file.exists(file)){
		Data <- read.delim(file, h=T, stringsAsFactors=F)
	}else{
		library(BgeeDB)
		setwd(paste0("./", rfolder))
		DrerBgee <- Bgee$new(species=sname, dataType="rna_seq")
		DrerBgeeData <- getData(DrerBgee)
		write.table(DrerBgeeData, file = bgeefile, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
		Data <- as.data.frame(cbind(unique(DrerBgeeData$Gene.ID)))
		colnames(Data) <- c("Gene")
		for(i in c(1:length(tissues[,1]))){
			print(tissues$Name[i])
			system_out <- system(paste0("cat ", bgeefile, " | sed 's/\"//g' | awk -F'\\t' -v tissue=\"", tissues$Drer[i], "\" '{if($7==tissue){print $5\"\\t\"$13\"\\t\"$16}}' | sort -k1,1 | awk '{if(g!=$1){if(NR!=1){mTPM=mTPM/num;} print g\"\\t\"mTPM\"\\t\"pres; g=$1; num=1;mTPM=$2;pres=$3}else{num=num+1;mTPM=mTPM+$2;if($3==\"present\"){pres=$3}}}END{mTPM=mTPM/num; print g\"\\t\"mTPM\"\\t\"pres;}' | tail -n +2"), intern=T)
			tData <- read.table(text=system_out, h=F, sep = "\t")
			Data[,paste0(tissues$Name[i],"TPM")] <- tData[match(Data$Gene, tData[,1]),2]
			Data[,paste0(tissues$Name[i],"Presence")] <- tData[match(Data$Gene, tData[,1]),3]
		}
		write.table(Data, file = file, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
		setwd(pwd)
	}
	return(Data)
}

prepare_BlanGeneData <- function(file, MetaInfoFile, tpmfile){
	if(file.exists(file)){
		Data <- read.delim(file, h=T, stringsAsFactors=F)
	}else{
		system_out <- system(paste("cat ", MetaInfoFile," | cut -f1,13,24", sep=""), intern=T)
		MetaRNA <- read.table(text=system_out, h=F, sep = "\t")
		colnames(MetaRNA) <- c("Sample","Age","Tissue")
		MetaRNA$Tissue <- sub(" ", ".", MetaRNA$Tissue)
		MetaRNA$Age <- sub(" ", ".", MetaRNA$Age)
		MetaRNA$Age <- sub("-", ".", MetaRNA$Age)

		TPM <- read.table(tpmfile, h=T)
		TPM <- TPM[,MetaRNA$Sample[which(MetaRNA$Tissue!="egg" & MetaRNA$Sample%in%colnames(TPM))]]
		MetaRNA <- MetaRNA[which(MetaRNA$Sample %in% colnames(TPM)),]

	Tissues <- unique(MetaRNA$Tissue[which(MetaRNA$Age=="adult")])
	Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads")
	EmbAges <- unique(MetaRNA$Age[which(MetaRNA$Age!="adult")])
	EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre.metamorphic.larvae")

		# TAU calculation per adult tissue
		GeneData <- data.frame("Gene" = rownames(TPM))
		for(t in c(1:length(Tissues))){
			samples <- MetaRNA$Sample[which(MetaRNA$Tissue==Tissues[t])]
			if(length(samples) == 1){
				GeneData <- cbind(GeneData, as.numeric(TPM[, samples]))
			}else{
				GeneData <- cbind(GeneData, as.numeric(rowMeans(TPM[, samples])))	
			}
		}
		for(a in c(1:length(EmbAges))){
			samples <- MetaRNA$Sample[which(MetaRNA$Age==EmbAges[a])]
			if(length(samples) == 1){
				GeneData <- cbind(GeneData, as.numeric(TPM[, samples]))
			}else{
				GeneData <- cbind(GeneData, as.numeric(rowMeans(TPM[, samples])))	
			}
		}
		colnames(GeneData) <- c("Gene", Tissues, EmbAges)
		m.GeneData <- as.matrix(GeneData[,Tissues])
		Tau <- apply(log2(m.GeneData+1), 1, compute.tau)
		MaxGroup <- apply(log2(m.GeneData+1), 1, max.col)
		GeneData$TauTissues <- Tau
		GeneData$MaxTissue <- MaxGroup
		m.GeneData <- as.matrix(GeneData[,EmbAges])
		Tau <- apply(log2(m.GeneData+1), 1, compute.tau)
		MaxGroup <- apply(log2(m.GeneData+1), 1, max.col)
		GeneData$TauEmbAge <- Tau
		GeneData$MaxEmbAge <- MaxGroup

		GeneData$MeanAdult <- rowMeans(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age=="adult")]])
		GeneData$MaxAdult <- apply(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age=="adult")]],1,max)
		GeneData$MeanEmbr <- rowMeans(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age!="adult")]])
		GeneData$MaxEmbr <- apply(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age!="adult")]],1,max)
		return(GeneData)
		write.table(Data, file = file, quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	}
}

compute.tau <- function(exp){
	if(max(exp)==0){
	  return(NA)
	}
	n=length(exp)
	newexp=exp/max(exp)
	tau=sum(1-newexp)/(n-1)
	return(tau)
}

max.col <- function(exp){
	if(max(exp)==0){
	  return(NA)
	}
	maxcollist <- paste(names(exp[which(exp==max(exp))]), sep=":")
	return(maxcollist)
}


BoxPlot_BlanTypes_VertTypes <- function(Genes, GCN, column, ylab, ylim, vcolors){
	dist <- .15
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0.5, 4.5), col=NA)
	mtext(ylab, side = 2, line = 5, cex=1.5)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="SC" & GCN$TypeVert=="M")]), column], 1-dist, vcolors[1], .7, .4)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="D" & GCN$TypeVert=="M")]), column], 1+dist, vcolors[1], .7, .4, 10)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="SC" & GCN$TypeVert=="SC")]), column], 2-dist, vcolors[2], .7, .4)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="D" & GCN$TypeVert=="SC")]), column], 2+dist, vcolors[2], .7, .4, 10)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="SC" & GCN$TypeVert=="O")]), column], 3-dist, vcolors[3], .7, .4)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="D" & GCN$TypeVert=="O")]), column], 3+dist, vcolors[3], .7, .4, 10)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="SC" & GCN$TypeVert=="D")]), column], 4-dist, vcolors[4], .7, .4)
	BoxPlot(Genes[which(Genes$Gene %in% GCN$Gene[which(GCN$TypeBlan=="D" & GCN$TypeVert=="D")]), column], 4+dist, vcolors[4], .7, .4, 10)
	#axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5), lwd.ticks=1, las=1, cex.axis=2)

	xlab <- c("Missing", "Single copy", "Ohnologs", "Duplicated")
	text(x = c(1:4),
		y = par("usr")[3] - (ylim[2]-ylim[1])/40,
		labels = xlab,
		xpd = NA,
		srt = 35,
		cex = 1.5,
		adj = .9)

}

BoxPlot <- function(values, pos, col, cextext=1, w=.8, den=NULL, text=FALSE){
	s <- boxplot(values, plot=FALSE)
	lines(c(pos, pos),c(s$stats[1], s$stats[5]), lwd=2)
	if(!is.null(den)){
		polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=modifColor(col, .3), border=col, lwd=3)		
	}
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=col, border=col, lwd=3, density=den)
	lines(c(pos-w/2, pos+w/2),c(s$stats[3], s$stats[3]), lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2), c(s$stats[2],s$stats[2],s$stats[4],s$stats[4]), col=NA, border=col, lwd=3)		
	if(text){
		par(xpd=TRUE) 
		text(pos, 0,  labels =length(values), pos=1, cex=cextext)
		par(xpd=FALSE) 		
	}
}


modifColor <- function(col, change){
	apply(sapply(col, col2rgb)/255, 2,
	function(x)
	return(rgb(max(0,min(x[1]+change,1)), max(0,min(x[2]+change,1)), max(0,min(x[3]+change,1)))))
}

modif_alpha <- function(col, alpha=.5){
  if(missing(col))
  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}


Hist_ExpressDomanis <- function(gp, og, t1, t2, tissues, lab1, lab2, main){
	breaks <- seq(-length(tissues[,1])-.5,length(tissues[,1])+.5,1)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(-length(tissues[,1]),length(tissues[,1])), col=NA)
	mtext(main, side = 3, line = 4, cex=1.5)
	mtext(lab1, side = 1, line = 6, cex=1.5)
	mtext(lab2, side = 2, line = 5, cex=1.5)
	add_relative_hist(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], breaks, "cornflowerblue")
	add_relative_hist(gp$DiffDom[which(gp$Type1==t1 & gp$Type2==t2)], breaks, "tomato3")
	add_relative_hist(og$DiffDom[which(og$Type1==t1 & og$Type2==t2)], breaks, "tomato3", 20)

	axis(1, at = seq(-length(tissues[,1]),length(tissues[,1]),2), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.5)
}

add_relative_hist <- function(vec, b, col, dens=NULL){
	h <- hist(vec, breaks=b, plot=FALSE)
	h$counts <- h$counts/length(vec)
	for(i in c(1:length(h$counts))){
		polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,h$counts[i],h$counts[i]), col=modif_alpha(col,.3), border=col, density=dens, lwd=2)
	}
	return(h$mids)
}


ScatterPlot_ExpressDomanis <- function(gp, og, tissues, lab1, lab2, lim){
	maxcex <- 10
	relatTableGP <- table(gp[,c("BlanSum","DrerSum")])/max(table(gp[,c("BlanSum","DrerSum")]))
	print(relatTableGP)
	print(relatTableGP[0,0])
	relatTableOG <- table(og[,c("BlanSum","DrerSum")])/max(table(og[,c("BlanSum","DrerSum")]))
	print(relatTableOG)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(lim[1]-1, lim[2]+1), xlim=c(lim[1]-.5, lim[2]+.5), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1.5)
	mtext(lab2, side = 2, line = 5, cex=1.5)
	for(t1 in c(0:length(tissues$Name))){
		for(t2 in c(0:length(tissues$Name))){
			points(t1, t2, col="darkred", bg=modif_alpha("darkred",.5), pch=21, cex=relatTableGP[t1+1,t2+1]*maxcex)
			points(t1, t2, col="forestgreen", bg=modif_alpha("forestgreen",.5), pch=21, cex=relatTableOG[t1+1,t2+1]*maxcex)
		}
	}
	axis(1, at = c(lim[1]:lim[2]), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = c(lim[1]:lim[2]), lwd.ticks=1, las=1, cex.axis=1.5)
}


ExprProfileLines <- function(gp, og, tissues, lab1, lab2, vtypes, cols){
	breaks <- seq(-.5, length(tissues$Name)+.5, 1)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(min(breaks), max(breaks)), col=NA)
	Draw_ExpProfileLine(gp$BlanSum[which(gp$Type1=="SC")], breaks, cols[2], 21)
	Draw_ExpProfileLine(gp$DrerSum[which(gp$Type2=="SC")], breaks, cols[2], 23)
	Draw_ExpProfileLine(gp$BlanSum[which(gp$Type1=="D")], breaks, cols[4], 21)
	Draw_ExpProfileLine(gp$DrerSum[which(gp$Type2=="D")], breaks, cols[4], 23)
	Draw_ExpProfileLine(gp$DrerSum[which(gp$Type2=="O")], breaks, cols[3], 23)
	Draw_ExpProfileLine(og$BlanSum[which(og$Type1=="D")], breaks, cols[4], 21, 2)
	Draw_ExpProfileLine(og$DrerSum[which(og$Type2=="D")], breaks, cols[4], 23, 2)
	Draw_ExpProfileLine(og$DrerSum[which(og$Type2=="O")], breaks, cols[3], 23, 2)
	axis(1, at = c(0:length(tissues$Name)), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.5)
	legend("topleft", c(vtypes[which(vtypes!="Missing")], "B. lanceolatum", "D. rerio", "Independent genes", "Union of duplicates"), lty=c(NA,NA,NA,NA,NA,1,2), pch=c(22,22,22,21,23,NA,NA), text.col="black", col=c(cols[which(vtypes!="Missing")], "black", "black", "black", "black"), pt.bg=modifColor(c(cols[which(vtypes!="Missing")], "black", "black", "black", "black"),.1), bty = "n", pt.cex=1.5, cex=1, xjust = 0, yjust = 0)
}

Draw_ExpProfileLine <- function(vec, b, col, pch, lty=1){
	h <- hist(vec, breaks=b, plot=FALSE)
	lines(h$mids, h$counts/length(vec), col=col, lwd=1.5, lty=lty)
	points(h$mids, h$counts/length(vec), col=col, bg=modifColor(col,.1), pch=pch, cex=1.5)
}








