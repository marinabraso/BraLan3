



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
	mtext(main, side = 3, line = 2, cex=1.2)
	mtext(lab1, side = 1, line = 5, cex=1.2)
	mtext(lab2, side = 2, line = 5, cex=1.2)
	add_relative_hist(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], breaks, "cornflowerblue")
	add_relative_hist(gp$DiffDom[which(gp$Type1==t1 & gp$Type2==t2)], breaks, "tomato3")
	add_relative_hist(og$DiffDom[which(og$Type1==t1 & og$Type2==t2)], breaks, "tomato3", 20)

	axis(1, at = seq(-length(tissues[,1]),length(tissues[,1]),2), lwd.ticks=1, las=1, cex.axis=1.5)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.5)
}

Separate_Hist_ExpressDomanis <- function(gp, og, t1, t2, tissues, lab1, lab2, main){
	breaks <- seq(-length(tissues[,1])-.5,length(tissues[,1])+.5,1)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(-length(tissues[,1]),length(tissues[,1])), col=NA)
	mtext(main, side = 3, line = 2, cex=0.7)
	mtext(lab1, side = 1, line = 3, cex=0.7)
	mtext(lab2, side = 2, line = 3, cex=0.7)
	add_relative_hist(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], breaks, "cornflowerblue")
	axis(1, at = seq(-length(tissues[,1]),length(tissues[,1]),2), lwd.ticks=1, las=1, cex.axis=0.8)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=0.8)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(-length(tissues[,1]),length(tissues[,1])), col=NA)
	mtext(main, side = 3, line = 2, cex=0.7)
	mtext(lab1, side = 1, line = 3, cex=0.7)
	mtext(lab2, side = 2, line = 3, cex=0.7)
	add_relative_hist(gp$DiffDom[which(gp$Type1==t1 & gp$Type2==t2)], breaks, "tomato3")
	axis(1, at = seq(-length(tissues[,1]),length(tissues[,1]),2), lwd.ticks=1, las=1, cex.axis=0.8)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=0.8)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(-length(tissues[,1]),length(tissues[,1])), col=NA)
	mtext(main, side = 3, line = 2, cex=0.7)
	mtext(lab1, side = 1, line = 3, cex=0.7)
	mtext(lab2, side = 2, line = 3, cex=0.7)
	add_relative_hist(og$DiffDom[which(og$Type1==t1 & og$Type2==t2)], breaks, "tomato3", 20)
	axis(1, at = seq(-length(tissues[,1]),length(tissues[,1]),2), lwd.ticks=1, las=1, cex.axis=0.8)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=0.8)
}

add_relative_hist <- function(vec, b, col, dens=NULL){
	h <- hist(vec, breaks=b, plot=FALSE)
	h$counts <- h$counts/length(vec)
	for(i in c(1:length(h$counts))){
		polygon(c(h$breaks[i],h$breaks[i+1],h$breaks[i+1],h$breaks[i]),c(0,0,h$counts[i],h$counts[i]), col=modif_alpha(col,.3), border=col, density=dens, lwd=2)
	}
	return(h$mids)
}

DifferenceWithSC <- function(gp, og, tissues, lab1, lab2){
	breaks <- seq(-length(tissues[,1])-.5,length(tissues[,1])+.5,1)
	width <- .7
	color <- "black"
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0-.5, 8+.5), col=NA)
	mtext(lab1, side = 1, line = 6, cex=1)
	mtext(lab2, side = 2, line = 5, cex=1)
	print("D - SC")
	Plot_Difference(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], gp$DiffDom[which(gp$Type1=="D" & gp$Type2=="SC")], 1, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], og$DiffDom[which(og$Type1=="D" & og$Type2=="SC")], 2, width, breaks, color, 20)
	print("SC - D")
	Plot_Difference(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="D")], 4, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], og$DiffDom[which(og$Type1=="SC" & og$Type2=="D")], 5, width, breaks, color, 20)
	print("SC - O")
	Plot_Difference(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="O")], 7, width, breaks, color)
	Plot_Difference(gp$DiffDom[which(gp$Type1=="SC" & gp$Type2=="SC")], og$DiffDom[which(og$Type1=="SC" & og$Type2=="O")], 8, width, breaks, color, 20)
	axis(1, at = c(1.5,4.5,7.5), labels= c("Small scale\nB. lanceolatum","Small scale\nD. rerio","Ohnologs\nD. rerio"), line=2, tick=FALSE, las=1, cex.axis=1.2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=1.2)

	polygon(c(4,4.5,4.5,4),c(0.9,0.9,0.95,0.95), col=color, border=color, lwd=2)
	polygon(c(4,4.5,4.5,4),c(0.8,0.8,0.85,0.85), col=modif_alpha(color,.3), border=color, lwd=2)
	polygon(c(4,4.5,4.5,4),c(0.8,0.8,0.85,0.85), col=color, border=color, density=20, lwd=2)
	text(4.5, 0.92,  labels ="Independent genes", pos=4, cex=1.2)
	text(4.5, 0.82,  labels ="Union of duplicates", pos=4, cex=1.2)
}


Plot_Difference <- function(vecSC, vec, pos, w, b, col, dens=NA){
	hSC <- hist(vecSC, breaks=b, plot=FALSE)
	h <- hist(vec, breaks=b, plot=FALSE)
	difference <- sum(abs(h$counts/length(vec)-hSC$counts/length(vecSC)))
	print(abs(h$counts/length(vec)-hSC$counts/length(vecSC)))
	print(difference)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2),c(0,0,difference, difference), col=modif_alpha(col,.3), border=col, lwd=2)
	polygon(c(pos-w/2, pos+w/2, pos+w/2, pos-w/2),c(0,0,difference, difference), col=col, border=col, density=dens, lwd=2)
}









