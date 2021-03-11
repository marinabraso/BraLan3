




HistogramGroupSize <- function(values, breaks, xlim, ylim, main){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=xlim, col=NA)
	mtext("Group size", side = 1, line = 3, cex=1.5)
	mtext("Counts", side = 2, line = 3, cex=1.5)
	mtext(main, side = 3, line = 1, cex=1.5)
	h <- hist(values, breaks=breaks, plot=F)
	h$lcounts <- log(h$counts)
	for(p in c(1:length(h$counts))){
		polygon(c(breaks[p],breaks[p],breaks[p+1],breaks[p+1]), c(0,h$lcounts[p],h$lcounts[p],0), col="darkred")
	}
	axis(1, at = seq(0, 200, 20), lwd.ticks=1, las=1, cex.axis=1)
	axis(2, at = log(c(0,1,10,100,1000,10000)), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
	box()
}



BlanSharedGenes <- function(Data, main){
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,12), xlim=c(0, 17), col=NA)
	mtext("Number of orthologous groups", side = 2, line = 4, cex=1)
	mtext(main, side = 3, line = 1, cex=1.5)
	Values <- c(length(Data[which(Data[,"Blan"]>0),1]))
	ComparingNames <- c("Total")
	for(s in c(1:length(Amphi))){
		if(Amphi[s] != "Blan"){
			ComparingNames <- c(ComparingNames, paste("Shared with\n", Amphi[s]))
			Values <- c(Values, length(Data[which(rowSums(Data[,c(Amphi[s],"Blan")])==Data$Sum),1]))
		}
	}
	ComparingNames <- c(ComparingNames, paste("Shared with\n", paste(Amphi[-which(Amphi=="Blan")], collapse="\n"), sep=""))
	Values <- c(Values, length(Data[which(rowSums(Data[,Amphi])==Data$Sum & apply(Data[,Amphi], 1, function(x){length(grep("0", x))})==0),1]))
	ComparingNames <- c(ComparingNames, "Amphioxus\nspecific")
	Values <- c(Values, length(Data[which(rowSums(Data[,Amphi])==Data$Sum & Data$Blan>0),1]))

	for(s in c(1:length(Verteb))){
		ComparingNames <- c(ComparingNames, paste("Shared with\n", Verteb[s]))
		Values <- c(Values, length(Data[which(rowSums(Data[,c(Verteb[s],"Blan")])==Data$Sum),1]))
	}

	ComparingNames <- c(ComparingNames, paste("Shared with\n", paste(c(Amphi[-which(Amphi=="Blan")], Verteb), collapse="\n"), sep=""))
	Values <- c(Values, length(Data[which(rowSums(Data[,c(Amphi, Verteb)])==Data$Sum & apply(Data[,c(Amphi, Verteb)], 1, function(x){length(grep("0", x))})==0),1]))
	ComparingNames <- c(ComparingNames, "Shared with\nvertebrate\nspecies")
	Values <- c(Values, length(Data[which(rowSums(Data[,c(Amphi, Verteb)])==Data$Sum & Data$Blan>0 & rowSums(Data[,c(Verteb)])>0) ,1]))

	for(s in c(1:length(OutDeut))){
		ComparingNames <- c(ComparingNames, paste("Shared with\n", OutDeut[s]))
		Values <- c(Values, length(Data[which(rowSums(Data[,c(OutDeut[s],"Blan")])==Data$Sum),1]))
	}
	ComparingNames <- c(ComparingNames, "Shared with\nall species")
	Values <- c(Values, length(Data[which(apply(Data[,c(Amphi, Verteb, OutDeut)], 1, function(x){length(grep("0", x))})==0),1]))
	ComparingNames <- c(ComparingNames, "Shared with\nvertebrate &\nnon-chord\nspecies")
	Values <- c(Values, length(Data[which(Data$Blan>0 & rowSums(Data[,c(Verteb)])>0 & rowSums(Data[,c(OutDeut)])>0),1]))

	width <- .5
	for(v in c(1:length(Values))){
		polygon(c(v-width/2, v-width/2, v+width/2, v+width/2), c(0,log(Values[v]),log(Values[v]),0), col="gold")
		text(v, log(Values[v]), labels=Values[v], pos=3)
	}
	axis(1, at = c(1:length(Values)), labels=ComparingNames, lwd.ticks=1, padj=1, las=1, cex.axis=.8)
	axis(2, at = log(c(0,1,10,100,1000,10000)), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
	box()
}










