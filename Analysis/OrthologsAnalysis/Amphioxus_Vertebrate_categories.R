


Amphioxus_Vertebrate_categories <- function(OG.df){
	# Print numbers
	Shared <- length(OG.df[which(OG.df$AmphiType != "Missing" & OG.df$VertebType != "Missing"),1])
	VertebTotal <- length(OG.df[which(OG.df$VertebType != "Missing"),1])
	AmphiTotal <- length(OG.df[which(OG.df$AmphiType != "Missing"),1])
	Control <- length(OG.df[which(OG.df$AmphiType == "Missing" & OG.df$VertebType == "Missing"),1])
	print("Shared VertebTotal AmphiTotal Control")
	print(paste(Shared, VertebTotal, AmphiTotal, Control))
	print("Shared/VertebTotal*100")
	print(Shared/VertebTotal*100)
	print("Shared/AmphiTotal*100")
	print(Shared/AmphiTotal*100)

	for(sp in c(Amphi, Verteb)){
		total <- sum(OG.df[,sp]>0)
		dup <- sum(OG.df[,sp]>1)
		totalG <- sum(OG.df[which(OG.df[,sp]>0),sp])
		dupG <- sum(OG.df[which(OG.df[,sp]>1),sp])
		print(paste(sp, total, dup, dup/total*100, totalG, dupG, dupG/totalG*100))
	}

	# Testing significance coexistence of amphioxus and vertebrate categories
	OG.blan <- OG.df[which(OG.df$BlanType!="Missing"),]

	HyperTests <- data.frame()
	for(vtype in c("Missing","SingleCopy","Ohnolog","Duplicated")){
		for(atype in c("SingleCopy","Duplicated")){
			vnum <- length(OG.blan[which(OG.blan$VertebType == vtype),1])
			anum <- length(OG.blan[which(OG.blan$BlanType == atype),1])
			onum <- length(OG.blan[which(OG.blan$VertebType == vtype & OG.blan$BlanType == atype),1])
			HyperTests <- rbind(HyperTests, unlist(HypergeometricTest(onum, anum, vnum, length(OG.blan[,1]), atype, vtype, PQvalThreshold)))
		}
	}
	colnames(HyperTests) <- c("VertebType", "BlanType", "FoldChange", "Observed", "Expected", "HpvalDepleted", "HpvalEnriched")
	HyperTests$HqvalDepleted <- qvalue(as.numeric(HyperTests$HpvalDepleted))$qvalues
	HyperTests$HqvalEnriched <- qvalue(as.numeric(HyperTests$HpvalEnriched))$qvalues
	HyperTests$HResult <- rep("NA", length(HyperTests[,1]))
	HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)] <- rep("D", length(HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)]))
	HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)] <- rep("E", length(HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)]))
	print(HyperTests)

	# Plotting
	pdf(paste(ResultsFolder, "/Amphioxus_Vertebrate_categories.pdf", sep=""), width=20, height=10)
	par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
	## Bar plot of Blan/Verteb categories
	layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 2.5), col=NA)
	mtext("% in each category", side = 2, line = 5, cex=2)
	PlotColumnVertebrateType(OG.blan[which(OG.blan$Blan==1),], 1, VTypeColors, 2)
	PlotColumnVertebrateType(OG.blan[which(OG.blan$Blan>=2),], 2, VTypeColors, 2, den=10)
	#abline(h=length(OG.blan[which(OG.blan$VertebType=="SingleCopy"),1])/length(OG.blan[,1])*100, lty=2, lwd=2, col="darkred")
	axis(1, at = c(1:2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2, las=1, cex.axis=1.5)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot.new()
	legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VTypeColors), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
	plot.new()
	plot.new()

	## 2nd Bar plot of Blan/Verteb categories
	layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,40), xlim=c(0.5, 4.5), col=NA)
	mtext("% in each category", side = 2, line = 5, cex=2)
	w <- .3
	for(vt in c(1:length(VTypes))){
		scvalue <- length(OG.blan[which(OG.blan$Blan==1 & OG.blan$VertebType==VTypes[vt]),1])/length(OG.blan[which(OG.blan$Blan==1),1])*100
		dvalue <- length(OG.blan[which(OG.blan$Blan>=2 & OG.blan$VertebType==VTypes[vt]),1])/length(OG.blan[which(OG.blan$Blan>=2),1])*100
		polygon(c(vt-w, vt, vt, vt-w),c(0,0,scvalue, scvalue), col=VTypeColors[vt], border=VTypeColors[vt], lwd=4)		
		polygon(c(vt, vt+w, vt+w, vt),c(0,0,dvalue, dvalue), col=modifColor(VTypeColors[vt], .1), border=VTypeColors[vt], lwd=4)		
		polygon(c(vt, vt+w, vt+w, vt),c(0,0,dvalue, dvalue), col=VTypeColors[vt], border=VTypeColors[vt], lwd=4, den=10)		
	}
	axis(2, at = seq(0,100,5), lwd.ticks=1, las=1, cex.axis=2)
	axis(1, at = seq(1,4,1), labels=VTypes, lwd.ticks=1, las=1, cex.axis=1.5)


	dev.off()
}




