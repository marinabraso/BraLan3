


GeneFunction <- function(Gene.df, OG.df, GOMF){
	OG.blan <- OG.df[which(OG.df$BlanType!="Missing"),]

	pdf(paste(ResultsFolder, "/GeneFunction.pdf", sep=""), width=20, height=10)
	par(mar=c(5,5,3,3),oma=c(1,1,1,1), yaxs='i', xaxs='i')

	layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(3), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,103), xlim=c(0.5,length(MFunctGO)+.5), col=NA)
	mtext("% of duplicated genes", side = 2, line = 7, cex=2)
	Tdprop <- length(Gene.df[which(Gene.df$BlanType.AV=="Duplicated"),1])/length(Gene.df[,1])*100
	abline(h=Tdprop, lty=2, lwd=2, col="darkred")
	for(l in c(1:length(MFunctGO))){
		polygon(c(l-.5, l+.5, l+.5, l-.5),c(-5,-5,0,0), col=GeneralMFColors[GeneralMFNum[l]], border=NA, xpd = NA)
		tmp.Gene.df <- Gene.df[which(Gene.df$Gene %in% unique(GOMF$Gene[GOMF$GO %in% MFunctGO[l]])),]
		adprop <- length(tmp.Gene.df[which(tmp.Gene.df$BlanType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		vdprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		voprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Ohnolog"),1])/length(tmp.Gene.df[,1])*100
		points(c(l,l,l,l), c(adprop, vdprop, voprop, vdprop+voprop), col=c("black", VTypeColors[c(4,3)], "darkred"), bg=c("black", VTypeColors[c(4,3)], "darkred"), pch=16, cex=2)
	}
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)
	text(x = GeneralMFMidp, y = par("usr")[3]-10, labels = GeneralMFNames, xpd = NA, srt = 25, adj=1, cex = 1.2)


	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,103), xlim=c(0.5,length(MFunctGO)+.5), col=NA)
	mtext("% of duplicated orthologous groups", side = 2, line = 7, cex=2)
	Tdprop <- length(OG.blan[which(OG.df$BlanType=="Duplicated"),1])/length(OG.blan[,1])*100
	abline(h=Tdprop, lty=2, lwd=2, col="darkred")
	for(l in c(1:length(MFunctGO))){
		polygon(c(l-.5, l+.5, l+.5, l-.5),c(-5,-5,0,0), col=GeneralMFColors[GeneralMFNum[l]], border=NA, xpd = NA)
		tmp.Gene.df <- Gene.df[which(Gene.df$Gene %in% unique(GOMF$Gene[GOMF$GO %in% MFunctGO[l]])),]
		tmp.OG.blan <- OG.blan[which(rownames(OG.blan) %in% tmp.Gene.df$OG.AV),]
		adprop <- length(tmp.OG.blan[which(tmp.OG.blan$BlanType=="Duplicated"),1])/length(tmp.OG.blan[,1])*100
		vdprop <- length(tmp.OG.blan[which(tmp.OG.blan$VertebType=="Duplicated"),1])/length(tmp.OG.blan[,1])*100
		voprop <- length(tmp.OG.blan[which(tmp.OG.blan$VertebType=="Ohnolog"),1])/length(tmp.OG.blan[,1])*100
		points(c(l,l,l,l), c(adprop, vdprop, voprop, vdprop+voprop), col=c("black", VTypeColors[c(4,3)], "darkred"), bg=c("black", VTypeColors[c(4,3)], "darkred"), pch=16, cex=2)
	}
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)
	text(x = GeneralMFMidp, y = par("usr")[3]-10, labels = GeneralMFNames, xpd = NA, srt = 25, adj=1, cex = 1.2)

	layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,103), xlim=c(0,103), col=NA)
	mtext("% of duplicated genes", side = 2, line = 7, cex=2)
	Tdprop <- length(Gene.df[which(Gene.df$BlanType.AV=="Duplicated"),1])/length(Gene.df[,1])*100
	abline(h=Tdprop, lty=2, lwd=2, col="darkred")
	for(l in c(1:length(MFunctGO))){
		tmp.Gene.df <- Gene.df[which(Gene.df$Gene %in% unique(GOMF$Gene[GOMF$GO %in% MFunctGO[l]])),]
		adprop <- length(tmp.Gene.df[which(tmp.Gene.df$BlanType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		vdprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		voprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Ohnolog"),1])/length(tmp.Gene.df[,1])*100
		points(vdprop, adprop, col=GeneralMFColors[GeneralMFNum[l]], bg=GeneralMFColors[GeneralMFNum[l]], pch=16, cex=2)
	}
	axis(1, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,103), xlim=c(0,103), col=NA)
	mtext("% of duplicated genes", side = 2, line = 7, cex=2)
	Tdprop <- length(Gene.df[which(Gene.df$BlanType.AV=="Duplicated"),1])/length(Gene.df[,1])*100
	abline(h=Tdprop, lty=2, lwd=2, col="darkred")
	for(l in c(1:length(MFunctGO))){
		tmp.Gene.df <- Gene.df[which(Gene.df$Gene %in% unique(GOMF$Gene[GOMF$GO %in% MFunctGO[l]])),]
		adprop <- length(tmp.Gene.df[which(tmp.Gene.df$BlanType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		vdprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		voprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Ohnolog"),1])/length(tmp.Gene.df[,1])*100
		points(voprop, adprop, col=GeneralMFColors[GeneralMFNum[l]], bg=GeneralMFColors[GeneralMFNum[l]], pch=16, cex=2)
	}
	axis(1, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,103), xlim=c(0,103), col=NA)
	mtext("% of duplicated genes", side = 2, line = 7, cex=2)
	Tdprop <- length(Gene.df[which(Gene.df$BlanType.AV=="Duplicated"),1])/length(Gene.df[,1])*100
	abline(h=Tdprop, lty=2, lwd=2, col="darkred")
	for(l in c(1:length(MFunctGO))){
		tmp.Gene.df <- Gene.df[which(Gene.df$Gene %in% unique(GOMF$Gene[GOMF$GO %in% MFunctGO[l]])),]
		adprop <- length(tmp.Gene.df[which(tmp.Gene.df$BlanType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		vdprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Duplicated"),1])/length(tmp.Gene.df[,1])*100
		voprop <- length(tmp.Gene.df[which(tmp.Gene.df$VertebType.AV=="Ohnolog"),1])/length(tmp.Gene.df[,1])*100
		points(vdprop+voprop, adprop, col=GeneralMFColors[GeneralMFNum[l]], bg=GeneralMFColors[GeneralMFNum[l]], pch=16, cex=2)
	}
	axis(1, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)


	plot.new()
	legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VTypeColors), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)


	dev.off()
}