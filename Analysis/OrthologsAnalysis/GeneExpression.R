
### Gene expression and specificity boxplots for Blan categories


GeneExpression <- function(Gene.df){
	pdf(paste(ResultsFolder, "/GeneExpression.pdf", sep=""), width=20, height=10)
	par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
	layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
	mtext("Mean adult TPM", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$BlanType.AV=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$BlanType.AV=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
	axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
	mtext("Mean embrionic TPM", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$BlanType.AV=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$BlanType.AV=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
	axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
	mtext("Tau among tissues", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$BlanType.AV=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$BlanType.AV=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
	axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
	mtext("Tau among developmental stages", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$BlanType.AV=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$BlanType.AV=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
	axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)


	### Gene specificity boxplots for tissues & developmental stages
	layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(Tissues)*2+20), col=NA)
	mtext("Tau", side = 2, line = 5, cex=2)
	for(t in c(1:length(Tissues))){
		BoxPlot(Gene.df$TauTissues[which(Gene.df$BlanType.AV=="SingleCopy" & Gene.df$MaxTissue==Tissues[t])], t, TissueColors[t], 1)
		BoxPlot(Gene.df$TauTissues[which(Gene.df$BlanType.AV=="Duplicated" & Gene.df$MaxTissue==Tissues[t])], t+length(Tissues)+3, TissueColors[t], 1)
	}
	axis(1, at = c(1+length(Tissues)/2, 1+length(Tissues)+3+length(Tissues)/2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
	legend("bottomright", Tissues, pch=19, text.col="black", col=TissueColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(EmbAges)*2+20), col=NA)
	mtext("Tau", side = 2, line = 5, cex=2)
	for(t in c(1:length(EmbAges))){
		BoxPlot(Gene.df$TauEmbAge[which(Gene.df$BlanType.AV=="SingleCopy" & Gene.df$MaxEmbAge==EmbAges[t])], t, EmbAgesColors[t], .7)
		BoxPlot(Gene.df$TauEmbAge[which(Gene.df$BlanType.AV=="Duplicated" & Gene.df$MaxEmbAge==EmbAges[t])], t+length(EmbAges)+3, EmbAgesColors[t], .7)
	}
	axis(1, at = c(1+length(EmbAges)/2, 1+length(EmbAges)+3+length(EmbAges)/2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
	legend("bottomright", EmbAges, pch=19, text.col="black", col=EmbAgesColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)

	### Gene expression & specificity boxplots for vertebrate & blan categories
	dist <- .15
	layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 4.5), col=NA)
	mtext("Mean adult TPM", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
	BoxPlot(Gene.df$MeanAdult[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
	axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 4.5), col=NA)
	mtext("Mean embrionic TPM", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
	BoxPlot(Gene.df$MeanEmbr[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
	axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 4.5), col=NA)
	mtext("Tau among tissues", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
	BoxPlot(Gene.df$TauTissues[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
	axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 4.5), col=NA)
	mtext("Tau among developmental stages", side = 2, line = 5, cex=2)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="Missing" & Gene.df$BlanType.AV=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="SingleCopy" & Gene.df$BlanType.AV=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="Ohnolog" & Gene.df$BlanType.AV=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
	BoxPlot(Gene.df$TauEmbAge[which(Gene.df$VertebType.AV=="Duplicated" & Gene.df$BlanType.AV=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
	axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
	axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

	### Expression per vertebrate, blan categories and for developmental stages and tissues
	layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
	plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 4), xlim=c(0.8, length(EmbAges)+length(Tissues)+.6), col=NA)
	mtext("Median log(TPM+1)", side = 2, line = 5, cex=2)
	for(v in c(1:length(VTypes))){
		for(a in c(1:length(BlanTypes))){
			vec <- c()
			for(age in c(1:length(EmbAges))){
				vec <- c(vec, median(log(Gene.df[which(Gene.df$VertebType.AV==VTypes[v] & Gene.df$BlanType.AV==BlanTypes[a]), EmbAges[age]]+1)))
			}
			lines(c(1:length(EmbAges)), vec, col=VTypeColors[v], lty=BlanType.AV.lty[a], lwd=3)
			vec <- c()
			for(t in c(1:length(Tissues))){
				vec <- c(vec, median(log(Gene.df[which(Gene.df$VertebType.AV==VTypes[v] & Gene.df$BlanType.AV==BlanTypes[a]), Tissues[t]]+1)))			
				#points(jitter(length(EmbAges)+t, amount=.2), vec[-1], pch=BlanType.AV.pch[a], cex=2, col=VTypeColors[v])
			}
			lines(length(EmbAges)+c(1:length(Tissues)), vec, col=VTypeColors[v], lty=BlanType.AV.lty[a], lwd=2)
			points(length(EmbAges)+c(1:length(Tissues)), vec, pch=BlanType.AV.pch[a], cex=2, col=VTypeColors[v])
		}
	}
	for(t in c(1:length(Tissues))){
		abline(v=length(EmbAges)+t-.5, col="grey90", lty=2)
	}
	abline(v=length(EmbAges)+t+.5, col="grey90", lty=2)
	axis(1, at = c((length(EmbAges)+1):(length(EmbAges)+length(Tissues))) , labels=Tissues, lwd.ticks=1, line=1, las=1, cex.axis=1)
	axis(1, at = c(1:length(EmbAges)), labels=rep("",length(EmbAges)), lwd.ticks=1, line=1, las=1, cex.axis=1)
	axis(2, at = seq(0,10,1), lwd.ticks=1, las=1, cex.axis=2)


	dev.off()
}
