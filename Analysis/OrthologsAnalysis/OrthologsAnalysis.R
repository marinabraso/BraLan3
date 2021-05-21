#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
library(qvalue)
library(ghibli)

######################################################################
# Files & folders

ResultsFolder <- "Plots/OrthologsAnalysis"
ExpressionDataFile <- paste(ResultsFolder, "/ExpressionDataProcessed.txt", sep ="")
OGCountsFile <- paste(ResultsFolder, "/OGCountsProcessed.txt", sep ="")
RNASeqMetadataFile <- paste(ResultsFolder, "/RNASeqMetadataProcessed.txt", sep ="")

######################################################################
# General parameters
QvalT <- 0.1 # q-value threshold for positive selection
PQvalThreshold <- 0.01 # for hypergeometric test

Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Blan", "Bflo", "Bbel")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Amphi, Verteb, OutDeut)

Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(Tissues))

EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre.metamorphic.larvae")
#EmbAges <- c("32cells.Blastula", "7h.8h", "10h.11h", "15h.18h", "21h.24h", "27h.36h", "50h.60h")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
EmbAgesColors <- colfunc(length(EmbAges))

VTypes <- c("Missing", "SingleCopy", "Ohnolog", "Duplicated")
BlanTypes <- c("SingleCopy", "Duplicated")
VTypeColors <- c("deepskyblue2", "gold", "olivedrab2", "olivedrab4")
VTypeColors <- c("#481567FF", "#FDE725FF", "#1F968BFF", "#39568CFF")
VTypeColors <- ghibli_palettes$MarnieMedium2[c(7,2,4,6)]
BaseColor <- ghibli_palettes$MarnieMedium2[5]
BlanType.pch <- c(16,18)
BlanType.lty <- c(1,2)

###########################################################################
###########################################################################
# Read data

# Expression data
ExpData <- read.table(ExpressionDataFile, h=T, sep = "\t", row.names=1)
colnames(ExpData) <- sub("^X", "", colnames(ExpData))
ExpData[,"32cells.Blastula"] <- rowMeans(ExpData[,c("32cells", "Blastula")])
ExpData[,"7h.8h"] <- rowMeans(ExpData[,c("7h","8h")])
ExpData[,"10h.11h"] <- rowMeans(ExpData[,c("10h","11h")])
ExpData[,"15h.18h"] <- rowMeans(ExpData[,c("15h", "18h")])
ExpData[,"21h.24h"] <- rowMeans(ExpData[,c("21h", "24h")])
ExpData[,"27h.36h"] <- rowMeans(ExpData[,c("27h", "36h")])
ExpData[,"50h.60h"] <- rowMeans(ExpData[,c("50h", "60h")])
print(head(ExpData))

OGData.AV <- read.table(OGCountsFile, h=T, sep = "\t", row.names=1)
print(head(OGData.AV))

MetaRNA <- read.table(RNASeqMetadataFile, h=T, sep = "\t", row.names=1)
print(head(MetaRNA))


###########################################################################
###########################################################################

Shared <- length(OGData.AV[which(OGData.AV$AmphiType != "Missing" & OGData.AV$VertebType != "Missing"),1])
VertebTotal <- length(OGData.AV[which(OGData.AV$VertebType != "Missing"),1])
AmphiTotal <- length(OGData.AV[which(OGData.AV$AmphiType != "Missing"),1])
Control <- length(OGData.AV[which(OGData.AV$AmphiType == "Missing" & OGData.AV$VertebType == "Missing"),1])
print(paste(Shared, VertebTotal, AmphiTotal, Control))
print(Shared/VertebTotal*100)
print(Shared/AmphiTotal*100)

# Testing
OGData.Blan <- OGData.AV[which(OGData.AV$BlanType!="Missing"),]

HyperTests <- data.frame()
for(vtype in c("Missing","SingleCopy","Ohnolog","Duplicated")){
	for(atype in c("SingleCopy","Duplicated")){
		vnum <- length(OGData.Blan[which(OGData.Blan$VertebType == vtype),1])
		anum <- length(OGData.Blan[which(OGData.Blan$BlanType == atype),1])
		onum <- length(OGData.Blan[which(OGData.Blan$VertebType == vtype & OGData.Blan$BlanType == atype),1])
		HyperTests <- rbind(HyperTests, unlist(HypergeometricTest(onum, anum, vnum, length(OGData.Blan[,1]), atype, vtype, PQvalThreshold)))
	}
}
colnames(HyperTests) <- c("VertebType", "BlanType", "FoldChange", "Observed", "Expected", "HpvalDepleted", "HpvalEnriched")
HyperTests$HqvalDepleted <- qvalue(as.numeric(HyperTests$HpvalDepleted))$qvalues
HyperTests$HqvalEnriched <- qvalue(as.numeric(HyperTests$HpvalEnriched))$qvalues
HyperTests$HResult <- rep("NA", length(HyperTests[,1]))
HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)] <- rep("D", length(HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)]))
HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)] <- rep("E", length(HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)]))
print(HyperTests)


#RTests.PosSel <- data.frame()
#for(vtype in c("Missing","SingleCopy","Ohnolog","Duplicated")){
#	print(vtype)
#	Typenum <- length(OGData.Blan[which(OGData.Blan$VertebType==vtype),1])
#	Dupnum <- length(OGData.Blan[which(OGData.Blan$VertebType==vtype & OGData.Blan$BlanType=="Duplicated"),1])
#	PSelnum <- length(OGData.Blan[which(OGData.Blan$VertebType==vtype & OGData.Blan$M8.Qval<=QvalT),1])
#	Overlap <- length(OGData.Blan[which(OGData.Blan$VertebType==vtype & OGData.Blan$BlanType=="Duplicated" & OGData.Blan$M8.Qval<=QvalT),1])
#	RTests.PosSel <- rbind(RTests.PosSel, c(vtype, unlist(RandomizationTestProportion(Overlap, Dupnum, PSelnum, Typenum, "Duplicated", "Under positive selection", 10000))))
#}
#colnames(RTests.PosSel) <- c("Univers", "Group1", "Group2", "pvalDepleted", "pvalEnriched", "MeanExpected", "SDExpected")
#print(RTests.PosSel)

for(sp in c(Amphi, Verteb)){
	total <- sum(OGData.Blan[,sp]>0)
	dup <- sum(OGData.Blan[,sp]>1)
	totalG <- sum(OGData.Blan[which(OGData.Blan[,sp]>0),sp])
	dupG <- sum(OGData.Blan[which(OGData.Blan[,sp]>1),sp])
	print(paste(sp, total, dup, dup/total*100, totalG, dupG, dupG/totalG*100))
}

quit()





###########################################################################
###########################################################################
# Plotting
###########################################################################
###########################################################################
pdf(paste(ResultsFolder, "/OrthologsAnalysis.pdf", sep=""), width=20, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')


## Bar plot of Blan/Verteb categories
layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 2.5), col=NA)
mtext("% in each category", side = 2, line = 5, cex=2)
PlotColumnVertebrateType(OGData.AV[which(OGData.AV$Blan==1),], 1, VTypeColors, 2)
PlotColumnVertebrateType(OGData.AV[which(OGData.AV$Blan>=2),], 2, VTypeColors, 2, den=10)
axis(1, at = c(1:2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2, las=1, cex.axis=1.5)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot.new()
legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VTypeColors), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
plot.new()
plot.new()


### Gene expression and specificity boxplots for Blan categories
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
mtext("Mean adult TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanAdult[which(ExpData$BlanType=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
BoxPlot(ExpData$MeanAdult[which(ExpData$BlanType=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
mtext("Mean embrionic TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanEmbr[which(ExpData$BlanType=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
BoxPlot(ExpData$MeanEmbr[which(ExpData$BlanType=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
mtext("Tau among tissues", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
mtext("Tau among developmental stages", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauEmbAge[which(ExpData$BlanType=="SingleCopy")], 1, BaseColor, 1.5, w=.6)
BoxPlot(ExpData$TauEmbAge[which(ExpData$BlanType=="Duplicated")], 2, BaseColor, 1.5, w=.6, den=10)
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)


### Gene specificity boxplots for tissues & developmental stages
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(Tissues)*2+20), col=NA)
mtext("Tau", side = 2, line = 5, cex=2)
for(t in c(1:length(Tissues))){
	BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="SingleCopy" & ExpData$MaxTissue==Tissues[t])], t, TissueColors[t], 1)
	BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="Duplicated" & ExpData$MaxTissue==Tissues[t])], t+length(Tissues)+3, TissueColors[t], 1)
}
axis(1, at = c(1+length(Tissues)/2, 1+length(Tissues)+3+length(Tissues)/2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
legend("bottomright", Tissues, pch=19, text.col="black", col=TissueColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(EmbAges)*2+20), col=NA)
mtext("Tau", side = 2, line = 5, cex=2)
for(t in c(1:length(EmbAges))){
	BoxPlot(ExpData$TauEmbAge[which(ExpData$BlanType=="SingleCopy" & ExpData$MaxEmbAge==EmbAges[t])], t, EmbAgesColors[t], .7)
	BoxPlot(ExpData$TauEmbAge[which(ExpData$BlanType=="Duplicated" & ExpData$MaxEmbAge==EmbAges[t])], t+length(EmbAges)+3, EmbAgesColors[t], .7)
}
axis(1, at = c(1+length(EmbAges)/2, 1+length(EmbAges)+3+length(EmbAges)/2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
legend("bottomright", EmbAges, pch=19, text.col="black", col=EmbAgesColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)

### Gene expression & specificity boxplots for vertebrate & blan categories
dist <- .15
layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 4.5), col=NA)
mtext("Mean adult TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 4.5), col=NA)
mtext("Mean embrionic TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 4.5), col=NA)
mtext("Tau among tissues", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 4.5), col=NA)
mtext("Tau among developmental stages", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="SingleCopy")], 1-dist, VTypeColors[1], .7, .4)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Missing" & ExpData$BlanType=="Duplicated")], 1+dist, VTypeColors[1], .7, .4, 10)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="SingleCopy")], 2-dist, VTypeColors[2], .7, .4)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="SingleCopy" & ExpData$BlanType=="Duplicated")], 2+dist, VTypeColors[2], .7, .4, 10)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="SingleCopy")], 3-dist, VTypeColors[3], .7, .4)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Ohnolog" & ExpData$BlanType=="Duplicated")], 3+dist, VTypeColors[3], .7, .4, 10)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="SingleCopy")], 4-dist, VTypeColors[4], .7, .4)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Duplicated" & ExpData$BlanType=="Duplicated")], 4+dist, VTypeColors[4], .7, .4, 10)
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
			vec <- c(vec, median(log(ExpData[which(ExpData$VertebType==VTypes[v] & ExpData$BlanType==BlanTypes[a]), EmbAges[age]]+1)))
		}
		lines(c(1:length(EmbAges)), vec, col=VTypeColors[v], lty=BlanType.lty[a], lwd=3)
		vec <- c()
		for(t in c(1:length(Tissues))){
			vec <- c(vec, median(log(ExpData[which(ExpData$VertebType==VTypes[v] & ExpData$BlanType==BlanTypes[a]), Tissues[t]]+1)))			
			#points(jitter(length(EmbAges)+t, amount=.2), vec[-1], pch=BlanType.pch[a], cex=2, col=VTypeColors[v])
		}
		lines(length(EmbAges)+c(1:length(Tissues)), vec, col=VTypeColors[v], lty=BlanType.lty[a], lwd=2)
		points(length(EmbAges)+c(1:length(Tissues)), vec, pch=BlanType.pch[a], cex=2, col=VTypeColors[v])
	}
}
for(t in c(1:length(Tissues))){
	abline(v=length(EmbAges)+t-.5, col="grey90", lty=2)
}
abline(v=length(EmbAges)+t+.5, col="grey90", lty=2)
axis(1, at = c((length(EmbAges)+1):(length(EmbAges)+length(Tissues))) , labels=Tissues, lwd.ticks=1, line=1, las=1, cex.axis=1)
axis(1, at = c(1:length(EmbAges)), labels=rep("",length(EmbAges)), lwd.ticks=1, line=1, las=1, cex.axis=1)
axis(2, at = seq(0,10,1), lwd.ticks=1, las=1, cex.axis=2)


layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
mtext("M8 D", side = 2, line = 5, cex=2)
PlotaJitterPoints(ExpData$M8.D[which(ExpData$BlanType=="SingleCopy")], 1, VTypeColors[2])
PlotaJitterPoints(ExpData$M8.D[which(ExpData$BlanType=="Duplicated")], 2, VTypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
mtext("M8 q-value", side = 2, line = 5, cex=2)
PlotaJitterPoints(ExpData$M8.Qval[which(ExpData$BlanType=="SingleCopy")], 1, VTypeColors[2])
PlotaJitterPoints(ExpData$M8.Qval[which(ExpData$BlanType=="Duplicated")], 2, VTypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=4, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)


layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
xlim <- 0.1
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-2, length(VTypes)+1), xlim=c(-xlim, xlim), col=NA)
mtext("Per gene", side = 3, line = 2, cex=2)
ypos <- 0.5
PlotRandomizeDifferenceDupSC(ExpData, ypos, BaseColor, 5)
for(v in c(1:length(VTypes))){
	ypos <- ypos+1
	PlotRandomizeDifferenceDupSC(ExpData[which(ExpData$VertebType==VTypes[v]),], ypos, VTypeColors[v], 5)
}
abline(v=0)
abline(h=0)
abline(h=1, lty=2, col="grey70")
arrows(x0=xlim/4*3, y0=-1, x1=xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
text(xlim/6, -1.5, pos=4, labels="+ in duplicated genes", cex=1.5)
arrows(x0=-xlim/4*3, y0=-1, x1=-xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
text(-xlim/6, -1.5, pos=2, labels="+ in single copy genes", cex=1.5)
text(rep(-xlim,4), c(0:length(VTypes))+.5, labels=c("Total", VTypes), pos=4, cex=1.5)

xlim <- 0.1
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-2, length(VTypes)+1), xlim=c(-xlim, xlim), col=NA)
mtext("Per orthologous group", side = 3, line = 2, cex=2)
ypos <- 0.5
PlotRandomizeDifferenceDupSC(OGData.Blan, ypos, BaseColor, 5)
for(v in c(1:length(VTypes))){
	ypos <- ypos+1
	PlotRandomizeDifferenceDupSC(OGData.Blan[which(OGData.Blan$VertebType==VTypes[v]),], ypos, VTypeColors[v], 5)
}
abline(v=0)
abline(h=0)
abline(h=1, lty=2, col="grey70")
arrows(x0=xlim/4*3, y0=-1, x1=xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
text(xlim/6, -1.5, pos=4, labels="+ in duplicated genes", cex=1.5)
arrows(x0=-xlim/4*3, y0=-1, x1=-xlim/4, y1=-1, code=1, length=0.1, col="black", lwd=3)
text(-xlim/6, -1.5, pos=2, labels="+ in single copy genes", cex=1.5)
text(rep(-xlim,4), c(0:length(VTypes))+.5, labels=c("Total", VTypes), pos=4, cex=1.5)
quit()

xlim <- 0.1
step <- 5
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-4, 40), xlim=c(-xlim, xlim), col=NA)
mtext("Mean adult expression", side = 2, line = 5, cex=2)
for(i in seq(0, 40, step)){
	tmpExpData <- ExpData[which(ExpData$MeanAdult>=i & ExpData$MeanAdult<(i+step)),]
	PlotRandomizeDifferenceDupSC(tmpExpData, i+step/2, BaseColor, 5)
}
abline(v=0)
abline(h=0)
arrows(x0=xlim/4*3, y0=-1.5, x1=xlim/4, y1=-1.5, code=1, length=0.1, col="black", lwd=3)
text(xlim/6, -3, pos=4, labels="+ in duplicated genes", cex=1.5)
arrows(x0=-xlim/4*3, y0=-1.5, x1=-xlim/4, y1=-1.5, code=1, length=0.1, col="black", lwd=3)
text(-xlim/6, -3, pos=2, labels="+ in single copy genes", cex=1.5)
axis(2, at = seq(0,40,step), lwd.ticks=1, las=1, cex.axis=2)


xlim <- 0.1
step <- 5
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(-2, 6), xlim=c(-xlim, xlim), col=NA)
mtext("Blan copy number", side = 2, line = 5, cex=2)
for(i in c(2:4)){
	tmpOGData.Blan <- OGData.Blan[which(OGData.Blan$Blan==i | OGData.Blan$BlanType=="SingleCopy"),]
	PlotRandomizeDifferenceDupSC(tmpOGData.Blan, i, BaseColor, 5)
}
tmpOGData.Blan <- OGData.Blan[which(OGData.Blan$Blan>=5 | OGData.Blan$BlanType=="SingleCopy"),]
PlotRandomizeDifferenceDupSC(tmpOGData.Blan, 5, BaseColor, 5)
abline(v=0)
abline(h=1)
arrows(x0=xlim/4*3, y0=0, x1=xlim/4, y1=0, code=1, length=0.1, col="black", lwd=3)
text(xlim/6, -1, pos=4, labels="+ in duplicated genes", cex=1.5)
arrows(x0=-xlim/4*3, y0=0, x1=-xlim/4, y1=0, code=1, length=0.1, col="black", lwd=3)
text(-xlim/6, -1, pos=2, labels="+ in single copy genes", cex=1.5)
axis(2, at = c(2:5), labels=c(c(2:4), "5+"), lwd.ticks=1, line=1, las=1, cex.axis=2)






















dev.off()












