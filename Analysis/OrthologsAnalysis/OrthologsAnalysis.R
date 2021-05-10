#!/usr/bin/env Rscript



######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))

######################################################################
# Files & folders

ResultsFolder <- "Plots/OrthologsAnalysis"
ExpressionDataFile <- paste(ResultsFolder, "/ExpressionDataProcessed.txt", sep ="")
OGCountsFile <- paste(ResultsFolder, "/OGCountsProcessed.txt", sep ="")
RNASeqMetadataFile <- paste(ResultsFolder, "/RNASeqMetadataProcessed.txt", sep ="")

######################################################################
# General parameters

Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Blan", "Bflo", "Bbel")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Amphi, Verteb, OutDeut)

TypeColors <- c("deepskyblue2", "gold", "olivedrab2", "olivedrab4")

Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic diverticulum", "muscle", "neural tube", "female gonads", "male gonads")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(Tissues))

EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre-metamorphic larvae")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
EmbAgesColors <- colfunc(length(EmbAges))

###########################################################################
###########################################################################
# Read data

# Expression data
ExpData <- read.table(ExpressionDataFile, h=T, sep = "\t", row.names=1)
print(head(ExpData))

CountsAV <- read.table(OGCountsFile, h=T, sep = "\t", row.names=1)
print(head(CountsAV))

MetaRNA <- read.table(RNASeqMetadataFile, h=T, sep = "\t", row.names=1)
print(head(MetaRNA))


###########################################################################
###########################################################################
# Testing

CountBlan <- CountsAV[which(CountsAV$BlanType!="Missing"),]
# SC genes in Blan tend to be SC in Vertebrates?
chisq.test(CountBlan$BlanType == "SingleCopy", CountBlan$VertebType == "SingleCopy")
sum(CountBlan$BlanType == "SingleCopy")*sum(CountBlan$VertebType == "SingleCopy")/length(CountBlan[,1])
sum(CountBlan$BlanType == "SingleCopy" & CountBlan$VertebType == "SingleCopy")

# Dup genes in Blan tend to be Ohn in Vertebrates?
chisq.test(CountBlan$BlanType == "Duplicated", CountBlan$VertebType == "Ohnolog")
sum(CountBlan$BlanType == "Duplicated")*sum(CountBlan$VertebType == "Ohnolog")/length(CountBlan[,1])
sum(CountBlan$BlanType == "Duplicated" & CountBlan$VertebType == "Ohnolog")

# Dup genes in Blan tend to be Dup in Vertebrates?
chisq.test(CountBlan$BlanType == "Duplicated", CountBlan$VertebType == "Duplicated")
sum(CountBlan$BlanType == "Duplicated")*sum(CountBlan$VertebType == "Duplicated")/length(CountBlan[,1])
sum(CountBlan$BlanType == "Duplicated" & CountBlan$VertebType == "Duplicated")

# Do dup Blan genes have more Onh orthologs than expected? 
CountBlanPVertebD <- CountBlan[which(CountBlan$VertebType!="Missing" & CountBlan$VertebType!="SingleCopy"),]
chisq.test(CountBlanPVertebD$BlanType == "Duplicated", CountBlanPVertebD$VertebType == "Ohnolog")
sum(CountBlanPVertebD$BlanType == "Duplicated")*sum(CountBlanPVertebD$VertebType == "Ohnolog")/length(CountBlanPVertebD[,1])
sum(CountBlanPVertebD$BlanType == "Duplicated" & CountBlanPVertebD$VertebType == "Ohnolog")


###########################################################################
###########################################################################
# Plotting
###########################################################################
###########################################################################
pdf(paste(ResultsFolder, "/OrthologsAnalysis.pdf", sep=""), width=20, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 3.5), col=NA)
PlotColumnVertebrateType(CountsAV[which(CountsAV$Blan==1),], 1, TypeColors)
PlotColumnVertebrateType(CountsAV[which(CountsAV$Blan>=2),], 2, TypeColors)
axis(1, at = c(1:2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,19,19,19,19), text.col="black", col=c(NA, TypeColors), bty = "n", cex=2, xjust = 0, yjust = 0)

#plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 5.5), col=NA)
#PlotColumnBlanCN(CountsAV[which(CountsAV$VertebType=="Missing"),], 1, TypeColors[c(1,2,4)])
#PlotColumnBlanCN(CountsAV[which(CountsAV$VertebType=="SingleCopy"),], 2, TypeColors[c(1,2,4)])
#PlotColumnBlanCN(CountsAV[which(CountsAV$VertebType=="Ohnolog"),], 3, TypeColors[c(1,2,4)])
#PlotColumnBlanCN(CountsAV[which(CountsAV$VertebType=="Duplicated"),], 4, TypeColors[c(1,2,4)])
#axis(1, at = c(1:4), labels=c("Missing\nin vertebrates", "Single copy\nin vertebrates", "Ohnolog\nin vertebrates", "Duplicated\nin vertebrates"), tick=FALSE, line=2.5, las=1, cex.axis=2)
#axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
#legend("bottomright", c("In Blan", "Missing", "Single copy", "Duplicated"), pch=c(NA,19,19,19), text.col="black", col=c(NA, TypeColors[c(1,2,4)]), bty = "n", cex=2, xjust = 0, yjust = 0)




layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
mtext("Mean adult TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanAdult[which(ExpData$BlanType=="SingleCopy")], 1, TypeColors[2])
BoxPlot(ExpData$MeanAdult[which(ExpData$BlanType=="Duplicated")], 2, TypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 200), xlim=c(0.5, 2.5), col=NA)
mtext("Max adult TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MaxAdult[which(ExpData$BlanType=="SingleCopy")], 1, TypeColors[2])
BoxPlot(ExpData$MaxAdult[which(ExpData$BlanType=="Duplicated")], 2, TypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,200,40), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
mtext("Mean embrionic TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanEmbr[which(ExpData$BlanType=="SingleCopy")], 1, TypeColors[2])
BoxPlot(ExpData$MeanEmbr[which(ExpData$BlanType=="Duplicated")], 2, TypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 200), xlim=c(0.5, 2.5), col=NA)
mtext("Max embrionic TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MaxEmbr[which(ExpData$BlanType=="SingleCopy")], 1, TypeColors[2])
BoxPlot(ExpData$MaxEmbr[which(ExpData$BlanType=="Duplicated")], 2, TypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,200,40), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
mtext("Tau among tissues", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="SingleCopy")], 1, TypeColors[2])
BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="Duplicated")], 2, TypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 2.5), col=NA)
mtext("Tau among developmental stages", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauEmbAge[which(ExpData$BlanType=="SingleCopy")], 1, TypeColors[2])
BoxPlot(ExpData$TauEmbAge[which(ExpData$BlanType=="Duplicated")], 2, TypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 2.5), col=NA)
mtext("M8 D", side = 2, line = 5, cex=2)
PlotaJitterPoints(ExpData$M8.D[which(ExpData$BlanType=="SingleCopy")], 1, TypeColors[2])
PlotaJitterPoints(ExpData$M8.D[which(ExpData$BlanType=="Duplicated")], 2, TypeColors[4])
axis(1, at = c(1,2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)






layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(Tissues)*2+20), col=NA)
mtext("Tau", side = 2, line = 5, cex=2)
for(t in c(1:length(Tissues))){
	BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="SingleCopy" & ExpData$MaxTissue==Tissues[t])], t, TissueColors[t])
	BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="Duplicated" & ExpData$MaxTissue==Tissues[t])], t+length(Tissues)+3, TissueColors[t])
}
axis(1, at = c(1+length(Tissues)/2, 1+length(Tissues)+3+length(Tissues)/2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
legend("bottomright", Tissues, pch=19, text.col="black", col=TissueColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(EmbAges)*2+20), col=NA)
mtext("Tau", side = 2, line = 5, cex=2)
for(t in c(1:length(EmbAges))){
	BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="SingleCopy" & ExpData$MaxEmbAge==EmbAges[t])], t, EmbAgesColors[t])
	BoxPlot(ExpData$TauTissues[which(ExpData$BlanType=="Duplicated" & ExpData$MaxEmbAge==EmbAges[t])], t+length(EmbAges)+3, EmbAgesColors[t])
}
axis(1, at = c(1+length(EmbAges)/2, 1+length(EmbAges)+3+length(EmbAges)/2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
legend("bottomright", EmbAges, pch=19, text.col="black", col=EmbAgesColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)






layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 4.5), col=NA)
mtext("Mean adult TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Missing")], 1, TypeColors[1])
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="SingleCopy")], 2, TypeColors[2])
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Ohnolog")], 3, TypeColors[3])
BoxPlot(ExpData$MeanAdult[which(ExpData$VertebType=="Duplicated")], 4, TypeColors[4])
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 200), xlim=c(0.5, 4.5), col=NA)
mtext("Max adult TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MaxAdult[which(ExpData$VertebType=="Missing")], 1, TypeColors[1])
BoxPlot(ExpData$MaxAdult[which(ExpData$VertebType=="SingleCopy")], 2, TypeColors[2])
BoxPlot(ExpData$MaxAdult[which(ExpData$VertebType=="Ohnolog")], 3, TypeColors[3])
BoxPlot(ExpData$MaxAdult[which(ExpData$VertebType=="Duplicated")], 4, TypeColors[4])
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,200,40), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 100), xlim=c(0.5, 4.5), col=NA)
mtext("Mean embrionic TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Missing")], 1, TypeColors[1])
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="SingleCopy")], 2, TypeColors[2])
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Ohnolog")], 3, TypeColors[3])
BoxPlot(ExpData$MeanEmbr[which(ExpData$VertebType=="Duplicated")], 4, TypeColors[4])
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 200), xlim=c(0.5, 4.5), col=NA)
mtext("Max embrionic TPM", side = 2, line = 5, cex=2)
BoxPlot(ExpData$MaxEmbr[which(ExpData$VertebType=="Missing")], 1, TypeColors[1])
BoxPlot(ExpData$MaxEmbr[which(ExpData$VertebType=="SingleCopy")], 2, TypeColors[2])
BoxPlot(ExpData$MaxEmbr[which(ExpData$VertebType=="Ohnolog")], 3, TypeColors[3])
BoxPlot(ExpData$MaxEmbr[which(ExpData$VertebType=="Duplicated")], 4, TypeColors[4])
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,200,40), lwd.ticks=1, las=1, cex.axis=2)


plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 4.5), col=NA)
mtext("Tau among tissues", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Missing")], 1, TypeColors[1])
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="SingleCopy")], 2, TypeColors[2])
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Ohnolog")], 3, TypeColors[3])
BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Duplicated")], 4, TypeColors[4])
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, 4.5), col=NA)
mtext("Tau among developmental stages", side = 2, line = 5, cex=2)
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Missing")], 1, TypeColors[1])
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="SingleCopy")], 2, TypeColors[2])
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Ohnolog")], 3, TypeColors[3])
BoxPlot(ExpData$TauEmbAge[which(ExpData$VertebType=="Duplicated")], 4, TypeColors[4])
axis(1, at = c(1:4), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)




#layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
#plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(Tissues)*4+20), col=NA)
#mtext("Tau", side = 2, line = 5, cex=2)
#for(t in c(1:length(Tissues))){
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Missing" & ExpData$MaxTissue==Tissues[t])], t, TissueColors[t])
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="SingleCopy" & ExpData$MaxTissue==Tissues[t])], t+(length(Tissues)+3)*1, TissueColors[t])
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Ohnolog" & ExpData$MaxTissue==Tissues[t])], t+(length(Tissues)+3)*2, TissueColors[t])
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Duplicated" & ExpData$MaxTissue==Tissues[t])], t+(length(Tissues)+3)*3, TissueColors[t])
#}
#axis(1, at = c(1+length(Tissues)/2, 1+length(Tissues)+3+length(Tissues)/2, 1+(length(Tissues)+3)*2+length(Tissues)/2, 1+(length(Tissues)+3)*3+length(Tissues)/2), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
#axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
#legend("bottomright", Tissues, pch=19, text.col="black", col=TissueColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)

#plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, 1), xlim=c(0.5, length(EmbAges)*4+20), col=NA)
#mtext("Tau", side = 2, line = 5, cex=2)
#for(t in c(1:length(EmbAges))){
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Missing" & ExpData$MaxEmbAge==EmbAges[t])], t, EmbAgesColors[t])
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="SingleCopy" & ExpData$MaxEmbAge==EmbAges[t])], t+(length(EmbAges)+3)*1, EmbAgesColors[t])
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Ohnolog" & ExpData$MaxEmbAge==EmbAges[t])], t+(length(EmbAges)+3)*2, EmbAgesColors[t])
#	BoxPlot(ExpData$TauTissues[which(ExpData$VertebType=="Duplicated" & ExpData$MaxEmbAge==EmbAges[t])], t+(length(EmbAges)+3)*3, EmbAgesColors[t])
#}
#axis(1, at = c(1+length(EmbAges)/2, 1+length(EmbAges)+3+length(EmbAges)/2, 1+(length(EmbAges)+3)*2+length(EmbAges)/2, 1+(length(EmbAges)+3)*3+length(EmbAges)/2), labels=c("M", "SC", "O", "D"), tick=FALSE, line=2.5, las=1, cex.axis=2)
#axis(2, at = seq(0,1,.2), lwd.ticks=1, las=1, cex.axis=2)
#legend("bottomright", EmbAges, pch=19, text.col="black", col=EmbAgesColors, bty = "n", cex=1.5, xjust = 0, yjust = 0)















dev.off()












