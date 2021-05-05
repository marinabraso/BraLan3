#!/usr/bin/env Rscript



######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))

######################################################################
# Files & folders

ResultsFolder <- "Plots/FindOrthologs"
CountsFile <- "Results/FindOrthologs/broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFile <- "Results/FindOrthologs/broccoli/dir_step3/table_OGs_protein_names.txt"
OnhSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"

######################################################################
# General parameters

Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Blan", "Bflo", "Bbel")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Amphi, Verteb, OutDeut)
colfunc <- colorRampPalette(c("darkred", "firebrick1"))
CNcolors <- c("grey80", "gold", colfunc(3))
VertebTypeColors <- modif_alpha(c("deepskyblue2", "gold", "firebrick2", "firebrick4"), .7)

###########################################################################
# Read data

# Orthologous groups counts
Counts <- read.table(CountsFile, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFile, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(Counts) <- as.character(unlist(lapply(Header[1,], as.character)))

OnhOG.S <- read.table(OnhSFile, h=F)[,1]

Counts$Sum <- rowSums(Counts)
Counts$SumVerteb <- rowSums(Counts[,Verteb])
Counts$SumAmphi <- rowSums(Counts[,Amphi])
Counts$SumOutDeut <- rowSums(Counts[,OutDeut])
Counts$VertebType <- rep("Missing", length(Counts[,1]))
Counts$VertebType[which(Counts$SumVerteb>0)] <- rep("SingleCopy", length(Counts$VertebType[which(Counts$SumVerteb>0)]))
Counts$VertebType[which(apply(Counts[,Verteb], 1, max)>1)] <- rep("Duplicated", length(Counts$VertebType[which(apply(Counts[,Verteb], 1, max)>1)]))
Counts$VertebType[which(rownames(Counts) %in% OnhOG.S)] <- rep("Ohnolog", length(Counts$VertebType[which(rownames(Counts) %in% OnhOG.S)]))
print(head(Counts))
Counts$BlanType <- rep("Missing", length(Counts[,1]))
Counts$BlanType[which(Counts$Blan>0)] <- rep("SingleCopy", length(Counts$BlanType[which(Counts$Blan>0)]))
Counts$BlanType[which(Counts$Blan>1)] <- rep("Duplicated", length(Counts$BlanType[which(Counts$Blan>1)]))


###########################################################################
###########################################################################
# Plotting
###########################################################################
###########################################################################
pdf(paste(ResultsFolder, "/GroupsOfOrthologsOverview.pdf", sep=""), width=20, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')

layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)
###########################################################################
# Number of orthologous groups per phylogenetic groups
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,10000), xlim=c(0, 6), col=NA)
mtext("Number of orthologous groups", side = 2, line = 3, cex=1.5)
mtext("Orthologous groups per phylogenetic groups", side = 3, line = 1, cex=1.5)

PhyGroups <- c("Vertebrate\nspecific", "Amphioxus\nspecific", "Chordate\nshared", "Non-chordate\nspecific", "Shared by\nall groups")
Values <- length(Counts[which(Counts$SumVerteb>0 & Counts$SumAmphi==0 & Counts$SumOutDeut==0),1])
Values <- c(Values, length(Counts[which(Counts$SumVerteb==0 & Counts$SumAmphi>0 & Counts$SumOutDeut==0),1]))
Values <- c(Values, length(Counts[which(Counts$SumVerteb>0 & Counts$SumAmphi>0 & Counts$SumOutDeut==0),1]))
Values <- c(Values, length(Counts[which(Counts$SumVerteb==0 & Counts$SumAmphi==0 & Counts$SumOutDeut>0),1]))
Values <- c(Values, length(Counts[which(Counts$SumVerteb>0 & Counts$SumAmphi>0 & Counts$SumOutDeut>0),1]))
width <- .5
for(v in c(1:length(Values))){
	polygon(c(v-width/2, v-width/2, v+width/2, v+width/2), c(0,Values[v],Values[v],0), col="gold")
	text(v, Values[v], labels=Values[v], pos=3)
}
axis(1, at = c(1:length(Values)), labels=NA, lwd.ticks=1, las=1, cex.axis=1)
axis(1, at = c(1:length(Values)), labels=PhyGroups, tick=FALSE, line=1, las=1, cex.axis=1)
axis(2, at = seq(0,10000,2000), lwd.ticks=1, las=1, cex.axis=1)
box()

CountBlanPresent <- Counts[which(Counts$BlanType!="Missing"),]
chisq.test(CountBlanPresent$BlanType == "SingleCopy", CountBlanPresent$VertebType == "SingleCopy")
sum(CountBlanPresent$BlanType == "SingleCopy")*sum(CountBlanPresent$VertebType == "SingleCopy")/length(CountBlanPresent[,1])
sum(CountBlanPresent$BlanType == "SingleCopy" & CountBlanPresent$VertebType == "SingleCopy")

chisq.test(CountBlanPresent$BlanType == "Duplicated", CountBlanPresent$VertebType == "Ohnolog")
sum(CountBlanPresent$BlanType == "Duplicated")*sum(CountBlanPresent$VertebType == "Ohnolog")/length(CountBlanPresent[,1])
sum(CountBlanPresent$BlanType == "Duplicated" & CountBlanPresent$VertebType == "Ohnolog")

chisq.test(CountBlanPresent$BlanType == "Duplicated", CountBlanPresent$VertebType == "Duplicated")
sum(CountBlanPresent$BlanType == "Duplicated")*sum(CountBlanPresent$VertebType == "Duplicated")/length(CountBlanPresent[,1])
sum(CountBlanPresent$BlanType == "Duplicated" & CountBlanPresent$VertebType == "Duplicated")


CountBlanPVertebD <- CountBlanPresent[which(CountBlanPresent$VertebType!="Missing" & CountBlanPresent$VertebType!="SingleCopy"),]
chisq.test(CountBlanPVertebD$BlanType == "Duplicated", CountBlanPVertebD$VertebType == "Ohnolog")
sum(CountBlanPVertebD$BlanType == "Duplicated")*sum(CountBlanPVertebD$VertebType == "Ohnolog")/length(CountBlanPVertebD[,1])
sum(CountBlanPVertebD$BlanType == "Duplicated" & CountBlanPVertebD$VertebType == "Ohnolog")


layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 3.5), col=NA)
mtext("", side = 2, line = 3, cex=1.5)

PlotColumnVertebrateType(Counts[which(Counts$Blan==1),], 1, VertebTypeColors)
PlotColumnVertebrateType(Counts[which(Counts$Blan>=2),], 2, VertebTypeColors)

axis(1, at = c(1:2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,19,19,19,19), text.col="black", col=c(NA, VertebTypeColors), bty = "n", cex=2, xjust = 0, yjust = 0)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 3.5), col=NA)
mtext("", side = 2, line = 3, cex=1.5)

PlotColumnVertebrateType(CountBlanPVertebD[which(CountBlanPVertebD$Blan==1),], 1, VertebTypeColors)
PlotColumnVertebrateType(CountBlanPVertebD[which(CountBlanPVertebD$Blan>=2),], 2, VertebTypeColors)

axis(1, at = c(1:2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,19,19,19,19), text.col="black", col=c(NA, VertebTypeColors), bty = "n", cex=2, xjust = 0, yjust = 0)


Counts <- Counts[which(Counts$SumOutDeut==0),]
CountBlanPresent <- Counts[which(Counts$BlanType!="Missing"),]
chisq.test(CountBlanPresent$BlanType == "SingleCopy", CountBlanPresent$VertebType == "SingleCopy")
sum(CountBlanPresent$BlanType == "SingleCopy")*sum(CountBlanPresent$VertebType == "SingleCopy")/length(CountBlanPresent[,1])
sum(CountBlanPresent$BlanType == "SingleCopy" & CountBlanPresent$VertebType == "SingleCopy")

chisq.test(CountBlanPresent$BlanType == "Duplicated", CountBlanPresent$VertebType == "Ohnolog")
sum(CountBlanPresent$BlanType == "Duplicated")*sum(CountBlanPresent$VertebType == "Ohnolog")/length(CountBlanPresent[,1])
sum(CountBlanPresent$BlanType == "Duplicated" & CountBlanPresent$VertebType == "Ohnolog")

chisq.test(CountBlanPresent$BlanType == "Duplicated", CountBlanPresent$VertebType == "Duplicated")
sum(CountBlanPresent$BlanType == "Duplicated")*sum(CountBlanPresent$VertebType == "Duplicated")/length(CountBlanPresent[,1])
sum(CountBlanPresent$BlanType == "Duplicated" & CountBlanPresent$VertebType == "Duplicated")


CountBlanPVertebD <- CountBlanPresent[which(CountBlanPresent$VertebType!="Missing" & CountBlanPresent$VertebType!="SingleCopy"),]
chisq.test(CountBlanPVertebD$BlanType == "Duplicated", CountBlanPVertebD$VertebType == "Ohnolog")
sum(CountBlanPVertebD$BlanType == "Duplicated")*sum(CountBlanPVertebD$VertebType == "Ohnolog")/length(CountBlanPVertebD[,1])
sum(CountBlanPVertebD$BlanType == "Duplicated" & CountBlanPVertebD$VertebType == "Ohnolog")


layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 3.5), col=NA)
mtext("", side = 2, line = 3, cex=1.5)

PlotColumnVertebrateType(Counts[which(Counts$Blan==1),], 1, VertebTypeColors)
PlotColumnVertebrateType(Counts[which(Counts$Blan>=2),], 2, VertebTypeColors)

axis(1, at = c(1:2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,19,19,19,19), text.col="black", col=c(NA, VertebTypeColors), bty = "n", cex=2, xjust = 0, yjust = 0)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,100), xlim=c(0.5, 3.5), col=NA)
mtext("", side = 2, line = 3, cex=1.5)

PlotColumnVertebrateType(CountBlanPVertebD[which(CountBlanPVertebD$Blan==1),], 1, VertebTypeColors)
PlotColumnVertebrateType(CountBlanPVertebD[which(CountBlanPVertebD$Blan>=2),], 2, VertebTypeColors)

axis(1, at = c(1:2), labels=c("Single copy\nin Bralan3", "Duplicated\nin Bralan3"), tick=FALSE, line=2.5, las=1, cex.axis=2)
axis(2, at = seq(0,100,20), lwd.ticks=1, las=1, cex.axis=1.5)
legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,19,19,19,19), text.col="black", col=c(NA, VertebTypeColors), bty = "n", cex=2, xjust = 0, yjust = 0)

	


