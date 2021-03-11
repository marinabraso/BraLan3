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

######################################################################
# General parameters

Species <- c("Spur", "Bflo", "Blan", "Bbel", "Drer", "Skow", "Arub", "Ggal", "Mmus", "Hsap")
Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Blan", "Bflo", "Bbel")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Amphi, Verteb, OutDeut)


###########################################################################
# Read data

# Orthologous groups counts
Counts <- read.table(CountsFile, h=F, sep = "\t", row.names=1)
colnames(Counts) <- c("Spur", "Bflo", "Blan", "Bbel", "Drer", "Skow", "Arub", "Ggal", "Mmus", "Hsap")
Counts$Sum <- rowSums(Counts)
Counts$SumVerteb <- rowSums(Counts[,Verteb])
Counts$SumAmphi <- rowSums(Counts[,Amphi])
Counts$SumOutDeut <- rowSums(Counts[,OutDeut])
print(head(Counts))







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

layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
###########################################################################
# Number of orthologous groups per phylogenetic groups
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1000), xlim=c(0, 13), col=NA)
mtext("Number of orthologous groups", side = 2, line = 3, cex=1.5)
mtext("Orthologous groups specific of each pair of species", side = 3, line = 1, cex=1.5)
Values <- c()
PairsNames <- c()
for(s1 in c(1:(length(Verteb)-1))){
	for(s2 in c((s1+1):length(Verteb))){
		PairsNames <- c(PairsNames, paste(Verteb[s1],"\n", Verteb[s2]))
		Values <- c(Values, length(Counts[which(rowSums(Counts[,c(Verteb[s1],Verteb[s2])])==Counts$Sum),1]))
	}
}
for(s1 in c(1:(length(Amphi)-1))){
	for(s2 in c((s1+1):length(Amphi))){
		PairsNames <- c(PairsNames, paste(Amphi[s1],"\n", Amphi[s2]))
		Values <- c(Values, length(Counts[which(rowSums(Counts[,c(Amphi[s1],Amphi[s2])])==Counts$Sum),1]))
	}
}
for(s1 in c(1:(length(OutDeut)-1))){
	for(s2 in c((s1+1):length(OutDeut))){
		PairsNames <- c(PairsNames, paste(OutDeut[s1],"\n", OutDeut[s2]))
		Values <- c(Values, length(Counts[which(rowSums(Counts[,c(OutDeut[s1],OutDeut[s2])])==Counts$Sum),1]))
	}
}
width <- .5
for(v in c(1:length(Values))){
	polygon(c(v-width/2, v-width/2, v+width/2, v+width/2), c(0,Values[v],Values[v],0), col="gold")
	text(v, Values[v], labels=Values[v], pos=3)
}
axis(1, at = c(1:length(Values)), labels=NA, lwd.ticks=1, las=1, cex.axis=1)
axis(1, at = c(1:length(Values)), labels=PairsNames, tick=FALSE, line=1, las=1, cex.axis=1)
axis(2, at = seq(0,10000,2000), lwd.ticks=1, las=1, cex.axis=1)
box()

layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)
###########################################################################
# Histogram group size
HistogramGroupSize(Counts$Sum, seq(0,10000,2), c(0, 200), c(0,10), "Histogram group size")

###########################################################################
# Histogram group size vertebrate specific groups
HistogramGroupSize(Counts$Sum[which(Counts$SumVerteb>0 & Counts$SumAmphi==0 & Counts$SumOutDeut==0)], seq(0,10000,1), c(0, 50), c(0,10), "Histogram group size vertebrate specific groups")

###########################################################################
# Histogram group size amphioxus specific groups
HistogramGroupSize(Counts$Sum[which(Counts$SumVerteb==0 & Counts$SumAmphi>0 & Counts$SumOutDeut==0)], seq(0,10000,1), c(0, 50), c(0,10), "Histogram group size amphioxus specific groups")

###########################################################################
# Histogram group size chordata shared groups
HistogramGroupSize(Counts$Sum[which(Counts$SumVerteb>0 & Counts$SumAmphi>0 & Counts$SumOutDeut==0)], seq(0,10000,1), c(0, 50), c(0,10), "Histogram group size chordata shared groups")



layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(2), heights=c(1), TRUE)
###########################################################################
# Number of orthologous groups with at least 1 Blan gene member
BlanSharedGenes(Counts[which(Counts[,"Blan"]>0),], "Orthologous groups with at least 1 Blan gene member")

###########################################################################
# Number of orthologous groups with at least 2 Blan gene members
BlanSharedGenes(Counts[which(Counts[,"Blan"]>1),], "Orthologous groups with at least 2 Blan gene members")

layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1,1), heights=c(1), TRUE)
###########################################################################
# Histogram group size
HistogramGroupSize(Counts$Sum[which(Counts$Blan>0)], seq(0,10000,1), c(0, 100), c(0,10), "Histogram group size of groups with at least 1 Blan gene member")














