#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
pwd <- system("pwd", intern=TRUE)
source(paste0(pwd,"/Analysis/ComparativeGeneExpression/ComparativeGeneExpression_functions.R"))
library(BgeeDB)

######################################################################
# Files & folders
RDataFolder <- "Plots"
ResultsFolder <- "Plots/ComparativeGeneExpression"
OrganOrthBgeeFile <- paste0(pwd, "/Data/AnatomicalOntology/Julien_AnatOrht_Bgee15.tsv")
CountsFileAV <- paste0(pwd, "/Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt")
NamesFileAV <- paste0(pwd, "/Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt")
BlanGeneDataFile <- paste0(pwd, "/", RDataFolder, "/GeneDataProcessed.txt")
OhnologsSFile <- paste0(pwd, "/Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt")

######################################################################
# General parameters
TPM.threshold <- 1

setwd(paste0("./", ResultsFolder))
Species1 <- "Blan"
Species2 <- "Drer"

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum")

SName1 <- SpeciesNames[which(ShortSpeciesNames==Species1)]
SName2 <- SpeciesNames[which(ShortSpeciesNames==Species2)]

BlanTissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads", "egg", "32cells", "Blastula", "MeanEmbr")
DrerTissues <- c("blastula", "embryo", "head", "tail", "granulocyte", "brain", "heart", "muscle tissue", "mesonephros", "intestine", "pharyngeal gill", "testis", "ovary", "bone element", "liver", "mature ovarian follicle", "zone of skin", "swim bladder", "head kidney", "spleen")
DrerTissuesIDs <- c("UBERON:0000307", "UBERON:0000922", "UBERON:0000033", "UBERON:0002415", "CL:0000094", "UBERON:0000955", "UBERON:0000948", "UBERON:0002385", "UBERON:0000080", "UBERON:0000160", "UBERON:0000206", "UBERON:0000473", "UBERON:0000992", "UBERON:0001474", "UBERON:0002107", "UBERON:0003982", "UBERON:0000014", "UBERON:0006860", "UBERON:0007132", "UBERON:0002106")

BlanMatchingTissues <- c("Blastula", "MeanEmbr", "male.gonads", "female.gonads", "muscle", "neural.tube", "gut", "gills", "hepatic.diverticulum", "epidermis")
DrerMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle tissue", "brain", "intestine", "pharyngeal gill", "liver", "zone of skin")
NamesMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle", "neural", "digestive", "gills", "hepatic", "skin")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(NamesMatchingTissues))

MatchingTissues <- as.data.frame(cbind(BlanMatchingTissues, DrerMatchingTissues, NamesMatchingTissues, TissueColors))
colnames(MatchingTissues) <- c("Blan", "Drer", "Name","Color")
print(head(MatchingTissues))
MatchingTissues <- MatchingTissues[MatchingTissues$Name!="blastula" & MatchingTissues$Name!="ovary" & MatchingTissues$Name!="skin",]

###########################################################################
###########################################################################
# Read data

# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]

if(file.exists("DrerGeneData.txt")){
	# Drer Gene information previously processed from Bgee
	DrerGeneData <- read.delim("DrerGeneData.txt", h=T, stringsAsFactors=F)
}else{
	DrerBgee <- Bgee$new(species=SName2, dataType="rna_seq")
	DrerBgeeData <- getData(DrerBgee)
	write.table(DrerBgeeData, file = "DrerBgeeData.txt", quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	DrerGeneData <- as.data.frame(cbind(unique(DrerBgeeData$Gene.ID)))
	colnames(DrerGeneData) <- c("Gene")
	for(i in c(1:length(MatchingTissues[,1]))){
		print(MatchingTissues$Name[i])
		system_out <- system(paste0("cat ", "DrerBgeeData.txt", " | sed 's/\"//g' | awk -F'\\t' -v tissue=\"", MatchingTissues$Drer[i], "\" '{if($7==tissue){print $5\"\\t\"$13\"\\t\"$16}}' | sort -k1,1 | awk '{if(g!=$1){if(NR!=1){mTPM=mTPM/num;} print g\"\\t\"mTPM\"\\t\"pres; g=$1; num=1;mTPM=$2;pres=$3}else{num=num+1;mTPM=mTPM+$2;if($3==\"present\"){pres=$3}}}END{mTPM=mTPM/num; print g\"\\t\"mTPM\"\\t\"pres;}' | tail -n +2"), intern=T)
		tDrerGeneData <- read.table(text=system_out, h=F, sep = "\t")
		DrerGeneData[,paste0(MatchingTissues$Name[i],"TPM")] <- tDrerGeneData[match(DrerGeneData$Gene, tDrerGeneData[,1]),2]
		DrerGeneData[,paste0(MatchingTissues$Name[i],"Presence")] <- tDrerGeneData[match(DrerGeneData$Gene, tDrerGeneData[,1]),3]
	}
	head(DrerGeneData)
	write.table(DrerGeneData, file = "DrerGeneData.txt", quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
}

BlanGeneData <- read.table(BlanGeneDataFile, h=T, sep = "\t", row.names=1)
colnames(BlanGeneData) <- sub("^X", "", colnames(BlanGeneData))

CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=NULL)
system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAV) <- c("OG", as.character(unlist(lapply(Header[1,], as.character))))
print(head(CountsAV))
CountsAV$Ohnologs <- rep(FALSE, length(CountsAV[,1]))
CountsAV$Ohnologs[which(CountsAV$OG %in% Ohnologs)] <- rep(TRUE, length(CountsAV[which(CountsAV$OG %in% Ohnologs),1]))

system_out <- system(paste("cat ",NamesFileAV ," | tail -n +2 | awk -F '\\t' '{split($2,bl,\" \"); split($3,dr,\" \"); for(i in bl){for(j in dr){print $1\"\\t\"bl[i]\"\\t\"dr[j]}}}'"), intern=T)
GenePairs <- read.table(text=system_out, h=F, sep = "\t")
colnames(GenePairs) <- c("OG", "Gene1", "Gene2")
GenePairs <- GenePairs[which(GenePairs$Gene1 %in% BlanGeneData$Gene & GenePairs$Gene2 %in% DrerGeneData$Gene),]
print(head(GenePairs))

for(i in c(1:length(MatchingTissues[,1]))){
	GenePairs[,paste0(MatchingTissues$Name[i],1)] <- BlanGeneData[match(GenePairs$Gene1, BlanGeneData$Gene), MatchingTissues$Blan[i]]
	GenePairs[,paste0(MatchingTissues$Name[i],2)] <- DrerGeneData[match(GenePairs$Gene2, DrerGeneData$Gene), paste0(MatchingTissues$Name[i],"TPM")]
	GenePairs[,paste0(MatchingTissues$Name[i],"2p")] <- DrerGeneData[match(GenePairs$Gene2, DrerGeneData$Gene), paste0(MatchingTissues$Name[i],"Presence")]
}
print(head(GenePairs))
GenePairs$SumDom1 <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(GenePairs[x,paste0(MatchingTissues$Name, 1)]>TPM.threshold)}))
GenePairs$SumDom2 <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(GenePairs[x,paste0(MatchingTissues$Name, 2)]>TPM.threshold)}))
GenePairs$DiffDom <- GenePairs$SumDom1-GenePairs$SumDom2
GenePairs$DiffAbs <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(abs((GenePairs[x,paste0(MatchingTissues$Name, 1)]>TPM.threshold)-(GenePairs[x,paste0(MatchingTissues$Name, 2)]>TPM.threshold)))}))
GenePairs$Ohnologs <- rep(FALSE, length(GenePairs[,1]))
GenePairs$Ohnologs[which(GenePairs$OG %in% Ohnologs)] <- rep(TRUE, length(GenePairs[which(GenePairs$OG %in% Ohnologs),1]))
GenePairs.121 <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer==1)]),]
GenePairs.12m <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer>=2)]),]
GenePairs.m21 <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan>=2 & CountsAV$Drer==1)]),]

SumUnionOfPatterns <- function(OG, spnum){
	return(sum(colSums(GenePairs[which(GenePairs$OG == OG),paste0(MatchingTissues$Name,spnum)]>TPM.threshold)>0))
}
DiffAbsUnionOfPatterns <- function(OG){
	pattern1 <- colSums(GenePairs[which(GenePairs$OG == OG),paste0(MatchingTissues$Name,1)]>TPM.threshold)>0
	pattern2 <- colSums(GenePairs[which(GenePairs$OG == OG),paste0(MatchingTissues$Name,2)]>TPM.threshold)>0
	return(sum(abs(pattern1-pattern2)))
}
CountsAV$SumUnionDom1 <- unlist(lapply(CountsAV$OG, SumUnionOfPatterns, 1))
CountsAV$SumUnionDom2 <- unlist(lapply(CountsAV$OG, SumUnionOfPatterns, 2))
CountsAV$DiffUnionDom <- CountsAV$SumUnionDom1-CountsAV$SumUnionDom2
CountsAV$DiffAbsUnion <- unlist(lapply(CountsAV$OG, DiffAbsUnionOfPatterns))
print(summary(CountsAV$DiffUnionDom))
print(summary(CountsAV$DiffAbsUnion))


pdf("ComparativeGeneExpression.pdf", width=20, height=20)
par(mar=c(5,5,2,2),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)


breaks <- seq(-.5,length(MatchingTissues[,1])+.5,1)
hist(GenePairs.121$SumDom1, breaks=breaks)
hist(GenePairs.121$SumDom2, breaks=breaks)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,log(3000)), xlim=c(-.5,length(MatchingTissues[,1])+.5), col=NA)
mtext("Frequency", side = 2, line = 3, cex=1.5)
mtext("Differences between patterns", side = 1, line = 3, cex=1.5)
h <- hist(GenePairs.121$DiffAbs, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "grey60")
h <- hist(GenePairs.12m$DiffAbs, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "royalblue1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan==1 & CountsAV$Drer>1)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "royalblue4")
h <- hist(GenePairs.m21$DiffAbs, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "orangered1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan>1 & CountsAV$Drer==1)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "orangered4")
h <- hist(GenePairs.12m$DiffAbs[which(GenePairs.12m$Ohnologs==TRUE)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "forestgreen1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan==1 & CountsAV$Drer>1 & GenePairs.12m$Ohnologs==TRUE)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "forestgreen4")
h <- hist(GenePairs.12m$DiffAbs[which(GenePairs.12m$Ohnologs==FALSE)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "gold1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan==1 & CountsAV$Drer>1 & GenePairs.12m$Ohnologs==FALSE)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "gold4")
axis(1, at = c(0:length(MatchingTissues[,1])), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = c(0,log(c(1,10,100,1000,10000))), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
box()

breaks <- seq(-length(MatchingTissues[,1])-.5,length(MatchingTissues[,1])+.5,1)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,log(3000)), xlim=c(-length(MatchingTissues[,1])+.5,length(MatchingTissues[,1])+.5), col=NA)
mtext("Frequency", side = 2, line = 3, cex=1.5)
mtext("Difference in number of domains", side = 1, line = 3, cex=1.5)
h <- hist(GenePairs.121$DiffDom, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "grey60")
h <- hist(GenePairs.12m$DiffDom, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "royalblue1")
h <- hist(GenePairs.m21$DiffDom, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "orangered1")
h <- hist(CountsAV$DiffUnionDom[which(CountsAV$Blan==1 & CountsAV$Drer>1)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "royalblue4")
h <- hist(CountsAV$DiffUnionDom[which(CountsAV$Blan>1 & CountsAV$Drer==1)], breaks=breaks, plot=FALSE)
AddHistlogYaxis(h, "orangered4")
axis(1, at = c(-length(MatchingTissues[,1]):length(MatchingTissues[,1])), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = c(0,log(c(1,10,100,1000,10000))), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
box()

dev.off()


plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,7000), xlim=c(0,length(MatchingTissues[,1])+1), col=NA)
mtext("Density", side = 2, line = 3, cex=1.5)
mtext("Tissues with shared expression", side = 1, line = 3, cex=1.5)
for(i in c(1:length(MatchingTissues[,1]))){
	print(MatchingTissues$Name[i])
	tmp <- rowSums(BlanGeneData[which(BlanGeneData[,MatchingTissues$Blan[i]]>TPM.threshold),MatchingTissues$Blan]>TPM.threshold)
	h <- hist(tmp, breaks=seq(0.5,length(MatchingTissues[,1])+.5,1), plot=FALSE)
	lines(h$mids, h$counts, col=MatchingTissues$Color[i],lwd=3)
	points(h$mids, h$counts, col=MatchingTissues$Color[i],pch=16)
}
axis(1, at = c(0:length(MatchingTissues[,1])), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = seq(0, 10000, 2000), lwd.ticks=1, las=1, cex.axis=1)
legend("topleft", MatchingTissues$Name, pch=15, col=MatchingTissues$Color, bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
box()


plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,7000), xlim=c(0,length(MatchingTissues[,1])+1), col=NA)
mtext("Density", side = 2, line = 3, cex=1.5)
mtext("Tissues with shared expression", side = 1, line = 3, cex=1.5)
for(i in c(1:length(MatchingTissues[,1]))){
	print(MatchingTissues$Name[i])
	tmp <- rowSums(DrerGeneData[which(DrerGeneData[,paste0(MatchingTissues$Name[i],"TPM")]>TPM.threshold),paste0(MatchingTissues$Name,"TPM")]>TPM.threshold)
	h <- hist(tmp, breaks=seq(0.5,length(MatchingTissues[,1])+.5,1), plot=FALSE)
	lines(h$mids, h$counts, col=MatchingTissues$Color[i],lwd=3)
	points(h$mids, h$counts, col=MatchingTissues$Color[i],pch=16)
}
axis(1, at = c(0:length(MatchingTissues[,1])), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = seq(0, 10000, 2000), lwd.ticks=1, las=1, cex.axis=1)
legend("topleft", MatchingTissues$Name, pch=15, col=MatchingTissues$Color, bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
box()



Corr <- c()
CorrSC <- c()
Corr21 <- c()
Corr12 <- c()
for(i in c(1:length(MatchingTissues[,1]))){
	Corr <- c(Corr, GetSpearmanCorr(GenePairs, MatchingTissues$Name[i]))
	CorrSC <- c(CorrSC, GetSpearmanCorr(GenePairsSC, MatchingTissues$Name[i]))
	Corr21 <- c(Corr21, GetSpearmanCorr(GenePairs.m21, MatchingTissues$Name[i]))
	Corr12 <- c(Corr12, GetSpearmanCorr(GenePairs.12m, MatchingTissues$Name[i]))
}
MatchingTissues$Corr <- Corr
MatchingTissues$CorrSC <- CorrSC
MatchingTissues$Corr21 <- Corr21
MatchingTissues$Corr12 <- Corr12





layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(2), heights=c(1), TRUE)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,1), xlim=c(0,length(MatchingTissues[,1])+1), col=NA)
mtext("Spearman correlation", side = 2, line = 5, cex=1.5)
points(c(1:length(MatchingTissues[,1])), MatchingTissues$Corr, col="black", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1])), MatchingTissues$CorrSC, col="grey60", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1])), MatchingTissues$Corr21, col="royalblue1", pch=16, cex=2)
points(c(1:length(MatchingTissues[,1])), MatchingTissues$Corr12, col="orangered4", pch=16, cex=2)
axis(1, at = c(1:length(MatchingTissues[,1])), labels=MatchingTissues$Name, lwd.ticks=1, las=1, cex.axis=1.2)
axis(2, at = seq(0,1,.05), lwd.ticks=1, las=1, cex.axis=1.2)



dev.off()












