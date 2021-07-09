#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
pwd <- system("pwd", intern=TRUE)
script <- sub(".*=", "", commandArgs()[4])
source(paste0(dirname(script), "/", substr(basename(script),1, nchar(basename(script))-2), "_functions.r", sep=""))

######################################################################
# Files & folders
RDataFolder <- "Plots"
ResultsFolder <- "Plots/ComparativeGeneExpression"
OrganOrthBgeeFile <- paste0(pwd, "/Data/AnatomicalOntology/Julien_AnatOrht_Bgee15.tsv")
CountsFileAV <- paste0(pwd, "/Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt")
NamesFileAV <- paste0(pwd, "/Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt")
BlanGeneDataFile <- paste0(pwd, "/", RDataFolder, "/GeneDataProcessed.txt")
OhnologsSFile <- paste0(pwd, "/Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt")
DrerGeneDataFile <- paste0(pwd, "/", ResultsFolder, "/DrerGeneData.txt")
DrerBgeeDataFile <- paste0(pwd, "/", ResultsFolder, "/DrerBgeeData.txt")

######################################################################
# General parameters
TPM.threshold <- 1
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
MatchingTissues <- MatchingTissues[MatchingTissues$Name!="blastula" & MatchingTissues$Name!="ovary" & MatchingTissues$Name!="skin",]

###########################################################################
###########################################################################
# Read data

# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]

# Read Drer gene information previously processed or precess it and write itinto a file
DrerGeneData <- prepare_DrerGeneData(DrerGeneDataFile, DrerBgeeDataFile, SName2, MatchingTissues, ResultsFolder, pwd)
print(head(DrerGeneData))

# Read Blan gene data file 
BlanGeneData <- read.table(BlanGeneDataFile, h=T, sep = "\t", row.names=1)
colnames(BlanGeneData) <- sub("^X", "", colnames(BlanGeneData))

# Read orthologs counts information
CountsAV <- prepare_CountsAV(CountsFileAV, Ohnologs)
print(head(CountsAV))

# Prepare Blan - Drer gene pairs
prepare_GenePairs(NamesFileAV, BlanGeneData, DrerGeneData, MatchingTissues, TPM.threshold)

prepare_GenePairs <- function(nfile, blangd, drergd, tissues, tpmthresh){
	system_out <- system(paste("cat ",nfile ," | tail -n +2 | awk -F '\\t' '{split($2,bl,\" \"); split($3,dr,\" \"); for(i in bl){for(j in dr){print $1\"\\t\"bl[i]\"\\t\"dr[j]}}}'"), intern=T)
	GenePairs <- read.table(text=system_out, h=F, sep = "\t")
	colnames(GenePairs) <- c("OG", "Gene1", "Gene2")
	GenePairs <- GenePairs[which(GenePairs$Gene1 %in% blangd$Gene & GenePairs$Gene2 %in% drergd$Gene),]

	# Getting tissue gene expression data for each pair 
	for(i in c(1:length(tissues[,1]))){
		GenePairs[,paste0(tissues$Name[i],1)] <- blangd[match(GenePairs$Gene1, blangd$Gene), tissues$Blan[i]]
		GenePairs[,paste0(tissues$Name[i],2)] <- drergd[match(GenePairs$Gene2, drergd$Gene), paste0(tissues$Name[i],"TPM")]
		GenePairs[,paste0(tissues$Name[i],"2p")] <- drergd[match(GenePairs$Gene2, drergd$Gene), paste0(tissues$Name[i],"Presence")]
	}

	# Calculate sum of expressed tissues (domains) in each species and its difference
	GenePairs$SumDom1 <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(GenePairs[x,paste0(tissues$Name, 1)]>TPM.threshold)}))
	GenePairs$SumDom2 <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(GenePairs[x,paste0(tissues$Name, 2)]>TPM.threshold)}))
	GenePairs$DiffDom <- GenePairs$SumDom1-GenePairs$SumDom2

	# Calculate absolute difference between species expression profiles
	GenePairs$DiffAbs <- unlist(lapply(c(1:length(GenePairs[,1])), function(x){sum(abs((GenePairs[x,paste0(tissues$Name, 1)]>TPM.threshold)-(GenePairs[x,paste0(tissues$Name, 2)]>tpmthresh)))}))	
}

# Label genepairs that are ohnologs
GenePairs$Ohnologs <- rep(FALSE, length(GenePairs[,1]))
GenePairs$Ohnologs[which(GenePairs$OG %in% Ohnologs)] <- rep(TRUE, length(GenePairs[which(GenePairs$OG %in% Ohnologs),1]))

# Subdivisions of gene pairs (1 to 1 orthologs, 1 to many, many to 1)
GenePairs.121 <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer==1)]),]
GenePairs.12m <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan==1 & CountsAV$Drer>=2)]),]
GenePairs.m21 <- GenePairs[which(GenePairs$OG %in% CountsAV$OG[which(CountsAV$Blan>=2 & CountsAV$Drer==1)]),]

# 
CountsAV$SumUnionDom1 <- unlist(lapply(CountsAV$OG, sum_union_of_patterns, 1, GenePairs, MatchingTissues, TPM.threshold))
CountsAV$SumUnionDom2 <- unlist(lapply(CountsAV$OG, sum_union_of_patterns, 2, GenePairs, MatchingTissues, TPM.threshold))
CountsAV$DiffUnionDom <- CountsAV$SumUnionDom1-CountsAV$SumUnionDom2
CountsAV$DiffAbsUnion <- unlist(lapply(CountsAV$OG, absolute_difference_union_of_patterns, GenePairs, MatchingTissues, TPM.threshold))
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
add_hist_logYaxis(h, "grey60")
h <- hist(GenePairs.12m$DiffAbs, breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "royalblue1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan==1 & CountsAV$Drer>1)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "royalblue4")
h <- hist(GenePairs.m21$DiffAbs, breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "orangered1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan>1 & CountsAV$Drer==1)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "orangered4")
h <- hist(GenePairs.12m$DiffAbs[which(GenePairs.12m$Ohnologs==TRUE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "olivedrab1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan==1 & CountsAV$Drer>1 & CountsAV$Ohnologs==TRUE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "olivedrab4")
h <- hist(GenePairs.12m$DiffAbs[which(GenePairs.12m$Ohnologs==FALSE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "gold1")
h <- hist(CountsAV$DiffAbsUnion[which(CountsAV$Blan==1 & CountsAV$Drer>1 & CountsAV$Ohnologs==FALSE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "gold4")
axis(1, at = c(0:length(MatchingTissues[,1])), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = c(0,log(c(1,10,100,1000,10000))), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
box()

breaks <- seq(-length(MatchingTissues[,1])-.5,length(MatchingTissues[,1])+.5,1)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,log(3000)), xlim=c(-length(MatchingTissues[,1])+.5,length(MatchingTissues[,1])+.5), col=NA)
mtext("Frequency", side = 2, line = 3, cex=1.5)
mtext("Difference in number of domains", side = 1, line = 3, cex=1.5)
h <- hist(GenePairs.121$DiffDom, breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "grey60")
h <- hist(GenePairs.12m$DiffDom, breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "royalblue1")
h <- hist(CountsAV$DiffUnionDom[which(CountsAV$Blan==1 & CountsAV$Drer>1)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "royalblue4")
h <- hist(GenePairs.m21$DiffDom, breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "orangered1")
h <- hist(CountsAV$DiffUnionDom[which(CountsAV$Blan>1 & CountsAV$Drer==1)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "orangered4")
h <- hist(GenePairs.12m$DiffDom[which(GenePairs.12m$Ohnologs==TRUE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "olivedrab1")
h <- hist(CountsAV$DiffUnionDom[which(CountsAV$Blan==1 & CountsAV$Drer>1 & CountsAV$Ohnologs==TRUE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "olivedrab4")
h <- hist(GenePairs.12m$DiffDom[which(GenePairs.12m$Ohnologs==FALSE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "gold1")
h <- hist(CountsAV$DiffUnionDom[which(CountsAV$Blan==1 & CountsAV$Drer>1 & CountsAV$Ohnologs==FALSE)], breaks=breaks, plot=FALSE)
add_hist_logYaxis(h, "gold4")
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
	Corr <- c(Corr, get_spearman_corr(GenePairs, MatchingTissues$Name[i]))
	CorrSC <- c(CorrSC, get_spearman_corr(GenePairsSC, MatchingTissues$Name[i]))
	Corr21 <- c(Corr21, get_spearman_corr(GenePairs.m21, MatchingTissues$Name[i]))
	Corr12 <- c(Corr12, get_spearman_corr(GenePairs.12m, MatchingTissues$Name[i]))
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












