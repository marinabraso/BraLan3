#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
library(viridis)

######################################################################
# Files & folders

ResultsFolder <- "Plots/GeneOntology"
GOlistsFolder <- "Results/GeneOntology"
CountsFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
HumanOhnologsSFile <- "Data/Ohnologs/hsapiens.Families.Strict.2R.ENS.list"
MainGOlistFile <- "Data/GeneOntology/MainGOMF_supclass.list"
GOlistFile <- "Results/GeneOntology/GOlist_MF.list"
Ens2GOFile <- "Data/GeneOntology/HumanGeneToGoTerm_Ensembl.txt"
OBOFile <- "Data/GeneOntology/go.obo"
HsapGenelist <- "Results/FilteringGeneSets/GTFs/Homo_sapiens.GRCh38.list"
BlanGenelist <- "Results/FilteringGeneSets/GTFs/Branchiostoma_lanceolatum.BraLan3.list"
ExtractGOGenesScript <- "Scripts/GeneOntology/GetAllGenesFromTerm_from_OBO_ENS.sh"

######################################################################
# General parameters
minGOsize <- 200

###########################################################################
# Read data

# List of main GO terms
mGO <- read.table(MainGOlistFile, h=F, sep = '\t')
colnames(mGO) <- c("GO", "Class", "Name")
colfunc <- colorRampPalette(c(viridis(3), "firebrick"))
GOClass <- as.data.frame(cbind(unique(mGO$Class), colfunc(length(unique(mGO$Class)))))
colnames(GOClass) <- c("Class", "Color")
mGO$ClassColors <- GOClass$Color[match(mGO$Class, GOClass$Class)]
print(head(mGO))
# Complete list of  GO terms
GO <- read.delim(GOlistFile, h=F, sep = '\t')
colnames(GO) <- c("GO", "Class", "Name")
GO$ClassColors <- GOClass$Color[match(GO$Class, GOClass$Class)]
print(head(GO))

# List of onhologs
HumanOhnologs <- read.table(HumanOhnologsSFile, h=F)[,1]

# List of Hsap genes
HsapGenes <- read.table(HsapGenelist, h=F)
colnames(HsapGenes) <- "HGene"

# List of Blan genes
BlanGenes <- read.table(BlanGenelist, h=F)
colnames(BlanGenes) <- "BGene"

# Orthologous groups counts (amphioxus + vertebrates)
CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]*_\\([a-z][a-z][a-z]\\)[A-Za-z\\.0-9_]*/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAV) <- as.character(unlist(lapply(Header[1,], as.character)))
print(head(CountsAV))

# OG to Blan gene
system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($2,a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
OG2BlanGene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2BlanGene) <- c("OG", "BGene")

# OG to Blan gene
system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($5,a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
OG2HsapGene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2HsapGene) <- c("OG", "HGene")

BlanGenes$DupOrNot <- rep("SC", length(BlanGenes[,1]))
BlanGenes$DupOrNot[which(BlanGenes$BGene %in% OG2BlanGene$BGene[which(OG2BlanGene$OG %in% rownames(CountsAV[which(CountsAV$Blan>=2),]))])] <- rep("D", length(BlanGenes$DupOrNot[which(BlanGenes$BGene %in% OG2BlanGene$BGene[which(OG2BlanGene$OG %in% rownames(CountsAV[which(CountsAV$Blan>=2),]))])]))
BlanGenes$OrthOhnOrNot <- rep(FALSE, length(BlanGenes[,1]))
BGenelist <- OG2BlanGene$BGene[which(OG2BlanGene$OG %in% OG2HsapGene$OG[which(OG2HsapGene$HGene %in% HumanOhnologs)])]
BlanGenes$OrthOhnOrNot[which(BlanGenes$BGene %in% OG2BlanGene$BGene[which(OG2BlanGene$OG %in% OG2HsapGene$OG[which(OG2HsapGene$HGene %in% HumanOhnologs)])])] <- rep(TRUE, length(BlanGenes$OrthOhnOrNot[which(BlanGenes$BGene %in% OG2BlanGene$BGene[which(OG2BlanGene$OG %in% OG2HsapGene$OG[which(OG2HsapGene$HGene %in% HumanOhnologs)])])]))
print(head(BlanGenes))
print(table(BlanGenes$DupOrNot))
print(table(BlanGenes$OrthOhnOrNot))

HsapGenes$DupOrNot <- rep("SC", length(HsapGenes[,1]))
HsapGenes$DupOrNot[which(HsapGenes$HGene %in% OG2HsapGene$HGene[which(OG2HsapGene$OG %in% rownames(CountsAV[which(CountsAV$Hsap>=2),]))])] <- rep("D", length(HsapGenes$DupOrNot[which(HsapGenes$HGene %in% OG2HsapGene$HGene[which(OG2HsapGene$OG %in% rownames(CountsAV[which(CountsAV$Hsap>=2),]))])]))
HsapGenes$DupOrNot[which(HsapGenes$HGene %in% HumanOhnologs)] <- rep("O", length(HsapGenes$DupOrNot[which(HsapGenes$HGene %in% HumanOhnologs)]))
print(head(HsapGenes))
print(table(HsapGenes$DupOrNot))

pdf(paste(ResultsFolder, "/BlanHsap_GeneOntology_Duplicability.pdf", sep=""), width=20, height=10)
par(mar=c(10,10,3,3),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)

HsapSize.GO <- c()
HsapD.GO <- c()
HsapO.GO <- c()
HsapT.GO <- c()
BlanD.GO <- c()
BlanO.GO <- c()
BlanDO.GO <- c()
BlanT.GO <- c()
for(g in c(1:length(GO[,1]))){
	fileinfo = file.info(paste(GOlistsFolder, "/", GO$GO[g], ".txt", sep=""))
	empty = rownames(fileinfo[fileinfo$size == 0, ])
	if(length(empty)==0){
		HGenelist <-  read.table(paste(GOlistsFolder, "/", GO$GO[g], ".txt", sep=""), h=F)[,1]
	}else{
		HGenelist <- c()
	}
	HsapSize.GO <- c(HsapSize.GO, length(HGenelist))
	if(length(HGenelist) >= minGOsize){
		HsapD.GO <- c(HsapD.GO, length(HsapGenes[which(HsapGenes$DupOrNot=="D" & HsapGenes$HGene %in% HGenelist),1]))
		HsapO.GO <- c(HsapO.GO, length(HsapGenes[which(HsapGenes$DupOrNot=="O" & HsapGenes$HGene %in% HGenelist),1]))
		HsapT.GO <- c(HsapT.GO, length(HsapGenes[which(HsapGenes$HGene %in% HGenelist),1]))
		BGenelist <- OG2BlanGene$BGene[which(OG2BlanGene$OG %in% OG2HsapGene$OG[which(OG2HsapGene$HGene %in% HGenelist)])]
		BlanD.GO <- c(BlanD.GO, length(BlanGenes[which(BlanGenes$DupOrNot=="D" & BlanGenes$BGene %in% BGenelist),1]))
		BlanDO.GO <- c(BlanDO.GO, length(BlanGenes[which(BlanGenes$DupOrNot=="D" & BlanGenes$OrthOhnOrNot==TRUE & BlanGenes$BGene %in% BGenelist),1]))
		BlanO.GO <- c(BlanO.GO, length(BlanGenes[which(BlanGenes$OrthOhnOrNot==TRUE & BlanGenes$BGene %in% BGenelist),1]))
		BlanT.GO <- c(BlanT.GO, length(BlanGenes[which(BlanGenes$BGene %in% BGenelist),1]))
	}else{
		HsapD.GO <- c(HsapD.GO, 0)
		HsapO.GO <- c(HsapO.GO, 0)
		HsapT.GO <- c(HsapT.GO, 0)
		BlanD.GO <- c(BlanD.GO, 0)
		BlanDO.GO <- c(BlanDO.GO, 0)
		BlanO.GO <- c(BlanO.GO, 0)
		BlanT.GO <- c(BlanT.GO, 0)
	}
}
GO$HsapSize <- HsapSize.GO
GO$HsapD <- HsapD.GO
GO$HsapO <- HsapO.GO
GO$HsapT <- HsapT.GO
GO$HsapDprop <- HsapD.GO/HsapT.GO*100
GO$HsapOprop <- HsapO.GO/HsapT.GO*100
GO$HsapDOprop <- (HsapD.GO+HsapO.GO)/HsapT.GO*100
GO$BlanD <- BlanD.GO
GO$BlanO <- BlanO.GO
GO$BlanDO <- BlanDO.GO
GO$BlanT <- BlanT.GO
GO$BlanDprop <- BlanD.GO/BlanT.GO*100
GO$BlanOprop <- BlanO.GO/BlanT.GO*100
GO$BlanDOprop <- BlanDO.GO/BlanT.GO*100
GO <- GO[which(!is.na(GO$BlanDprop) & !is.na(GO$HsapDprop) & !is.na(GO$HsapOprop)),]
tableClass <- table(GO$Class)
print(tableClass)
GO$NumInClass <- tableClass[match(GO$Class, names(tableClass))]
GO <- unique(GO[order(GO$NumInClass, decreasing = TRUE),])
GOClass <- as.data.frame(cbind(unique(GO$Class), colfunc(length(unique(GO$Class)))))
colnames(GOClass) <- c("Class", "Color")
GO$ClassColors <- GOClass$Color[match(GO$Class, GOClass$Class)]
print(GO[order(GO$HsapDprop),c("Name", "Class", "HsapDprop", "BlanDprop")])

ScatterPercentagePlot(GO$HsapDprop, GO$BlanDprop, GO$ClassColors, "% of small scale duplicated genes\nH. sapiens", "B. lanceolatum\n% of duplicated genes", c(0,100), c(0,100))
ScatterPercentagePlot(GO$HsapOprop, GO$BlanDprop, GO$ClassColors, "% of ohnolog genes\nH. sapiens", "B. lanceolatum\n% of duplicated genes", c(0,30), c(0,100))
ScatterPercentagePlot(GO$HsapDOprop, GO$BlanDprop, GO$ClassColors, "% of duplicated + ohnolog genes in\nsapiens", "B. lanceolatum\n% of duplicated genes", c(0,100), c(0,100))
ScatterPercentagePlot(GO$HsapOprop, GO$HsapDprop, GO$ClassColors, "% of ohnolog genes in\nsapiens", "% of duplicated genes\nH. sapiens", c(0,30), c(0,100))

plot.new()
legendvec <- unique(GO[,c("Class", "ClassColors")])
legend("center", legendvec[,1], pch=15, col=legendvec[,2], bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
plot.new()

print(GO[which(GO$BlanDprop>95),c("GO", "Name", "Class", "HsapDprop", "BlanDprop")])
print(GO[which(GO$BlanDprop<30),c("GO", "Name", "Class", "HsapDprop", "BlanDprop")])
write.table(GO[order(GO$BlanDprop),], file = paste(ResultsFolder, "/GOterms_PropData.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)

xlim=c(0,100)
ylim=c(0,100)
ClassOrderDens <- unique(GO$Class)
ClassOrderDens <- c(ClassOrderDens[length(ClassOrderDens)], ClassOrderDens[c(1:(length(ClassOrderDens)-1))])
layout(matrix(c(4,4,7,7,1,1,3,3,6,6,5,2,3,3,6,6,5,2),nrow=3,ncol=6,byrow=T), widths=c(1), heights=c(1), TRUE)

par(mar=c(7,7,1,1),oma=c(1,1,1,1), yaxs='i', xaxs='i')
plot.new()
plot.new()
ScatterPercentagePlot(GO$HsapDprop, GO$BlanDprop, GO$ClassColors, "% of small scale duplicated genes\nH. sapiens", "B. lanceolatum\n% of duplicated genes", xlim, ylim)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,.1), xlim=xlim, col=NA)
mtext("Density", side = 2, line = 5, cex=1.5)
for(go in c(1:length(ClassOrderDens))){
	d <- density(GO$HsapDprop[which(GO$Class==ClassOrderDens[go])], adjust=2)
	polygon(c(0,d$x), c(0,d$y), col=modif_alpha(unique(GO$ClassColors[which(GO$Class==ClassOrderDens[go])]),.5), border=unique(GO$ClassColors[which(GO$Class==ClassOrderDens[go])]))
}
axis(2, at = seq(0, 1, .05), lwd.ticks=1, las=1, cex.axis=1.5)
axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=ylim, xlim=c(0,.25), col=NA)
mtext("Density", side = 3, line = 4, cex=1.5)
for(go in c(1:length(ClassOrderDens))){
	d <- density(GO$BlanDprop[which(GO$Class==ClassOrderDens[go])], adjust=2)
	polygon(c(0,d$y), c(0,d$x), col=modif_alpha(unique(GO$ClassColors[which(GO$Class==ClassOrderDens[go])]),.5), border=unique(GO$ClassColors[which(GO$Class==ClassOrderDens[go])]))
}
axis(3, at = seq(0, .2, .1), lwd.ticks=1, las=1, cex.axis=1.5)
axis(2, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)

xlim=c(0,30)
ylim=c(0,100)
ScatterPercentagePlot(GO$HsapOprop, GO$BlanDprop, GO$ClassColors, "% of ohnologs\nH. sapiens", "", xlim, ylim)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,.3), xlim=xlim, col=NA)
mtext("Density", side = 2, line = 5, cex=1.5)
for(go in c(1:length(ClassOrderDens))){
	d <- density(GO$HsapOprop[which(GO$Class==ClassOrderDens[go])], adjust=2)
	polygon(c(0,d$x), c(0,d$y), col=modif_alpha(unique(GO$ClassColors[which(GO$Class==ClassOrderDens[go])]),.5), border=unique(GO$ClassColors[which(GO$Class==ClassOrderDens[go])]))
}
axis(2, at = seq(0, 1, .1), lwd.ticks=1, las=1, cex.axis=1.5)
axis(1, at = seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/5), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)





plot.new()
box()
plot.new()
plot.new()
box()
plot.new()
box()




quit()

HsapD.GO <- c()
HsapO.GO <- c()
HsapT.GO <- c()
BlanD.GO <- c()
BlanT.GO <- c()
for(g in c(1:length(mGO[,1]))){
	HGenelist <- read.table(paste(GOlistsFolder, "/", mGO$GO[g], ".txt", sep=""), h=F)[,1]
	HsapD.GO <- c(HsapD.GO, length(HsapGenes[which(HsapGenes$DupOrNot=="D" & HsapGenes$HGene %in% HGenelist),1]))
	HsapO.GO <- c(HsapO.GO, length(HsapGenes[which(HsapGenes$DupOrNot=="O" & HsapGenes$HGene %in% HGenelist),1]))
	HsapT.GO <- c(HsapT.GO, length(HsapGenes[which(HsapGenes$HGene %in% HGenelist),1]))
	BGenelist <- OG2BlanGene$BGene[which(OG2BlanGene$OG %in% OG2HsapGene$OG[which(OG2HsapGene$HGene %in% HGenelist)])]
	BlanD.GO <- c(BlanD.GO, length(BlanGenes[which(BlanGenes$DupOrNot=="D" & BlanGenes$BGene %in% BGenelist),1]))
	BlanT.GO <- c(BlanT.GO, length(BlanGenes[which(BlanGenes$BGene %in% BGenelist),1]))
}
mGO$HsapD <- HsapD.GO
mGO$HsapO <- HsapO.GO
mGO$HsapT <- HsapT.GO
mGO$HsapDprop <- HsapD.GO/HsapT.GO*100
mGO$HsapOprop <- HsapO.GO/HsapT.GO*100
mGO$HsapDOprop <- (HsapD.GO+HsapO.GO)/HsapT.GO*100
mGO$BlanD <- BlanD.GO
mGO$BlanT <- BlanT.GO
mGO$BlanDprop <- BlanD.GO/BlanT.GO*100
print(mGO[,c("HsapDOprop","BlanDprop")])

ScatterPercentagePlot(mGO$HsapDprop, mGO$BlanDprop, mGO$ClassColors, "% of small scale duplicated genes in H. sapiens", "% of duplicated genes in B. lanceolatum", c(0,100), c(0,100))
ScatterPercentagePlot(mGO$HsapOprop, mGO$BlanDprop, mGO$ClassColors, "% of ohnolog genes in H. sapiens", "% of duplicated genes in B. lanceolatum", c(0,30), c(0,100))
ScatterPercentagePlot(mGO$HsapDOprop, mGO$BlanDprop, mGO$ClassColors, "% of duplicated + ohnolog genes in H. sapiens", "% of duplicated genes in B. lanceolatum", c(0,100), c(0,100))
ScatterPercentagePlot(mGO$HsapOprop, mGO$HsapDprop, mGO$ClassColors, "% of ohnolog genes in H. sapiens", "% of duplicated genes in H. sapiens", c(0,30), c(0,100))


plot.new()
legend("bottomright", unique(mGO$Class), pch=15, col=unique(mGO$ClassColors), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)
plot.new()
legend("bottomright", mGO$Name, pch=15, col=mGO$ClassColors, bty = "n", pt.cex=2, cex=1, xjust = 0, yjust = 0)




layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(3), heights=c(1), TRUE)

plot.new()
legend("center", GO$Name, pch=15, col=GO$ClassColors, bty = "n", pt.cex=2, cex=1, xjust = 0, yjust = 0)











dev.off()






