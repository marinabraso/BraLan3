#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
#library(viridis)

######################################################################
# Files & folders

ResultsFolder <- "Plots/Synteny"
CountsFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
OhnologsSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"
GTFfolder <- "Results/FilteringGeneSets/GTFs"
ChrLenFolder <- "Results/AssemblyStatistics"

######################################################################
# General parameters
Species1 <- "Blan"
Species2 <- "Bflo"

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan", "Bflo")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum", "Branchiostoma_floridae")
SpeciesGRefs <- c("Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Homo_sapiens.GRCh38", "Mus_musculus.GRCm39", "Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK")

print(paste0("Species1: ", Species1, " ", SpeciesNames[which(ShortSpeciesNames==Species1)], " ", SpeciesGRefs[which(ShortSpeciesNames==Species1)]))
print(paste0("Species2: ", Species2, " ", SpeciesNames[which(ShortSpeciesNames==Species2)], " ", SpeciesGRefs[which(ShortSpeciesNames==Species2)]))

TandemMaxDist <- 10000

###########################################################################
# Read data
# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]

# Orthologous groups counts (amphioxus + vertebrates)
CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]*_\\([a-z][a-z][a-z]\\)[A-Za-z\\.0-9_]*/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAV) <- as.character(unlist(lapply(Header[1,], as.character)))
CountsAV$Ohnologs <- rep(FALSE, length(CountsAV[,1]))
CountsAV$Ohnologs[which(CountsAV$OG %in% Ohnologs)] <- rep(TRUE, length(CountsAV[which(CountsAV$OG %in% Ohnologs),1]))
print(head(CountsAV))

# Read gene location information Species 1
system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], ".gtf.gz | sed 's/gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1/g' | cut -f1,4,5,9 | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"), intern=T)
GeneCoord1 <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord1) <- c("Gene","Chr","Start","End","Index")
GeneCoord1$Midp <- GeneCoord1$Start + (GeneCoord1$End - GeneCoord1$Start)/2
print(head(GeneCoord1))

# Read gene location information Species 1
print(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], ".gtf.gz | sed 's/gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1/g' | cut -f1,4,5,9 | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"))
quit()
system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], ".gtf.gz | sed 's/gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1/g' | cut -f1,4,5,9 | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"), intern=T)
GeneCoord2 <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord2) <- c("Gene","Chr","Start","End","Index")
GeneCoord2$Midp <- GeneCoord2$Start + (GeneCoord2$End - GeneCoord2$Start)/2
print(head(GeneCoord2))

# Read chomosomal length information Species 1
system_out <- system(paste0("cat ", ChrLenFolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], "_lengths.txt | sed 's/>//g' | sed 's/chr//g' | grep -P '^[0-9XY]+'"), intern=T)
ChrLen1 <- read.table(text=system_out, h=F, sep = "\t")
colnames(ChrLen1) <- c("Chr","Length","Num","CummLength")
ChrLen1 <- ChrLen1[order(ChrLen1$Chr),]
numericchr <- as.numeric(ChrLen1$Chr[grep("[0-9]", ChrLen1$Chr)])
XYchr <- ChrLen1$Chr[grep("[0-9]", ChrLen1$Chr, invert=TRUE)]
ChrLen1 <- ChrLen1[match(c(numericchr[order(numericchr)],XYchr), ChrLen1$Chr),]
ChrLen1$NumGenes <- unlist(lapply(ChrLen1$Chr, function(x){max(c(0,GeneCoord1$Index[which(GeneCoord1$Chr==x)]))}))
ChrLen1$CummNumGenes <- cumsum(as.numeric(ChrLen1$NumGenes))
ChrLen1$CummLength <- cumsum(as.numeric(ChrLen1$Length))
print(ChrLen1)

# Read chomosomal length information Species 2
system_out <- system(paste0("cat ", ChrLenFolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "_lengths.txt | sed 's/>//g' | sed 's/chr//g' | grep -P '^[0-9XY]+'"), intern=T)
ChrLen2 <- read.table(text=system_out, h=F, sep = "\t")
colnames(ChrLen2) <- c("Chr","Length","Num","CummLength")
numericchr <- as.numeric(ChrLen2$Chr[grep("[0-9]", ChrLen2$Chr)])
XYchr <- ChrLen2$Chr[grep("[0-9]", ChrLen2$Chr, invert=TRUE)]
ChrLen2 <- ChrLen2[match(c(numericchr[order(numericchr)],XYchr), ChrLen2$Chr),]
ChrLen2$NumGenes <- unlist(lapply(ChrLen2$Chr, function(x){max(c(0,GeneCoord2$Index[which(GeneCoord2$Chr==x)]))}))
ChrLen2$CummNumGenes <- cumsum(as.numeric(ChrLen2$NumGenes))
ChrLen2$CummLength <- cumsum(as.numeric(ChrLen2$Length))
print(ChrLen2)

# OG to Sp1 gene
system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($", which(colnames(CountsAV)==Species1)+1, ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
OG2Sp1Gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2Sp1Gene) <- c("OG", "Gene")
print(head(OG2Sp1Gene))

# OG to Sp2 gene
system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($", which(colnames(CountsAV)==Species2)+1, ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
OG2Sp2Gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2Sp2Gene) <- c("OG", "Gene")
print(head(OG2Sp2Gene))

## Calc num chromosomes, max distance between pairs of copies and num genes between copies
GeneData1 <- merge(OG2Sp1Gene, GeneCoord1, by="Gene")
CountsAV$NumChrs1 <- unlist(lapply(rownames(CountsAV), function(x){length(unique(GeneData1$Chr[which(GeneData1$OG==x)]))}))
CountsAV$MaxDist1 <- unlist(lapply(rownames(CountsAV), CalcMaxDist, GeneData=GeneData1))
CountsAV$BetwGenes1 <- unlist(lapply(rownames(CountsAV), CalcNumBetwGenes, GeneData=GeneData1))
CountsAV$DupType1 <- rep(NA, length(CountsAV[,1]))
CountsAV$DupType1[which(CountsAV$NumChrs1==1 & CountsAV[,Species1]>1)] <- rep("Intra",sum(CountsAV$NumChrs1==1 & CountsAV[,Species1]>1))
CountsAV$DupType1[which(CountsAV$NumChrs1>1 & CountsAV[,Species1]>1)] <- rep("Inter",sum(CountsAV$NumChrs1>1 & CountsAV[,Species1]>1))
CountsAV$DupType1[which(CountsAV$NumChrs1==1 & CountsAV$BetwGenes1==0 & CountsAV$MaxDist1<=TandemMaxDist & CountsAV[,Species1]>1)] <- rep("Tandem",sum(CountsAV$NumChrs1==1 & CountsAV$BetwGenes1==0 & CountsAV$MaxDist1<=TandemMaxDist & CountsAV[,Species1]>1))

GeneData2 <- merge(OG2Sp2Gene, GeneCoord2, by="Gene")
CountsAV$NumChrs2 <- unlist(lapply(rownames(CountsAV), function(x){length(unique(GeneData2$Chr[which(GeneData2$OG==x)]))}))
CountsAV$MaxDist2 <- unlist(lapply(rownames(CountsAV), CalcMaxDist, GeneData=GeneData2))
CountsAV$BetwGenes2 <- unlist(lapply(rownames(CountsAV), CalcNumBetwGenes, GeneData=GeneData2))
CountsAV$DupType2 <- rep(NA, length(CountsAV[,1]))
CountsAV$DupType2[which(CountsAV$NumChrs2==1 & CountsAV[,Species2]>1)] <- rep("Intra",sum(CountsAV$NumChrs2==1 & CountsAV[,Species2]>1))
CountsAV$DupType2[which(CountsAV$NumChrs2>1 & CountsAV[,Species2]>1)] <- rep("Inter",sum(CountsAV$NumChrs2>1 & CountsAV[,Species2]>1))
CountsAV$DupType2[which(CountsAV$NumChrs2==1 & CountsAV$BetwGenes2==0 & CountsAV$MaxDist2<=TandemMaxDist & CountsAV[,Species2]>1)] <- rep("Tandem",sum(CountsAV$NumChrs2==1 & CountsAV$BetwGenes2==0 & CountsAV$MaxDist2<=TandemMaxDist & CountsAV[,Species2]>1))

table(CountsAV$DupType1)
table(CountsAV[,Species1]>0)
table(CountsAV[,Species1]>1)
table(CountsAV$DupType2)
table(CountsAV[,Species2]>0)
table(CountsAV[,Species2]>1)
print(head(CountsAV[which(!is.na(CountsAV$DupType1) & !is.na(CountsAV$DupType2)),]))

GeneData1 <- as.matrix(GeneData1)
GeneData2 <- as.matrix(GeneData2)
tmpGenePairs1 <- unlist(lapply(GeneData1[,"Gene"], RetrieveGenesInOtherSpecies, GeneData=GeneData1, otherGeneData=GeneData2))
tmpGenePairs2 <- unlist(lapply(GeneData2[,"Gene"], RetrieveGenesInOtherSpecies, GeneData=GeneData2, otherGeneData=GeneData1))
GeneData1 <- as.data.frame(GeneData1)
GeneData2 <- as.data.frame(GeneData2)
splitGenePairs1 <- do.call(rbind, strsplit(tmpGenePairs1,'_'))
splitGenePairs2 <- do.call(rbind, strsplit(tmpGenePairs2,'_'))
GenePairs <- as.data.frame(unique(rbind(splitGenePairs1, splitGenePairs2[,c(2,1)])))
colnames(GenePairs) <- c(Species1, Species2)

GenePairs$Chr1 <- GeneData1$Chr[match(GenePairs[,Species1], GeneData1$Gene)]
GenePairs$Midp1 <- as.numeric(GeneData1$Midp[match(GenePairs[,Species1], GeneData1$Gene)])
GenePairs$Index1 <- as.numeric(GeneData1$Index[match(GenePairs[,Species1], GeneData1$Gene)])

GenePairs$Chr2 <- GeneData2$Chr[match(GenePairs[,Species2], GeneData2$Gene)]
GenePairs$Midp2 <- as.numeric(GeneData2$Midp[match(GenePairs[,Species2], GeneData2$Gene)])
GenePairs$Index2 <- as.numeric(GeneData2$Index[match(GenePairs[,Species2], GeneData2$Gene)])
print(head(GenePairs))



print(head(ChrLen1$CummLength[match(GenePairs$Chr1, ChrLen1$Chr)]))
print(head(match(GenePairs$Chr1, ChrLen1$Chr)))

GenePairs <- GenePairs[which(GenePairs$Chr1 %in% ChrLen1$Chr & GenePairs$Chr2 %in% ChrLen2$Chr),]

pdf(paste0(ResultsFolder, "/", Species1, Species2, "_Synteny.pdf"), width=20, height=20)
par(mar=c(10,10,3,3),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(1), heights=c(1), TRUE)


plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,ChrLen2$CummLength[length(ChrLen2[,1])]), xlim=c(0,ChrLen1$CummLength[length(ChrLen1[,1])]), col=NA)
mtext(SpeciesNames[which(ShortSpeciesNames==Species1)], side = 1, line = 5, cex=1.5)
mtext(SpeciesNames[which(ShortSpeciesNames==Species2)], side = 2, line = 5, cex=1.5)
Pos1 <- GenePairs$Midp1+ChrLen1$CummLength[match(GenePairs$Chr1, ChrLen1$Chr)]-ChrLen1$Length[match(GenePairs$Chr1, ChrLen1$Chr)]
Pos2 <- GenePairs$Midp2+ChrLen2$CummLength[match(GenePairs$Chr2, ChrLen2$Chr)]-ChrLen2$Length[match(GenePairs$Chr2, ChrLen2$Chr)]
points(Pos1, Pos2, col="darkred", pch=16, cex=.5)
axis(1, at = c(0,ChrLen1$CummLength), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
axis(2, at = c(0,ChrLen2$CummLength), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
axis(1, at = ChrLen1$CummLength-ChrLen1$Length/2, labels=ChrLen1$Chr, lwd.ticks=NA, las=1, cex.axis=1.5)
axis(2, at = ChrLen2$CummLength-ChrLen2$Length/2, labels=ChrLen2$Chr, lwd.ticks=NA, las=1, cex.axis=1.5)
box()

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,ChrLen2$CummNumGenes[length(ChrLen2[,1])]), xlim=c(0,ChrLen1$CummNumGenes[length(ChrLen1[,1])]), col=NA)
mtext(SpeciesNames[which(ShortSpeciesNames==Species1)], side = 1, line = 5, cex=1.5)
mtext(SpeciesNames[which(ShortSpeciesNames==Species2)], side = 2, line = 5, cex=1.5)
abline(v=ChrLen1$CummNumGenes, lty=2, lwd=.5, col="grey60")
abline(h=ChrLen2$CummNumGenes, lty=2, lwd=.5, col="grey60")
Pos1 <- GenePairs$Index1+ChrLen1$CummNumGenes[match(GenePairs$Chr1, ChrLen1$Chr)]-ChrLen1$NumGenes[match(GenePairs$Chr1, ChrLen1$Chr)]
Pos2 <- GenePairs$Index2+ChrLen2$CummNumGenes[match(GenePairs$Chr2, ChrLen2$Chr)]-ChrLen2$NumGenes[match(GenePairs$Chr2, ChrLen2$Chr)]
points(Pos1, Pos2, col="darkred", pch=16, cex=.5)
axis(1, at = c(0,ChrLen1$CummNumGenes), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
axis(2, at = c(0,ChrLen2$CummNumGenes), labels=NA, lwd.ticks=1, las=1, cex.axis=1.5)
axis(1, at = ChrLen1$CummNumGenes-ChrLen1$NumGenes/2, labels=ChrLen1$Chr, lwd.ticks=NA, las=1, cex.axis=1.5)
axis(2, at = ChrLen2$CummNumGenes-ChrLen2$NumGenes/2, labels=ChrLen2$Chr, lwd.ticks=NA, las=1, cex.axis=1.5)
box()














dev.off()