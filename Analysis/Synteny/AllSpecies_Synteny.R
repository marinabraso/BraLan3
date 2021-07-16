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

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan", "Bflo")
#ShortSpeciesNames <- c("Blan", "Bflo")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum", "Branchiostoma_floridae")
SpeciesGRefs <- c("Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Homo_sapiens.GRCh38", "Mus_musculus.GRCm39", "Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK")
colfunc <- colorRampPalette(c("royalblue4", "deepskyblue1"))
SpeciesColors <- c(colfunc(4), "orangered4", "orangered1")

TandemMaxDist <- 10000
window <- 1000000
step <- 100000

Dens <- c()
ini <- rep(NA_real_, length(ShortSpeciesNames))
OGTypes <- data.frame("Inter"=ini, "Intra"=ini, "SingleCopy"=ini, "Tandem"=ini)

###########################################################################
# Read data
# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]

# Orthologous groups counts (amphioxus + vertebrates)
CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=NULL)
system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]*_\\([a-z][a-z][a-z]\\)[A-Za-z\\.0-9_]*/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAV) <- c("OG", as.character(unlist(lapply(Header[1,], as.character))))
CountsAV$Ohnologs <- rep(FALSE, length(CountsAV[,1]))
CountsAV$Ohnologs[which(CountsAV$OG %in% Ohnologs)] <- rep(TRUE, length(CountsAV[which(CountsAV$OG %in% Ohnologs),1]))

for(Species in ShortSpeciesNames){
	# Read or produce Species dataframes
	if(file.exists(paste0(ResultsFolder, "/ChrLen_", Species, ".txt")) & file.exists(paste0(ResultsFolder, "/GeneData_", Species, ".txt"))){
		# Chr length information
		ChrLen <- read.delim(paste0(ResultsFolder, "/ChrLen_", Species, ".txt"), h=T, stringsAsFactors=F)
		# Per gene data 
		GeneData <- read.delim(paste0(ResultsFolder, "/GeneData_", Species, ".txt"), h=T, stringsAsFactors=F)
	}else{
		# Read gene location information Species 1
		system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species)], ".gtf.gz | sed 's/\\(.*\\)\\t.*\\t.*\\t\\(.*\\)\\t\\(.*\\)\\t.*\\t.*\\t.*\\t.*gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1\\t\\2\\t\\3\\t\\4/g' | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"), intern=T)
		GeneCoord <- read.table(text=system_out, h=F, sep = "\t")
		colnames(GeneCoord) <- c("Gene","Chr","Start","End","Index")
		GeneCoord$Midp <- GeneCoord$Start + (GeneCoord$End - GeneCoord$Start)/2
		print(head(GeneCoord))

		# OG to Sp1 gene
		system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($", which(colnames(CountsAV)==Species), ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
		OG2Sp1Gene <- read.table(text=system_out, h=F, sep = "\t")
		colnames(OG2Sp1Gene) <- c("OG", "Gene")
		print(head(OG2Sp1Gene))
		GeneData <- merge(OG2Sp1Gene, GeneCoord, by="Gene")
		print(head(GeneData))
		 
		# Read chomosomal length information Species 1
		system_out <- system(paste0("cat ", ChrLenFolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species)], "/", SpeciesGRefs[which(ShortSpeciesNames==Species)], "_lengths.txt | sed 's/>//g' | sed 's/chr//g' | grep -P '^[0-9XY]+'"), intern=T)
		ChrLen <- read.table(text=system_out, h=F, sep = "\t")
		colnames(ChrLen) <- c("Chr","Length","Num","CummLength")
		ChrLen <- ChrLen[order(ChrLen$Chr),]
		numericchr <- as.numeric(ChrLen$Chr[grep("[0-9]", ChrLen$Chr)])
		XYchr <- ChrLen$Chr[grep("[0-9]", ChrLen$Chr, invert=TRUE)]
		ChrLen <- ChrLen[match(c(numericchr[order(numericchr)],XYchr), ChrLen$Chr),]
		ChrLen$NumGenes <- unlist(lapply(ChrLen$Chr, function(x){max(c(0,GeneCoord$Index[which(GeneCoord$Chr==x)]))}))
		ChrLen$CummNumGenes <- cumsum(as.numeric(ChrLen$NumGenes))
		ChrLen$CummLength <- cumsum(as.numeric(ChrLen$Length))

		## Calc num chromosomes, max distance between pairs of copies and num genes between copies
		CountsAV$NumChrs <- unlist(lapply(CountsAV$OG, function(x){length(unique(GeneData$Chr[which(GeneData$OG==x)]))}))
		CountsAV$MaxDist <- unlist(lapply(CountsAV$OG, CalcMaxDist, GeneData=GeneData))
		CountsAV$BetwGenes <- unlist(lapply(CountsAV$OG, CalcNumBetwGenes, GeneData=GeneData))
		CountsAV$DupType <- rep("Missing", length(CountsAV[,1]))
		CountsAV$DupType[which(CountsAV[,Species]==1)] <- rep("SingleCopy",sum(CountsAV[,Species]==1))
		CountsAV$DupType[which(CountsAV$NumChrs==1 & CountsAV[,Species]>1)] <- rep("Intra",sum(CountsAV$NumChrs==1 & CountsAV[,Species]>1))
		CountsAV$DupType[which(CountsAV$NumChrs>1 & CountsAV[,Species]>1)] <- rep("Inter",sum(CountsAV$NumChrs>1 & CountsAV[,Species]>1))
		CountsAV$DupType[which(CountsAV$NumChrs==1 & CountsAV$BetwGenes==0 & CountsAV$MaxDist<=TandemMaxDist & CountsAV[,Species]>1)] <- rep("Tandem",sum(CountsAV$NumChrs==1 & CountsAV$BetwGenes==0 & CountsAV$MaxDist<=TandemMaxDist & CountsAV[,Species]>1))

		# Dup type to Gene Data
		GeneData$DupType <- CountsAV$DupType[match(GeneData$OG, CountsAV$OG)]

		# Write files
		cat(paste0("Write files ", Species, "\n"))
		write.table(ChrLen, file = paste0(ResultsFolder, "/ChrLen_", Species, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
		write.table(GeneData, file = paste0(ResultsFolder, "/GeneData_", Species, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	}
	# Restrict analysis to genes in chomosomal sequences 
	GeneData <- GeneData[which(GeneData$Chr %in% ChrLen$Chr & !is.na(GeneData[,1])),]
	# Sort data
	GeneData <- GeneData[order(GeneData$Chr, GeneData$Midp),]
	print(head(ChrLen))
	print(head(GeneData))
	# Subsets of Gene Data
	GeneDataSC <- GeneData[which(GeneData$DupType=="SingleCopy"),]
	GeneDataD <- GeneData[which(GeneData$DupType!="SingleCopy"),]
	GeneDataIe <- GeneData[which(GeneData$DupType=="Inter"),]
	GeneDataIa <- GeneData[which(GeneData$DupType=="Intra"),]
	GeneDataT <- GeneData[which(GeneData$DupType=="Tandem"),]
	spDens <- c()
	spDens.D <- c()
	spDens.SC <- c()
	spDens.Ie <- c()
	spDens.Ia <- c()
	spDens.T <- c()
	for(chr in ChrLen$Chr){
		chrGeneData <- GeneData[which(GeneData$Chr==chr),]
		spDens <- c(spDens, na.omit(unlist(lapply(seq(0,ChrLen$Length[which(ChrLen$Chr==chr)],step), function(x){length(unique(chrGeneData$Gene[which(chrGeneData$Midp>=x & chrGeneData$Midp<(x+window))]))}))))

		chrGeneDataD <- GeneDataD[which(GeneDataD$Chr==chr),]
		spDens.D <- c(spDens.D, na.omit(unlist(lapply(seq(0,ChrLen$Length[which(ChrLen$Chr==chr)],step), function(x){length(unique(chrGeneDataD$Gene[which(chrGeneDataD$Midp>=x & chrGeneDataD$Midp<(x+window))]))}))))

		chrGeneDataSC <- GeneDataSC[which(GeneDataSC$Chr==chr),]
		spDens.SC <- c(spDens.SC, na.omit(unlist(lapply(seq(0,ChrLen$Length[which(ChrLen$Chr==chr)],step), function(x){length(unique(chrGeneDataSC$Gene[which(chrGeneDataSC$Midp>=x & chrGeneDataSC$Midp<(x+window))]))}))))

		chrGeneDataIe <- GeneDataIe[which(GeneDataIe$Chr==chr),]
		spDens.Ie <- c(spDens.Ie, na.omit(unlist(lapply(seq(0,ChrLen$Length[which(ChrLen$Chr==chr)],step), function(x){length(unique(chrGeneDataIe$Gene[which(chrGeneDataIe$Midp>=x & chrGeneDataIe$Midp<(x+window))]))}))))

		chrGeneDataIa <- GeneDataIa[which(GeneDataIa$Chr==chr),]
		spDens.Ia <- c(spDens.Ia, na.omit(unlist(lapply(seq(0,ChrLen$Length[which(ChrLen$Chr==chr)],step), function(x){length(unique(chrGeneDataIa$Gene[which(chrGeneDataIa$Midp>=x & chrGeneDataIa$Midp<(x+window))]))}))))

		chrGeneDataT <- GeneDataT[which(GeneDataT$Chr==chr),]
		spDens.T <- c(spDens.T, na.omit(unlist(lapply(seq(0,ChrLen$Length[which(ChrLen$Chr==chr)],step), function(x){length(unique(chrGeneDataT$Gene[which(chrGeneDataT$Midp>=x & chrGeneDataT$Midp<(x+window))]))}))))
	}
	names(spDens) <- rep(paste0("G", Species), length(spDens))
	names(spDens.D) <- rep(paste0("D", Species), length(spDens.D))
	names(spDens.SC) <- rep(paste0("SC", Species), length(spDens.SC))
	names(spDens.Ie) <- rep(paste0("Ie", Species), length(spDens.Ie))
	names(spDens.Ia) <- rep(paste0("Ia", Species), length(spDens.Ia))
	names(spDens.T) <- rep(paste0("T", Species), length(spDens.T))
	Dens <- c(Dens, spDens, spDens.D, spDens.SC, spDens.Ie, spDens.Ia, spDens.T)
	OGType <- unique(GeneData[,c("OG", "DupType")])
	OGTypes[which(ShortSpeciesNames==Species),] <- table(OGType$DupType)
}
OGTypes$Total <- rowSums(OGTypes)
OGTypes$TotalD <- rowSums(OGTypes[,c("Inter","Intra","Tandem")])
OGTypes$Species <- ShortSpeciesNames
print(OGTypes)
print(summary(Dens))



pdf(paste0(ResultsFolder, "/AllSpecies_Synteny.pdf"), width=20, height=20)
par(mar=c(10,10,3,3),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4,byrow=T), widths=c(1), heights=c(1), TRUE)

breaks <- seq(0,200,2)
plotDensAllSpecies(Dens, "SC", 3, c(0,.20), c(0,60), "Gene density", "Single copy genes", ShortSpeciesNames, SpeciesColors)
plotDensAllSpecies(Dens, "D", 3, c(0,.20), c(0,60), "Gene density", "Duplicated genes", ShortSpeciesNames, SpeciesColors)

plotVerticalDensAllSpecies(Dens, "SC", 3, c(0,60), paste0("Gene density in ", window/1000, "kbp windows"), "Single copy genes", ShortSpeciesNames, SpeciesColors)
plotVerticalDensAllSpecies(Dens, "D", 3, c(0,60), paste0("Gene density in ", window/1000, "kbp windows"), "Duplicated genes", ShortSpeciesNames, SpeciesColors)
plotVerticalDensAllSpecies(Dens, "G", 3, c(0,120), paste0("Gene density in ", window/1000, "kbp windows"), "All genes", ShortSpeciesNames, SpeciesColors)
plotVerticalDensAllSpeciesSCD(Dens, 3, c(0,80), paste0("Gene density in ", window/1000, "kbp windows"), "All genes", ShortSpeciesNames, SpeciesColors)

plotOGTypePerSpecies(OGTypes)



dev.off()













