#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
#library(viridis)
library(qvalue)

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
PQvalThreshold <- 0.01 # for hypergeometric test
Species1 <- "Blan"
Species2 <- "Hsap"

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

if(file.exists(paste0(ResultsFolder, "/CountsAV_", Species1, "_", Species2, ".txt"))){
	# Orthologous groups information
	CountsAV <- read.delim(paste0(ResultsFolder, "/CountsAV_", Species1, "_", Species2, ".txt"), h=T, stringsAsFactors=F)
}else{
	# Orthologous groups counts (amphioxus + vertebrates)
	CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=NULL)
	system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]*_\\([a-z][a-z][a-z]\\)[A-Za-z\\.0-9_]*/\\1\\2/g'"), intern=T)
	Header <- read.table(text=system_out, h=F, sep = "\t")
	colnames(CountsAV) <- c("OG", as.character(unlist(lapply(Header[1,], as.character))))
	CountsAV$Ohnologs <- rep(FALSE, length(CountsAV[,1]))
	CountsAV$Ohnologs[which(CountsAV$OG %in% Ohnologs)] <- rep(TRUE, length(CountsAV[which(CountsAV$OG %in% Ohnologs),1]))
}

# Read or produce Species1 dataframes
if(file.exists(paste0(ResultsFolder, "/ChrLen_", Species1, ".txt")) & file.exists(paste0(ResultsFolder, "/GeneData_", Species1, ".txt"))){
	# Chr length information
	ChrLen1 <- read.delim(paste0(ResultsFolder, "/ChrLen_", Species1, ".txt"), h=T, stringsAsFactors=F)
	# Per gene data 
	GeneData1 <- read.delim(paste0(ResultsFolder, "/GeneData_", Species1, ".txt"), h=T, stringsAsFactors=F)
}else{
	# Read gene location information Species 1
	system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], ".gtf.gz | sed 's/\\(.*\\)\\t.*\\t.*\\t\\(.*\\)\\t\\(.*\\)\\t.*\\t.*\\t.*\\t.*gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1\\t\\2\\t\\3\\t\\4/g' | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"), intern=T)
	GeneCoord1 <- read.table(text=system_out, h=F, sep = "\t")
	colnames(GeneCoord1) <- c("Gene","Chr","Start","End","Index")
	GeneCoord1$Midp <- GeneCoord1$Start + (GeneCoord1$End - GeneCoord1$Start)/2

	# OG to Sp1 gene
	system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($", which(colnames(CountsAV)==Species1), ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
	OG2Sp1Gene <- read.table(text=system_out, h=F, sep = "\t")
	colnames(OG2Sp1Gene) <- c("OG", "Gene")
	print(head(OG2Sp1Gene))
	GeneData1 <- merge(OG2Sp1Gene, GeneCoord1, by="Gene")
	 
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

	# Write files
	cat(paste0("Write files ", Species1, "\n"))
	write.table(ChrLen1, file = paste0(ResultsFolder, "/ChrLen_", Species1, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	write.table(GeneData1, file = paste0(ResultsFolder, "/GeneData_", Species1, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
}
print(head(ChrLen1))
print(head(GeneData1))

# Read or produce Species2 dataframes
if(file.exists(paste0(ResultsFolder, "/ChrLen_", Species2, ".txt")) & file.exists(paste0(ResultsFolder, "/GeneData_", Species2, ".txt"))){
	# Chr length information
	ChrLen2 <- read.delim(paste0(ResultsFolder, "/ChrLen_", Species2, ".txt"), h=T, stringsAsFactors=F)
	# Per gene data 
	GeneData2 <- read.delim(paste0(ResultsFolder, "/GeneData_", Species2, ".txt"), h=T, stringsAsFactors=F)
}else{
	# Read gene location information Species 2
	system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], ".gtf.gz | sed 's/\\(.*\\)\\t.*\\t.*\\t\\(.*\\)\\t\\(.*\\)\\t.*\\t.*\\t.*\\t.*gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1\\t\\2\\t\\3\\t\\4/g' | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"), intern=T)
	GeneCoord2 <- read.table(text=system_out, h=F, sep = "\t")
	colnames(GeneCoord2) <- c("Gene","Chr","Start","End","Index")
	GeneCoord2$Midp <- GeneCoord2$Start + (GeneCoord2$End - GeneCoord2$Start)/2

	# OG to Sp2 gene
	system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($", which(colnames(CountsAV)==Species2), ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
	OG2Sp2Gene <- read.table(text=system_out, h=F, sep = "\t")
	colnames(OG2Sp2Gene) <- c("OG", "Gene")
	print(head(OG2Sp2Gene))
	GeneData2 <- merge(OG2Sp2Gene, GeneCoord2, by="Gene")

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

	# Write files
	cat(paste0("Write files ", Species2, "\n"))
	write.table(ChrLen2, file = paste0(ResultsFolder, "/ChrLen_", Species2, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
	write.table(GeneData2, file = paste0(ResultsFolder, "/GeneData_", Species2, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
}
print(head(ChrLen2))
print(head(GeneData2))

if(!file.exists(paste0(ResultsFolder, "/CountsAV_", Species1, "_", Species2, ".txt"))){
	## Calc num chromosomes, max distance between pairs of copies and num genes between copies
	CountsAV$NumChrs1 <- unlist(lapply(CountsAV$OG, function(x){length(unique(GeneData1$Chr[which(GeneData1$OG==x)]))}))
	CountsAV$MaxDist1 <- unlist(lapply(CountsAV$OG, CalcMaxDist, GeneData=GeneData1))
	CountsAV$BetwGenes1 <- unlist(lapply(CountsAV$OG, CalcNumBetwGenes, GeneData=GeneData1))
	CountsAV$DupType1 <- rep("Missing", length(CountsAV[,1]))
	CountsAV$DupType1[which(CountsAV[,Species1]==1)] <- rep("SingleCopy",sum(CountsAV[,Species1]==1))
	CountsAV$DupType1[which(CountsAV$NumChrs1==1 & CountsAV[,Species1]>1)] <- rep("Intra",sum(CountsAV$NumChrs1==1 & CountsAV[,Species1]>1))
	CountsAV$DupType1[which(CountsAV$NumChrs1>1 & CountsAV[,Species1]>1)] <- rep("Inter",sum(CountsAV$NumChrs1>1 & CountsAV[,Species1]>1))
	CountsAV$DupType1[which(CountsAV$NumChrs1==1 & CountsAV$BetwGenes1==0 & CountsAV$MaxDist1<=TandemMaxDist & CountsAV[,Species1]>1)] <- rep("Tandem",sum(CountsAV$NumChrs1==1 & CountsAV$BetwGenes1==0 & CountsAV$MaxDist1<=TandemMaxDist & CountsAV[,Species1]>1))

	CountsAV$NumChrs2 <- unlist(lapply(CountsAV$OG, function(x){length(unique(GeneData2$Chr[which(GeneData2$OG==x)]))}))
	CountsAV$MaxDist2 <- unlist(lapply(CountsAV$OG, CalcMaxDist, GeneData=GeneData2))
	CountsAV$BetwGenes2 <- unlist(lapply(CountsAV$OG, CalcNumBetwGenes, GeneData=GeneData2))
	CountsAV$DupType2 <- rep("Missing", length(CountsAV[,1]))
	CountsAV$DupType2[which(CountsAV[,Species2]==1)] <- rep("SingleCopy",sum(CountsAV[,Species2]==1))
	CountsAV$DupType2[which(CountsAV$NumChrs2==1 & CountsAV[,Species2]>1)] <- rep("Intra",sum(CountsAV$NumChrs2==1 & CountsAV[,Species2]>1))
	CountsAV$DupType2[which(CountsAV$NumChrs2>1 & CountsAV[,Species2]>1)] <- rep("Inter",sum(CountsAV$NumChrs2>1 & CountsAV[,Species2]>1))
	CountsAV$DupType2[which(CountsAV$NumChrs2==1 & CountsAV$BetwGenes2==0 & CountsAV$MaxDist2<=TandemMaxDist & CountsAV[,Species2]>1)] <- rep("Tandem",sum(CountsAV$NumChrs2==1 & CountsAV$BetwGenes2==0 & CountsAV$MaxDist2<=TandemMaxDist & CountsAV[,Species2]>1))
	write.table(CountsAV, file = paste0(ResultsFolder, "/CountsAV_", Species1, "_", Species2, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
}
print(head(CountsAV[which(!is.na(CountsAV$DupType1) & !is.na(CountsAV$DupType2)),]))

if(file.exists(paste0(ResultsFolder, "/GenePairs_", Species1, "_", Species2, ".txt"))){
	# Gene paris data
	GenePairs <- read.delim(paste0(ResultsFolder, "/GenePairs_", Species1, "_", Species2, ".txt"), h=T, stringsAsFactors=F)
}else{
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
	GenePairs$OG <- GeneData1$OG[match(GenePairs[,Species1], GeneData1$Gene)]
	GenePairs$Chr1 <- GeneData1$Chr[match(GenePairs[,Species1], GeneData1$Gene)]
	GenePairs$Midp1 <- as.numeric(GeneData1$Midp[match(GenePairs[,Species1], GeneData1$Gene)])
	GenePairs$Index1 <- as.numeric(GeneData1$Index[match(GenePairs[,Species1], GeneData1$Gene)])
	GenePairs$DupType1 <- CountsAV$DupType1[match(GenePairs$OG, CountsAV$OG)]
	GenePairs$Chr2 <- GeneData2$Chr[match(GenePairs[,Species2], GeneData2$Gene)]
	GenePairs$Midp2 <- as.numeric(GeneData2$Midp[match(GenePairs[,Species2], GeneData2$Gene)])
	GenePairs$Index2 <- as.numeric(GeneData2$Index[match(GenePairs[,Species2], GeneData2$Gene)])
	GenePairs$DupType2 <- CountsAV$DupType2[match(GenePairs$OG, CountsAV$OG)]
	write.table(GenePairs, file = paste0(ResultsFolder, "/GenePairs_", Species1, "_", Species2, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
}
print(head(GenePairs))

# Printing tables
table(CountsAV$DupType1)
table(CountsAV[,Species1]>0)
table(CountsAV[,Species1]>1)
table(CountsAV$DupType2)
table(CountsAV[,Species2]>0)
table(CountsAV[,Species2]>1)

# Dup type to Gene Data
GeneData1$DupType <- CountsAV$DupType1[match(GeneData1$OG, CountsAV$OG)]
GeneData2$DupType <- CountsAV$DupType2[match(GeneData2$OG, CountsAV$OG)]

# Restrict plots to genes in chomosomal sequences 
GenePairs <- GenePairs[which(GenePairs$Chr1 %in% ChrLen1$Chr & GenePairs$Chr2 %in% ChrLen2$Chr),]
GeneData1 <- GeneData1[which(GeneData1$Chr %in% ChrLen1$Chr),]
GeneData2 <- GeneData2[which(GeneData2$Chr %in% ChrLen2$Chr),]

# Subsets of GenePairs
GenePairsSC <- GenePairs[which(GenePairs$DupType1=="SingleCopy" & GenePairs$DupType2=="SingleCopy"),]
GenePairsD <- GenePairs[which(GenePairs$DupType1!="SingleCopy" & GenePairs$DupType2!="SingleCopy"),]
GenePairsD.Ie <- GenePairsD[which(GenePairsD$DupType1=="Inter" & GenePairsD$DupType1=="Inter"),]
GenePairsD.Ia <- GenePairsD[which(GenePairsD$DupType1=="Intra" & GenePairsD$DupType1=="Intra"),]
GenePairsD.T <- GenePairsD[which(GenePairsD$DupType1=="Tandem" & GenePairsD$DupType1=="Tandem"),]

# Subsets of Gene Data
GeneData1SC <- GeneData1[which(GeneData1$DupType=="SingleCopy"),]
GeneData1D <- GeneData1[which(GeneData1$DupType!="SingleCopy"),]
GeneData1D.Ie <- GeneData1[which(GeneData1$DupType=="Inter"),]
GeneData1D.Ia <- GeneData1[which(GeneData1$DupType=="Intra"),]
GeneData1D.T <- GeneData1[which(GeneData1$DupType=="Tandem"),]
GeneData2SC <- GeneData2[which(GeneData2$DupType=="SingleCopy"),]
GeneData2D <- GeneData2[which(GeneData2$DupType!="SingleCopy"),]
GeneData2D.Ie <- GeneData2[which(GeneData2$DupType=="Inter"),]
GeneData2D.Ia <- GeneData2[which(GeneData2$DupType=="Intra"),]
GeneData2D.T <- GeneData2[which(GeneData2$DupType=="Tandem"),]

pdf(paste0(ResultsFolder, "/", Species1, Species2, "_Synteny.pdf"), width=20, height=20)
par(mar=c(10,10,3,3),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(1), heights=c(1), TRUE)

#SyntenyAlongChrPlot(GenePairs, "Midp", "Gene mid point coordinates")
#SyntenyAlongChrPlot(GenePairs, "Index", "Gene order")
#SyntenyAlongChrPlot(GenePairsSC, "Index", "Gene order: single copy genes")
#SyntenyAlongChrPlot(GenePairsD, "Index", "Gene order: duplicated genes")
#SyntenyAlongChrPlot(GenePairsD.Ie, "Index", "Gene order: interchromosomal duplicated genes")
#SyntenyAlongChrPlot(GenePairsD.Ia, "Index", "Gene order: intrachromosomal duplicated genes")
#SyntenyAlongChrPlot(GenePairsD.T, "Index", "Gene order: tandem duplicated genes")

print(table(GenePairs$DupType1))
print(table(GenePairs$DupType2))
SyntenyAlongChrPlot("Midp", "Gene mid point coordinates", GenePairsSC, GenePairsD)
SyntenyAlongChrPlot("Index", "Gene order", GenePairsSC, GenePairsD)
quit()

GeneData1 <- GeneData1[order(GeneData1$Chr, GeneData1$Midp),]
GeneData2 <- GeneData2[order(GeneData2$Chr, GeneData2$Midp),]

layout(matrix(c(1,2,3,4,5),nrow=5,ncol=1,byrow=T), widths=c(10), heights=c(1), TRUE)
window <- 1000000
step <- 100000
GDensity1 <- c()
DDensity1 <- c()
InterGDistance1 <- c()
InterDDistance1 <- c()
for(chr in ChrLen1$Chr){
	chrGeneData1 <- GeneData1[which(GeneData1$Chr==chr),]
	chrGeneData1D <- GeneData1D[which(GeneData1D$Chr==chr),]
	chrInterGDistance1 <- unlist(lapply(c(1:(length(chrGeneData1[,1])-1)), function(x){chrGeneData1$Midp[x+1]-chrGeneData1$Midp[x]}))
	chrInterDDistance1 <- unlist(lapply(c(1:(length(chrGeneData1D[,1])-1)), function(x){chrGeneData1D$Midp[x+1]-chrGeneData1D$Midp[x]}))
	chrGDensity1 <- unlist(lapply(seq(0,ChrLen1$Length[which(ChrLen1$Chr==chr)],step), function(x){length(unique(GeneData1$Gene[which(GeneData1$Chr==chr & GeneData1$Midp>=x & GeneData1$Midp<(x+window))]))}))
	chrDDensity1 <- unlist(lapply(seq(0,ChrLen1$Length[which(ChrLen1$Chr==chr)],step), function(x){length(unique(GeneData1$Gene[which(GeneData1$DupType!="SingleCopy" & GeneData1$Chr==chr & GeneData1$Midp>=x & GeneData1$Midp<(x+window))]))}))
	GDensity1 <- c(GDensity1, chrGDensity1)
	DDensity1 <- c(DDensity1, chrDDensity1)
	InterGDistance1 <- c(InterGDistance1, chrInterGDistance1)
	InterDDistance1 <- c(InterDDistance1, chrInterDDistance1)
	plot(c(1:10), c(1:10), axes=F, las=1, xlab="", ylab="", ylim=c(0,80), xlim=c(0,max(ChrLen1$Length)), col=NA)
	mtext(chr, side = 2, line = 3, cex=1.5, las=1)
	lines(seq(0,ChrLen1$Length[which(ChrLen1$Chr==chr)],step)+window/2, chrGDensity1, lwd=2, col="royalblue1")
	lines(seq(0,ChrLen1$Length[which(ChrLen1$Chr==chr)],step)+window/2, chrDDensity1, lwd=2, col="royalblue4")
}
layout(matrix(c(1,2,3,4,5),nrow=5,ncol=1,byrow=T), widths=c(10), heights=c(1), TRUE)
GDensity2 <- c()
DDensity2 <- c()
InterGDistance2 <- c()
InterDDistance2 <- c()
for(chr in ChrLen2$Chr){
	chrGeneData2 <- GeneData2[which(GeneData2$Chr==chr),]
	chrGeneData2D <- GeneData2D[which(GeneData2D$Chr==chr),]
	chrInterGDistance2 <- unlist(lapply(c(1:(length(chrGeneData2[,1])-1)), function(x){chrGeneData2$Midp[x+1]-chrGeneData2$Midp[x]}))
	chrInterDDistance2 <- unlist(lapply(c(1:(length(chrGeneData2D[,1])-1)), function(x){chrGeneData2D$Midp[x+1]-chrGeneData2D$Midp[x]}))
	chrGDensity2 <- unlist(lapply(seq(0,ChrLen2$Length[which(ChrLen2$Chr==chr)],step), function(x){length(unique(GeneData2$Gene[which(GeneData2$Chr==chr & GeneData2$Midp>=x & GeneData2$Midp<(x+window))]))}))
	chrDDensity2 <- unlist(lapply(seq(0,ChrLen2$Length[which(ChrLen2$Chr==chr)],step), function(x){length(unique(GeneData2$Gene[which(GeneData2$DupType!="SingleCopy" & GeneData2$Chr==chr & GeneData2$Midp>=x & GeneData2$Midp<(x+window))]))}))
	GDensity2 <- c(GDensity2, chrGDensity2)
	DDensity2 <- c(DDensity2, chrDDensity2)
	InterGDistance2 <- c(InterGDistance2, chrInterGDistance2)
	InterDDistance2 <- c(InterDDistance2, chrInterDDistance2)
	plot(c(1:10), c(1:10), axes=F, las=1, xlab="", ylab="", ylim=c(0,80), xlim=c(0,max(ChrLen2$Length)), col=NA)
	mtext(chr, side = 2, line = 3, cex=1.5, las=1)
	lines(seq(0,ChrLen2$Length[which(ChrLen2$Chr==chr)],step)+window/2, chrGDensity2, lwd=2, col="orangered1")
	lines(seq(0,ChrLen2$Length[which(ChrLen2$Chr==chr)],step)+window/2, chrDDensity2, lwd=2, col="orangered4")
}

layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)


breaks <- seq(0,5000,1)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,log(10000)), xlim=c(0,150), col=NA)
mtext("Absolute frequency", side = 2, line = 3, cex=1.5)
mtext("Gene density", side = 1, line = 3, cex=1.5)
h1 <- hist(GDensity1, breaks=breaks, plot=FALSE)
h2 <- hist(GDensity2, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h1, "royalblue1")
AddHistlogYaxis(h2, "orangered1")
axis(1, at = seq(0, 1000,20), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = c(0,log(c(1,10,100,1000,10000))), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
box()
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,log(10000)), xlim=c(0,150), col=NA)
mtext("Absolute frequency", side = 2, line = 3, cex=1.5)
mtext("Gene density", side = 1, line = 3, cex=1.5)
h1 <- hist(DDensity1, breaks=breaks, plot=FALSE)
h2 <- hist(DDensity2, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h1, "royalblue4")
AddHistlogYaxis(h2, "orangered4")
axis(1, at = seq(0, 1000,20), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = c(0,log(c(1,10,100,1000,10000))), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
box()





breaks <- seq(0,100000000,1000)
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,log(10000)), xlim=c(0,500000), col=NA)
mtext("Absolute frequency", side = 2, line = 3, cex=1.5)
mtext("Intergenic distance", side = 1, line = 3, cex=1.5)
h1 <- hist(InterGDistance1, breaks=breaks, plot=FALSE)
h2 <- hist(InterGDistance2, breaks=breaks, plot=FALSE)
AddHistlogYaxis(h1, "royalblue4")
AddHistlogYaxis(h2, "orangered4")
axis(1, at = seq(0, 1000,20), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = c(0,log(c(1,10,100,1000,10000))), labels=c(0,1,10,100,1000,10000), lwd.ticks=1, las=1, cex.axis=1)
box()

plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,.1), xlim=c(0,100), col=NA)
mtext("Density", side = 2, line = 3, cex=1.5)
mtext("Gene density", side = 1, line = 3, cex=1.5)
lines(density(GDensity1, adjust=5), lwd=2, col="royalblue1")
lines(density(GDensity2, adjust=5), lwd=2, col="orangered1")
lines(density(DDensity1, adjust=5), lwd=2, col="royalblue4")
lines(density(DDensity2, adjust=5), lwd=2, col="orangered4")
axis(1, at = seq(0, 100, 5), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = seq(0, 1, .05), lwd.ticks=1, las=1, cex.axis=1)
box()

print(summary(InterGDistance2))
print(head(density(InterGDistance1)$y))
plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0,.00005), xlim=c(0,500000), col=NA)
mtext("Density", side = 2, line = 3, cex=1.5)
mtext("Gene density", side = 1, line = 3, cex=1.5)
lines(density(InterGDistance1, adjust=1), lwd=2, col="royalblue1")
lines(density(InterGDistance2, adjust=1), lwd=2, col="orangered1")
lines(density(InterDDistance1, adjust=1), lwd=2, col="royalblue4")
lines(density(InterDDistance2, adjust=1), lwd=2, col="orangered4")
axis(1, at = seq(0, 500000, 50000), lwd.ticks=1, las=1, cex.axis=1)
axis(2, at = seq(0, 1, .00001), lwd.ticks=1, las=1, cex.axis=1)
box()



HyperTests <- data.frame()
for(vtype in c("SingleCopy","Inter","Intra", "Tandem")){
	for(atype in c("SingleCopy","Inter","Intra", "Tandem")){
		vnum <- length(GenePairs[which(GenePairs$DupType2 == vtype),1])
		anum <- length(GenePairs[which(GenePairs$DupType1 == atype),1])
		onum <- length(GenePairs[which(GenePairs$DupType2 == vtype & GenePairs$DupType1 == atype),1])
		HyperTests <- rbind(HyperTests, unlist(HypergeometricTest(onum, anum, vnum, length(GenePairs[,1]), atype, vtype, PQvalThreshold)))
	}
}
colnames(HyperTests) <- c("VertebType", "BlanType", "FoldChange", "Observed", "Expected", "HpvalDepleted", "HpvalEnriched")
HyperTests$HqvalDepleted <- qvalue(as.numeric(HyperTests$HpvalDepleted))$qvalues
HyperTests$HqvalEnriched <- qvalue(as.numeric(HyperTests$HpvalEnriched))$qvalues
HyperTests$HResult <- rep("NA", length(HyperTests[,1]))
HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)] <- rep("D", length(HyperTests$HResult[which(HyperTests$HqvalDepleted <= PQvalThreshold)]))
HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)] <- rep("E", length(HyperTests$HResult[which(HyperTests$HqvalEnriched <= PQvalThreshold)]))
print(HyperTests)

write.table(HyperTests, file = paste0(ResultsFolder, "/GenePairs_DupTypeCoocurrence_", Species1, "_", Species2, ".txt"), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)




dev.off()













