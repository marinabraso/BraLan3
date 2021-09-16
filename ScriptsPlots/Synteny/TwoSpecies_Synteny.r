#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))

######################################################################
# Files & folders

ResultsFolder <- "Plots/Synteny"
CountsFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
GTFfolder <- "Results/FilteringGeneSets/GTFs"
ChrLenFolder <- "Results/AssemblyStatistics"

######################################################################
# General parameters
Margin <- 100000000
Species1 <- "Blan"
Species2 <- "Ggal"
HoxGenesList <- c("BLAG12000123", "BLAG12000124", "BLAG12000125", "BLAG12000129", "BLAG12000130", "BLAG12000132", "BLAG12000134", "BLAG12000135", "BLAG12000136", "BLAG12000137", "BLAG12000138", "BLAG12000139", "BLAG12000140", "BLAG12000141", "BLAG12000142")

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan", "Bflo")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum", "Branchiostoma_floridae")
SpeciesGRefs <- c("Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Homo_sapiens.GRCh38", "Mus_musculus.GRCm39", "Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK")

print(paste0("Species1: ", Species1, " ", SpeciesNames[which(ShortSpeciesNames==Species1)], " ", SpeciesGRefs[which(ShortSpeciesNames==Species1)]))
print(paste0("Species2: ", Species2, " ", SpeciesNames[which(ShortSpeciesNames==Species2)], " ", SpeciesGRefs[which(ShortSpeciesNames==Species2)]))

###########################################################################
# Read data

# Read gene pairs and the gene copy number in each species
head <- unlist(strsplit(system(paste0("head -1 ", NamesFileAV, " | sed 's/#//g'"), intern=T), "\t"))
columnS1 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species1)],".fa"))
columnS2 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species2)],".fa"))
system_out <- system(paste0("tail -n +2 ", NamesFileAV, " | awk -F '\t' '{print $1\"\t\"$", columnS1, "\"\t\"$", columnS2,"}' | awk -F '\t' '{split($2,a,\" \");split($3,b,\" \"); for(i in a){for(j in b){print $1\"\t\"a[i]\"\t\"b[j]\"\t\"length(a)\"\t\"length(b)}}}'"), intern=T)
GenePairs <- read.table(text=system_out, h=F, sep = "\t")
colnames(GenePairs) <- c("OG","Gene1","Gene2","CopyNumber1","CopyNumber2")
GenePairs$Type1 <- c("SC", "D")[as.numeric(GenePairs$CopyNumber1>1)+1]
GenePairs$Type2 <- c("SC", "D")[as.numeric(GenePairs$CopyNumber2>1)+1]

# Read gene location information Species 1
system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], ".gtf.gz | sed 's/\\(.*\\)\\t.*\\t.*\\t\\(.*\\)\\t\\(.*\\)\\t.*\\t.*\\t.*\\t.*gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1\\t\\2\\t\\3\\t\\4/g' | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sed 's/chr//g'"), intern=T)
GeneCoord1 <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord1) <- c("Gene1","Chr1","Start1","End1")
GeneCoord1$Midp1 <- GeneCoord1$Start + (GeneCoord1$End - GeneCoord1$Start)/2
# Read gene location information Species 2
system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], ".gtf.gz | sed 's/\\(.*\\)\\t.*\\t.*\\t\\(.*\\)\\t\\(.*\\)\\t.*\\t.*\\t.*\\t.*gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1\\t\\2\\t\\3\\t\\4/g' | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sed 's/chr//g'"), intern=T)
GeneCoord2 <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord2) <- c("Gene2","Chr2","Start2","End2")
GeneCoord2$Midp2 <- GeneCoord2$Start + (GeneCoord2$End - GeneCoord2$Start)/2

# Merge dataframes
GenePairs <- merge(GenePairs, GeneCoord1, by="Gene1")
GenePairs <- merge(GenePairs, GeneCoord2, by="Gene2")

# Chromosome length species 1
system_out <- system(paste0("cat ", ChrLenFolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], "_lengths.txt | sed 's/>//g' | sed 's/chr//g' | grep -P '^[0-9XY]+'"), intern=T)
ChrLen1 <- read.table(text=system_out, h=F, sep = "\t")
colnames(ChrLen1) <- c("Chr","Length","Num","CummLength")
numericchr <- as.numeric(ChrLen1$Chr[grep("[0-9]", ChrLen1$Chr)])
XYchr <- ChrLen1$Chr[grep("[0-9]", ChrLen1$Chr, invert=TRUE)]
ChrLen1 <- ChrLen1[match(c(numericchr[order(numericchr)],XYchr), ChrLen1$Chr),]
ChrLen1$CummLength <- cumsum(as.numeric(ChrLen1$Length))
ChrLen1$CummLengthMargin <- cumsum(as.numeric(ChrLen1$Length)+Margin)

# Chromosome length species 2
system_out <- system(paste0("cat ", ChrLenFolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "_lengths.txt | sed 's/>//g' | sed 's/chr//g' | grep -P '^[0-9XY]+'"), intern=T)
ChrLen2 <- read.table(text=system_out, h=F, sep = "\t")
colnames(ChrLen2) <- c("Chr","Length","Num","CummLength")
numericchr <- as.numeric(ChrLen2$Chr[grep("[0-9]", ChrLen2$Chr)])
XYchr <- ChrLen2$Chr[grep("[0-9]", ChrLen2$Chr, invert=TRUE)]
ChrLen2 <- ChrLen2[match(c(numericchr[order(numericchr)],XYchr), ChrLen2$Chr),]
ChrLen2$CummLength <- cumsum(as.numeric(ChrLen2$Length))
ChrLen2$CummLengthMargin <- cumsum(as.numeric(ChrLen2$Length)+Margin)


GenePairs$CummMidp1 <- GenePairs$Midp1 + ChrLen1$CummLength[match(GenePairs$Chr1, ChrLen1$Chr)] - ChrLen1$Length[match(GenePairs$Chr1, ChrLen1$Chr)]
GenePairs$CummMidp2 <- GenePairs$Midp2 + ChrLen2$CummLength[match(GenePairs$Chr2, ChrLen2$Chr)] - ChrLen2$Length[match(GenePairs$Chr2, ChrLen2$Chr)]
GenePairs$CummMidp1Margin <- GenePairs$Midp1 + ChrLen1$CummLengthMargin[match(GenePairs$Chr1, ChrLen1$Chr)] - ChrLen1$Length[match(GenePairs$Chr1, ChrLen1$Chr)]-Margin
GenePairs$CummMidp2Margin <- GenePairs$Midp2 + ChrLen2$CummLengthMargin[match(GenePairs$Chr2, ChrLen2$Chr)] - ChrLen2$Length[match(GenePairs$Chr2, ChrLen2$Chr)]-Margin

# Printing tables
table(GenePairs[,c("Type1", "Type2")])

pdf(paste0(ResultsFolder, "/", Species1, Species2, "_Synteny.pdf"), width=20, height=20)
par(mar=c(10,10,3,3),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1),nrow=1,ncol=1,byrow=T), widths=c(1), heights=c(1), TRUE)

SyntenyScatterPlot(SpeciesNames[which(ShortSpeciesNames==Species1)], SpeciesNames[which(ShortSpeciesNames==Species2)], ChrLen1, ChrLen2, GenePairs[which(GenePairs$Type1=="SC" & GenePairs$Type2=="SC"),], GenePairs[which(GenePairs$Type1!="SC" | GenePairs$Type2!="SC"),], HoxGenesList)
#layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T), widths=c(3), heights=c(1), TRUE)
#SyntenyLinearPlot(SpeciesNames[which(ShortSpeciesNames==Species1)], SpeciesNames[which(ShortSpeciesNames==Species2)], ChrLen1, ChrLen2, GenePairs[which(GenePairs$Type1=="SC" & GenePairs$Type2=="SC"),], GenePairs[which(GenePairs$Type1!="SC" | GenePairs$Type2!="SC"),])

dev.off()













