#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
#source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))
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
#TandemMaxDist <- 10000
Species1 <- "Blan"
Species2 <- "Ggal"

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
print(head(GenePairs))

# Read gene location information Species 1
system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species1)], ".gtf.gz | sed 's/\\(.*\\)\\t.*\\t.*\\t\\(.*\\)\\t\\(.*\\)\\t.*\\t.*\\t.*\\t.*gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1\\t\\2\\t\\3\\t\\4/g' | sort -k4,4 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2 | sort -k2,2 -k3,3V | awk '{if(chr!=$2){count=1;chr=$2} print $0\"\\t\"count; count=count+1}' | sed 's/chr//g'"), intern=T)
GeneCoord1 <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord1) <- c("Gene1","Chr1","Start1","End1","Index1")
GeneCoord1$Midp <- GeneCoord1$Start + (GeneCoord1$End - GeneCoord1$Start)/2
print(head(GeneCoord1))

GenePairs <- merge(GenePairs, GeneCoord1, by="Gene1")
print(head(GenePairs))


quit()





ExistsFileOhn <- file.exists(OhnologsSFile)
ExistsFileCountsS1S2 <- file.exists(paste0(ResultsFolder, "/CountsAV_", Species1, "_", Species2, ".txt"))
ExistsFileChrLenS1 <- file.exists(paste0(ResultsFolder, "/ChrLen_", Species1, ".txt"))
ExistsFileGeneDataS1 <- file.exists(paste0(ResultsFolder, "/GeneData_", Species1, ".txt"))
ExistsFileChrLenS2 <- file.exists(paste0(ResultsFolder, "/ChrLen_", Species2, ".txt"))
ExistsFileGeneDataS2 <- file.exists(paste0(ResultsFolder, "/GeneData_", Species2, ".txt"))

if(!(ExistsFileOhn & ExistsFileCountsS1S2 & ExistsFileChrLenS1 & ExistsFileGeneDataS1 & ExistsFileChrLenS2 & ExistsFileGeneDataS2)){
	ProduceFilesS1S2(Species1, Species2, OhnologsSFile, ResultsFolder, CountsFileAV)
}
# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]
# Orthologous groups information
CountsAV <- read.delim(paste0(ResultsFolder, "/CountsAV_", Species1, "_", Species2, ".txt"), h=T, stringsAsFactors=F)
# Species 1 Chr length information
ChrLen1 <- read.delim(paste0(ResultsFolder, "/ChrLen_", Species1, ".txt"), h=T, stringsAsFactors=F)
# Species 1 Per gene data 
GeneData1 <- read.delim(paste0(ResultsFolder, "/GeneData_", Species1, ".txt"), h=T, stringsAsFactors=F)
# Species 2 Chr length information
ChrLen2 <- read.delim(paste0(ResultsFolder, "/ChrLen_", Species2, ".txt"), h=T, stringsAsFactors=F)
# Species 2 Per gene data 
GeneData2 <- read.delim(paste0(ResultsFolder, "/GeneData_", Species2, ".txt"), h=T, stringsAsFactors=F)



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

SyntenyAlongChrPlot("Midp", "Gene mid point coordinates", GenePairsSC, GenePairsD)
SyntenyAlongChrPlot("Index", "Gene order", GenePairsSC, GenePairsD)

dev.off()













