#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))
# source("ScriptsPlots/Synteny/TwoSpecies_Synteny_functions.r")

######################################################################
# Files & folders

ResultsFolder <- "Plots/Synteny"
NamesFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
GTFfolder <- "Results/FilteringGeneSets/GTFs"
ChrLenFolder <- "Results/AssemblyStatistics"

######################################################################
# General parameters
Margin <- 100000000
Species1 <- "Blan"
Species2 <- "Blan"
HoxGenesList <- c("BLAG12000123", "BLAG12000124", "BLAG12000125", "BLAG12000129", "BLAG12000130", "BLAG12000132", "BLAG12000134", "BLAG12000135", "BLAG12000136", "BLAG12000137", "BLAG12000138", "BLAG12000139", "BLAG12000140", "BLAG12000141", "BLAG12000142")

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan", "Bflo", "Bbel")
SpeciesNames <- c("D. rerio", "G. gallus", "H. sapiens", "M. musculus", "B.lanceolatum", "B. floridae", "B. belcheri")
SpeciesGRefs <- c("Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Homo_sapiens.GRCh38", "Mus_musculus.GRCm39", "Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK", "Branchiostoma_belcheri.Haploidv18h27")

print(paste0("Species1: ", Species1, " ", SpeciesNames[which(ShortSpeciesNames==Species1)], " ", SpeciesGRefs[which(ShortSpeciesNames==Species1)]))
print(paste0("Species2: ", Species2, " ", SpeciesNames[which(ShortSpeciesNames==Species2)], " ", SpeciesGRefs[which(ShortSpeciesNames==Species2)]))

###########################################################################
# Read data

# Read gene pairs and the gene copy number in each species
head <- unlist(strsplit(system(paste0("head -1 ", NamesFile, " | sed 's/#//g'"), intern=T), "\t"))
columnS1 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species1)],".fa"))
columnS2 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species2)],".fa"))
system_out <- system(paste0("tail -n +2 ", NamesFile, " | awk -F '\t' '{print $1\"\t\"$", columnS1, "\"\t\"$", columnS2,"}' | awk -F '\t' '{split($2,a,\" \");split($3,b,\" \"); for(i in a){for(j in b){print $1\"\t\"a[i]\"\t\"b[j]\"\t\"length(a)\"\t\"length(b)}}}'"), intern=T)
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
system_out <- system(paste0("cat ", ChrLenFolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "_lengths.txt | sed 's/>//g' | sed 's/\\.\\S\\+\\t/\\t/g' | awk -F '\\t' '{if($2 >= 1000000){print $0}}'"), intern=T)
if(Species2 !="Bbel"){
	system_out <- system(paste0("cat ", ChrLenFolder, "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "/", SpeciesGRefs[which(ShortSpeciesNames==Species2)], "_lengths.txt | sed 's/>//g' | sed 's/chr//g' | grep -P '^[0-9XY]+'"), intern=T)
}
ChrLen2 <- read.table(text=system_out, h=F, sep = "\t")
colnames(ChrLen2) <- c("Chr","Length","Num","CummLength")
if(Species2 !="Bbel"){
	numericchr <- as.numeric(ChrLen2$Chr[grep("[0-9]", ChrLen2$Chr)])
	XYchr <- ChrLen2$Chr[grep("[0-9]", ChrLen2$Chr, invert=TRUE)]
	ChrLen2 <- ChrLen2[match(c(numericchr[order(numericchr)],XYchr), ChrLen2$Chr),]
}
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
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)

SyntenyScatterPlotColors(
	SpeciesNames[which(ShortSpeciesNames==Species1)], 
	SpeciesNames[which(ShortSpeciesNames==Species2)], 
	ChrLen1, 
	ChrLen2, 
	GenePairs[which(GenePairs$Type1=="SC" & GenePairs$Type2=="SC"),]
)

SyntenyScatterPlotColors(
	SpeciesNames[which(ShortSpeciesNames==Species1)], 
	SpeciesNames[which(ShortSpeciesNames==Species2)], 
	ChrLen1, 
	ChrLen2, 
	GenePairs[which(GenePairs$Type1!="SC" | GenePairs$Type2!="SC"),]
)

dev.off()



CountsFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
Counts <- read.table(CountsFile, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFile, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
colnames(Counts) <- as.character(unlist(lapply(read.table(text=system_out, h=F, sep = "\t")[1,], as.character)))

BlanGTF <- "Data/Transcriptomes/Branchiostoma_lanceolatum.BraLan3.strong.gtf"
system_out <- system(paste("cat ", BlanGTF, " | awk '{if($3 == \"CDS\"){print $0}}' |  sed 's/.*gene_id \"\\([A-Z0-9]\\+\\)\"; tr.*gene_name \"\\(\\S\\+\\)\".*/\\1\\t\\2/g' | sort -u"), intern=T)
GeneToName <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneToName) <- c("Gene", "Name")
GenePairs$Name1 <- GeneToName$Name[match(GenePairs$Gene1, GeneToName$Gene)]

HsGeneDescFile <- "Data/Human_GeneIDNameDescription.txt"
HsGeneDesc <- read.table(HsGeneDescFile, h=F, sep = "\t", row.names=1)
colnames(HsGeneDesc) <- c("Name", "Desc")
GenePairs$HDesc <- HsGeneDesc$Desc[match(GenePairs$Name1, HsGeneDesc$Name)]




bfloCN <- c(101,94,64,49,39,36,36,33,57,33,50,64,38)
blanCN <- c(63,49,49,46,167,67,111,49,44,43,41,39,40)

for(og in c(1:length(Counts[,1]))){
	if(Counts$Bflo[og]>10 & Counts$Blan[og]>10 & Counts$Hsap[og]==1 & Counts$Ggal[og]==1 & Counts$Mmus[og]==1 & Counts$Drer[og]==1){
		OG <- rownames(Counts[og,])
		print(OG)
		write.table(unique(GenePairs[which(GenePairs$OG==OG),c("OG", "CopyNumber1", "CopyNumber2", "Name1", "HDesc")]), 
			file = paste0(ResultsFolder, "/BlanBflo_GeneExpansionsDescription_SCVert.txt"),
			append = TRUE, col.names = FALSE, quote = F, sep="\t")
	}
}


Counts[which(rownames(Counts) %in% GenePairs$OG[grep("Hypp", GenePairs$Name1)]),]



