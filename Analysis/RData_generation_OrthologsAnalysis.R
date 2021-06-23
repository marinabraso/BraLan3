#!/usr/bin/env Rscript



######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
library(qvalue)

######################################################################
# Files & folders

ResultsFolder <- "Plots"
CountsFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
CountsFileAVO <- "Results/FindOrthologs/AmphVertebOutDeut_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAVO <- "Results/FindOrthologs/AmphVertebOutDeut_broccoli/dir_step3/table_OGs_protein_names.txt"
OnhSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"
RNASeqMetadataFile <- "Metadata/Marletaz2018_RNAseq_SRA.txt"
TPMFile <- "Results/GeneExpression/Gene_TPM_kallisto_tximport.tab"
dNdSFolder <- "Results/dNdSBetweenParalogs/Godon_M8"
BlanGTF <- "Data/Transcriptomes/Branchiostoma_lanceolatum.BraLan3.gtf.gz"
MainGOMFFile <- "Data/GeneOntology/Human_GenesMainGOtermsMF_ENS.txt"
Blan2HsapFile <- "Results/FindOrthologs/AmphVerteb_broccoli/Blan2Hsap_genes.txt"

######################################################################
# General parameters

Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Blan", "Bflo", "Bbel")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Amphi, Verteb, OutDeut)

###########################################################################
# Read data

# List of onhologs
OnhOG.S <- read.table(OnhSFile, h=F)[,1]

# Orthologous groups counts (amphioxus + vertebrates)
CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAV) <- as.character(unlist(lapply(Header[1,], as.character)))
CountsAV$Sum <- rowSums(CountsAV)
CountsAV$SumVerteb <- rowSums(CountsAV[,Verteb])
CountsAV$SumAmphi <- rowSums(CountsAV[,Amphi])
CountsAV$VertebType <- rep("Missing", length(CountsAV[,1]))
CountsAV$VertebType[which(CountsAV$SumVerteb>0)] <- rep("SingleCopy", length(CountsAV$VertebType[which(CountsAV$SumVerteb>0)]))
CountsAV$VertebType[which(apply(CountsAV[,Verteb], 1, max)>1)] <- rep("Duplicated", length(CountsAV$VertebType[which(apply(CountsAV[,Verteb], 1, max)>1)]))
CountsAV$VertebType[which(rownames(CountsAV) %in% OnhOG.S)] <- rep("Ohnolog", length(CountsAV$VertebType[which(rownames(CountsAV) %in% OnhOG.S)]))
CountsAV$AmphiType <- rep("Missing", length(CountsAV[,1]))
CountsAV$AmphiType[which(CountsAV$SumAmphi>0)] <- rep("SingleCopy", length(CountsAV$AmphiType[which(CountsAV$SumAmphi>0)]))
CountsAV$AmphiType[which(apply(CountsAV[,Amphi], 1, max)>1)] <- rep("Duplicated", length(CountsAV$Amp[which(apply(CountsAV[,Amphi], 1, max)>1)]))
CountsAV$BlanType <- rep("Missing", length(CountsAV[,1]))
CountsAV$BlanType[which(CountsAV$Blan>0)] <- rep("SingleCopy", length(CountsAV$BlanType[which(CountsAV$Blan>0)]))
CountsAV$BlanType[which(CountsAV$Blan>1)] <- rep("Duplicated", length(CountsAV$BlanType[which(CountsAV$Blan>1)]))

# Orthologous groups counts (amphioxus + vertebrates + out deuterosome species)
CountsAVO <- read.table(CountsFileAVO, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFileAVO, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAVO) <- as.character(unlist(lapply(Header[1,], as.character)))
CountsAVO$Sum <- rowSums(CountsAVO)
CountsAVO$SumVerteb <- rowSums(CountsAVO[,Verteb])
CountsAVO$SumAmphi <- rowSums(CountsAVO[,Amphi])
CountsAVO$SumOutDeut <- rowSums(CountsAVO[,OutDeut])
CountsAVO$VertebType <- rep("Missing", length(CountsAVO[,1]))
CountsAVO$VertebType[which(CountsAVO$SumVerteb>0)] <- rep("SingleCopy", length(CountsAVO$VertebType[which(CountsAVO$SumVerteb>0)]))
CountsAVO$VertebType[which(apply(CountsAVO[,Verteb], 1, max)>1)] <- rep("Duplicated", length(CountsAVO$VertebType[which(apply(CountsAVO[,Verteb], 1, max)>1)]))
CountsAVO$VertebType[which(rownames(CountsAVO) %in% OnhOG.S)] <- rep("Ohnolog", length(CountsAVO$VertebType[which(rownames(CountsAVO) %in% OnhOG.S)]))
CountsAVO$AmphiType <- rep("Missing", length(CountsAVO[,1]))
CountsAVO$AmphiType[which(CountsAVO$SumAmphi>0)] <- rep("SingleCopy", length(CountsAVO$AmphiType[which(CountsAVO$SumAmphi>0)]))
CountsAVO$AmphiType[which(apply(CountsAVO[,Amphi], 1, max)>1)] <- rep("Duplicated", length(CountsAVO$Amp[which(apply(CountsAVO[,Amphi], 1, max)>1)]))
CountsAVO$BlanType <- rep("Missing", length(CountsAVO[,1]))
CountsAVO$BlanType[which(CountsAVO$Blan>0)] <- rep("SingleCopy", length(CountsAVO$BlanType[which(CountsAVO$Blan>0)]))
CountsAVO$BlanType[which(CountsAVO$Blan>1)] <- rep("Duplicated", length(CountsAVO$BlanType[which(CountsAVO$Blan>1)]))
CountsAVO$OutDeutType <- rep("Missing", length(CountsAVO[,1]))
CountsAVO$OutDeutType[which(CountsAVO$SumOutDeut>0)] <- rep("SingleCopy", length(CountsAVO$OutDeutType[which(CountsAVO$SumOutDeut>0)]))
CountsAVO$OutDeutType[which(apply(CountsAVO[,OutDeut], 1, max)>1)] <- rep("Duplicated", length(CountsAVO$OutDeutType[which(apply(CountsAVO[,OutDeut], 1, max)>1)]))

# dNdS data
system_out <- system(paste("grep '^Final' ", dNdSFolder, "/OG_*_M8_likelihood_ratio.txt | sed 's/Results\\/dNdSBetweenParalogs\\/Godon_M8\\///g' | sed 's/_M8_likelihood_ratio.txt:Final D=/\t/g'", sep=""), intern=T)
GodonM8.D <- read.table(text=system_out, h=F, sep = "\t", row.names=1)
colnames(GodonM8.D) <- c("D")
GodonM8.D$Pval <- unlist(lapply(GodonM8.D$D, p.value.fromM8D))
GodonM8.D$Qval <- qvalue(GodonM8.D$Pval, pi0.method="bootstrap")$qvalue
pdf(paste(ResultsFolder, "/M8_PvalQvalFDR_dist.pdf", sep=""), width=20, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
hist(GodonM8.D$Pval)
hist(GodonM8.D$Qval)
dev.off()

CountsAV$M8.D <- rep(NA, length(CountsAV[,1]))
CountsAV$M8.Pval <- rep(NA, length(CountsAV[,1]))
CountsAV$M8.Qval <- rep(NA, length(CountsAV[,1]))
CountsAV$M8.D[match(rownames(GodonM8.D), rownames(CountsAV))] <- GodonM8.D$D
CountsAV$M8.Pval[match(rownames(GodonM8.D), rownames(CountsAV))] <- GodonM8.D$Pval
CountsAV$M8.Qval[match(rownames(GodonM8.D), rownames(CountsAV))] <- GodonM8.D$Qval

# RNA-seq Metadata
system_out <- system(paste("cat ", RNASeqMetadataFile," | cut -f1,13,24", sep=""), intern=T)
MetaRNA <- read.table(text=system_out, h=F, sep = "\t")
colnames(MetaRNA) <- c("Sample","Age","Tissue")
MetaRNA$Tissue <- sub(" ", ".", MetaRNA$Tissue)
MetaRNA$Age <- sub(" ", ".", MetaRNA$Age)
MetaRNA$Age <- sub("-", ".", MetaRNA$Age)

# RNA-seq gene expression TPM data
TPM <- read.table(TPMFile, h=T)
TPM <- TPM[,MetaRNA$Sample[which(MetaRNA$Tissue!="egg" & MetaRNA$Sample%in%colnames(TPM))]]
MetaRNA <- MetaRNA[which(MetaRNA$Sample %in% colnames(TPM)),]

# Read gene location information
system_out <- system(paste("zcat ", BlanGTF, " | sed 's/gene_id \"\\([A-Z0-9]\\+\\)\"; transcript_id.*gene_name \"\\([A-Za-z0-9\\.-]\\+\\)\".*/\\1\\t\\2/g' | awk '{if($3==\"CDS\"){print $1\"\\t\"$4\"\\t\"$5\"\\t\"$9}}' | sort -k1,1 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2", sep=""), intern=T)
GeneCoord <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord) <- c("Gene","Chr","Start","End")

Tissues <- unique(MetaRNA$Tissue[which(MetaRNA$Age=="adult")])
Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads")
EmbAges <- unique(MetaRNA$Age[which(MetaRNA$Age!="adult")])
EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre.metamorphic.larvae")

###########################################################################
# Data Processing

# TAU calculation per adult tissue
GeneData <- data.frame("Gene" = rownames(TPM))
for(t in c(1:length(Tissues))){
	samples <- MetaRNA$Sample[which(MetaRNA$Tissue==Tissues[t])]
	if(length(samples) == 1){
		GeneData <- cbind(GeneData, as.numeric(TPM[, samples]))
	}else{
		GeneData <- cbind(GeneData, as.numeric(rowMeans(TPM[, samples])))	
	}
}
for(a in c(1:length(EmbAges))){
	samples <- MetaRNA$Sample[which(MetaRNA$Age==EmbAges[a])]
	if(length(samples) == 1){
		GeneData <- cbind(GeneData, as.numeric(TPM[, samples]))
	}else{
		GeneData <- cbind(GeneData, as.numeric(rowMeans(TPM[, samples])))	
	}
}
colnames(GeneData) <- c("Gene", Tissues, EmbAges)
m.GeneData <- as.matrix(GeneData[,Tissues])
Tau <- apply(log2(m.GeneData+1), 1, compute.tau)
MaxGroup <- apply(log2(m.GeneData+1), 1, max.col)
GeneData$TauTissues <- Tau
GeneData$MaxTissue <- MaxGroup
m.GeneData <- as.matrix(GeneData[,EmbAges])
Tau <- apply(log2(m.GeneData+1), 1, compute.tau)
MaxGroup <- apply(log2(m.GeneData+1), 1, max.col)
GeneData$TauEmbAge <- Tau
GeneData$MaxEmbAge <- MaxGroup

GeneData$MeanAdult <- rowMeans(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age=="adult")]])
GeneData$MaxAdult <- apply(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age=="adult")]],1,max)
GeneData$MeanEmbr <- rowMeans(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age!="adult")]])
GeneData$MaxEmbr <- apply(TPM[GeneData$Gene,MetaRNA$Sample[which(MetaRNA$Age!="adult")]],1,max)

# Ortologous groups to gene information 
system_out <- system(paste("cat ", NamesFileAV, " | tail -n +2 | awk '{for(i=2;i<=NF;i++){print $1\"\t\"$i}}' | grep 'BLAG' | sort -u | awk '{if(a[$2]){a[$2]=\"NA\"}else{a[$2]=$1}}END{for(i in a){print a[i]\"\t\"i}}'", sep=""), intern=T)
OG2gene.AV <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2gene.AV) <- c("OG.AV", "Gene")
system_out <- system(paste("cat ", NamesFileAVO, " | tail -n +2 | awk '{for(i=2;i<=NF;i++){print $1\"\t\"$i}}' | grep 'BLAG' | sort -u | awk '{if(a[$2]){a[$2]=\"NA\"}else{a[$2]=$1}}END{for(i in a){print a[i]\"\t\"i}}'", sep=""), intern=T)
OG2gene.AVO <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2gene.AVO) <- c("OG.AVO", "Gene")

# Join OG2gene data with GeneData
GeneData <- merge(GeneData, OG2gene.AV, by="Gene")
GeneData$BlanType.AV <- CountsAV$BlanType[match(GeneData$OG.AV, rownames(CountsAV))]
GeneData$VertebType.AV <- CountsAV$VertebType[match(GeneData$OG.AV, rownames(CountsAV))]
GeneData$M8.D <- CountsAV$M8.D[match(GeneData$OG.AV, rownames(CountsAV))]
GeneData$M8.Qval <- CountsAV$M8.Qval[match(GeneData$OG.AV, rownames(CountsAV))]
GeneData <- merge(GeneData, OG2gene.AVO, by="Gene")
GeneData$BlanType.AVO <- CountsAVO$BlanType[match(GeneData$OG.AV, rownames(CountsAVO))]
GeneData$VertebType.AVO <- CountsAVO$VertebType[match(GeneData$OG.AV, rownames(CountsAVO))]

# Join GeneCoord data with GeneData
GeneData <- merge(GeneData, GeneCoord, by="Gene")
CountsAV$NumChrs <- unlist(lapply(rownames(CountsAV), function(x){length(unique(GeneData$Chr[which(GeneData$OG.AV==x)]))}))
CountsAV$MaxDist <- unlist(lapply(rownames(CountsAV), CalcMaxDist))

GOMFHuman <- read.table(MainGOMFFile, h=F, sep = "\t")
colnames(GOMFHuman) <- c("HGene", "GO")

Blan2Hsap <- read.table(Blan2HsapFile, h=F, sep = "\t")
colnames(Blan2Hsap) <- c("Gene", "HGene")

GOMF <- merge(GOMFHuman, Blan2Hsap, by="HGene")

######################################################################
######################################################################
# Write files
cat("Write files\n")

write.table(GOMF, file = paste(ResultsFolder, "/GeneOntologyMainMF_FromHsap2BlanProcessed.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(GeneData, file = paste(ResultsFolder, "/GeneDataProcessed.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(CountsAV, file = paste(ResultsFolder, "/OGCountsProcessed.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(MetaRNA, file = paste(ResultsFolder, "/RNASeqMetadataProcessed.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)





